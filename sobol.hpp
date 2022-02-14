#ifndef SOBOL_SEQUENCE_HEADER_INCLUDED
#define SOBOL_SEQUENCE_HEADER_INCLUDED

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include "constants.hpp"

namespace Sobol
{

namespace detail
{

#ifdef __GNUC__

int count_trailing_ones(long unsigned i) { return __builtin_ctzl(~i); }

#else

int count_trailing_ones(long unsigned i)
{
    int pos = 0;
    while (i & 1)
    {
        i = i >> 1;
        pos += 1;
    }
    return pos;
}

#endif

} // namespace detail

/*
 * This struct is responsible for computing direction numbers given the
 * * dimension
 * * underlying integer type
 * 
 * These are the v_{jk} from Joe & Kuo scaled by 2^nbits.
 *
 * It is capable of performing the computation at compile time (with C++17) OR can be
 * used at run time. At compile time, a snippet like:
 *
 *   constexpr auto direction_numbers = DirectionNumbers<uint64_t>(dim).numbers();
 *
 * will return an array of direction numbers 1-64 as computed by the compiler.
 */
template <typename IntType>
struct DirectionNumbers
{
    static constexpr unsigned nbits = sizeof(IntType) * 8;

#if __cplusplus >= 201703L
    constexpr
#endif
    DirectionNumbers(unsigned dim) : m_numbers{0}
    {
        if (dim == 1)
        {
            for (unsigned i = 0; i < nbits; ++i)
            {
                m_numbers[i] = IntType(1) << (nbits - i - 1);
            }
        }
        else
        {
            const auto degree = degrees[dim - 2];
            const auto a = coeffs[dim - 2];
            for (uint8_t i = 0; i < degree; ++i)
            {
                m_numbers[i] = ms[dim - 2][i] << (nbits - i - 1);
            }
            for (unsigned i = degree; i < nbits; ++i)
            {
                m_numbers[i] = m_numbers[i - degree] ^ (m_numbers[i - degree] >> degree);
                for (uint8_t j = 0; j < degree - 1; ++j)
                {
                    m_numbers[i] ^= (((a >> (degree - 2 - j)) & 1) * m_numbers[i - j - 1]);
                }
            }
        }
    }

    constexpr const std::array<IntType, nbits> &numbers() noexcept { return m_numbers; }

  private:
    std::array<IntType, nbits> m_numbers;
};

namespace detail
{

#if __cplusplus >= 201703L

template <class IntType, int maxdims>
constexpr auto get_direction_numbers_array()
{
    static_assert(maxdims <= sizeof(coeffs) / sizeof(coeffs[0]) + 1);
    std::array<std::array<IntType, DirectionNumbers<IntType>::nbits>, maxdims> direction_numbers{
        0};

    for (unsigned i = 0; i < direction_numbers.size(); ++i)
    {
        direction_numbers[i] = DirectionNumbers<IntType>(i + 1).numbers();
    }

    return direction_numbers;
}

#endif

template <class IntType>
auto get_direction_numbers(unsigned ndims)
{
    assert(ndims <= sizeof(coeffs) / sizeof(coeffs[0]) + 1);

    std::vector<std::array<IntType, DirectionNumbers<IntType>::nbits>> direction_numbers;
    direction_numbers.reserve(ndims);
    for (unsigned i = 0; i < ndims; ++i)
    {
        direction_numbers.emplace_back(DirectionNumbers<IntType>(i + 1).numbers());
    }
    return direction_numbers;
}

} // namespace detail

namespace detail
{

/*
 * This is the unsafe interface to get the next Sobol point.
 * If index is not correct neither will be the sequence.
 * Returns the index for the next point (index + 1).
 */
template <class IntType>
constexpr IntType advance_sequence(
    unsigned ndims, IntType *x, IntType index,
    const std::vector<std::array<IntType, DirectionNumbers<IntType>::nbits>> &dnums)
{
    assert(ndims <= dnums.size());
    for (unsigned j = 0; j < ndims; ++j)
    {
        x[j] = x[j] ^ dnums[j][count_trailing_ones(index - 1)];
    }
    return index + 1;
}

#if __cplusplus >= 201703L

template <class IntType, size_t Dim>
constexpr IntType
advance_sequence(std::array<IntType, Dim> &x, IntType index)
{
    constexpr auto direction_numbers = detail::get_direction_numbers_array<IntType, Dim>();
    for (unsigned i = 0; i < Dim; ++i)
    {
        x[i] = x[i] ^ direction_numbers[i][count_trailing_ones(index - 1)];
    }
    return index + 1;
}

#endif

} // namespace detail

/*
 * Sequence acts as an iterator for a Sobol sequence (although it does not implement
 * a C++ iterator interface). Construct a sequence given its dimension (>= 1).
 * To retrieve the current point, call `get_point(T *dest)` where T is a
 * floating point type. To advance the sequence, call `advance()` or
 * `advance(T *dest)` to advance the sequence and retrieve the next point in
 * a single call.
 * 
 * The template argument is the integer type to use for computations. For
 * performance, 32-bit is the default. For even faster computations, 16 or
 * even 8-bit integers could be used; HOWEVER, if more than 2^nbits points
 * are generated you will encounter undefined behavior, and for large
 * dimensions direction numbers will begin to overflow the integer type for
 * these smaller types. For generating a few points in a few dimensions,
 * 8 or 16-bit computations will be fastest.
 */
template <class IntType = uint32_t>
class Sequence
{
    static_assert(std::is_integral<IntType>::value && std::is_unsigned<IntType>::value,
        "IntType should be an unsigned integer");

  public:
    Sequence(int N)
        : m_x(new IntType[N]),
          m_dnums(detail::get_direction_numbers<IntType>(N)),
          m_index{1}
    {
        std::fill(m_x.get(), m_x.get() + N, 0);
    }

    unsigned dimension() const noexcept { return m_dnums.size(); }

    void advance() noexcept
    {
        m_index =
            detail::advance_sequence(dimension(), m_x.get(), m_index, m_dnums);
    }

    template <class T>
    void advance(T *dest) noexcept
    {
        static_assert(std::is_floating_point<T>::value, "Point type should be floating point");
        advance();
        get_point(dest);
    }

    template <class T>
    void get_point(T* dest) const noexcept
    {
        static_assert(std::is_floating_point<T>::value, "Point type should be floating point");
        std::transform(m_x.get(), m_x.get() + dimension(), dest, 
            [](IntType x) { return x / std::pow(static_cast<T>(2), DirectionNumbers<IntType>::nbits); });
    }

  private:
    std::unique_ptr<IntType[]> m_x;
    std::vector<std::array<IntType, DirectionNumbers<IntType>::nbits>> m_dnums;
    IntType m_index;
};

#if __cplusplus >= 201703L

/*
 * The CompileTimeSequence has the same interface as Sequence, but does
 * not allocate and all methods are constexpr. Direction numbers are
 * computed at compile time as well. The CompileTimeSequence can then be
 * used to generate a list of integration points that will be embedded in
 * machine code, for example.
 * 
 * For small dimensions, it may be *extremely* fast to use the
 * CompileTimeSequence for run-time computations. This is highly dependent
 * on compiler optimizations, however, and I've also observed while benchmarking
 * that sometimes the CompileTimeSequence is slower, despite having direction
 * numbers precomputed.
 */
template <unsigned Dim, class IntType = uint32_t>
class CompileTimeSequence
{
    static_assert(Dim >= 1 && std::is_integral<IntType>::value && std::is_unsigned<IntType>::value);

  public:
    constexpr CompileTimeSequence() noexcept
        : m_x{0}, m_index{1} 
    {
    }

    constexpr void advance() noexcept
    {
        m_index = detail::advance_sequence(m_x, m_index);
    }

    template <class T>
    constexpr void get_point(T *dest) const noexcept
    {
        static_assert(std::is_floating_point<T>::value);
        for (unsigned i = 0; i < Dim; ++i) {
            dest[i] = m_x[i] / std::pow(static_cast<T>(2), DirectionNumbers<IntType>::nbits);
        }
    }
    
    template <class T>
    void advance(T *dest) noexcept
    {
        static_assert(std::is_floating_point<T>::value);
        advance();
        get_point(dest);
    }

  private:
    std::array<IntType, Dim> m_x;
    IntType m_index;
};

#endif

} // namespace Sobol

/********************************************************************************
 * Below this point is just tests, which can be run by compiling "testsuite.cpp"
 * Reading the tests may be instructive if you're trying to learn how the code
 * functions.
 *******************************************************************************/
#ifdef DOCTEST_LIBRARY_INCLUDED

#endif // DOCTEST_LIBRARY_INCLUDED

#endif // SOBOL_SEQUENCE_HEADER_INCLUDED
