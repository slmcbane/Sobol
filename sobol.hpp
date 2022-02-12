#ifndef SOBOL_SEQUENCE_HEADER_INCLUDED
#define SOBOL_SEQUENCE_HEADER_INCLUDED

#include <array>
#include <cassert>
#include <cmath>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

#include "constants.hpp"

namespace Sobol
{

namespace detail
{

template <auto Val>
constexpr auto type_constant = std::integral_constant<std::decay_t<decltype(Val)>, Val>{};

template <int Start, class F, int... Is>
void compile_time_for(F &&f, std::integer_sequence<int, Is...>)
{
    (static_cast<void>(f(type_constant<Is>)), ...);
}

template <int Start, int Stop, class F>
void compile_time_for(F &&f)
{
    static_assert(Stop >= Start);
    constexpr auto index_sequence = std::make_integer_sequence<int, Stop - Start>();
    compile_time_for<Start>(std::forward<F>(f), index_sequence);
}

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
 * It is capable of performing the computation at compile time OR can be
 * used at run time. At compile time, a snippet like:
 *
 *   constexpr auto direction_numbers = DirectionNumbers<uint64_t>(dim).numbers();
 *
 * will return an array of direction numbers 1-64 as computed by the compiler.
 *
 * Note that the array returned by numbers() is the direction numbers scaled by
 * 2^nbits, where nbits is the number of bits in IntType.
 */
template <typename IntType>
struct DirectionNumbers
{
    static constexpr unsigned nbits = sizeof(IntType) * 8;

    constexpr DirectionNumbers(unsigned dim) : m_numbers{0}
    {
        if (dim == 1)
        {
            for (int i = 0; i < nbits; ++i)
            {
                m_numbers[i] = IntType(1) << (nbits - i - 1);
            }
        }
        else
        {
            const auto degree = degrees[dim - 2];
            const auto a = coeffs[dim - 2];
            for (int i = 0; i < degree; ++i)
            {
                m_numbers[i] = ms[dim - 2][i] << (nbits - i - 1);
            }
            for (int i = degree; i < nbits; ++i)
            {
                m_numbers[i] = m_numbers[i - degree + 1] ^ (m_numbers[i - degree + 1] >> degree);
                for (int j = 0; j < degree - 1; ++j)
                {
                    m_numbers[i] ^= (((a >> (degree - 2 - j)) & 1) * m_numbers[i - j - 1]);
                }
            }
        }
    }

    constexpr std::array<IntType, nbits> numbers() noexcept { return m_numbers; }

  private:
    std::array<IntType, nbits> m_numbers;
};

namespace detail
{

template <class IntType, int maxdims>
constexpr auto get_direction_numbers_array()
{
    static_assert(maxdims <= sizeof(coeffs) / sizeof(coeffs[0]) + 1);
    std::array<std::array<IntType, DirectionNumbers<IntType>::nbits>, maxdims> direction_numbers{
        0};

    for (int i = 0; i < direction_numbers.size(); ++i)
    {
        direction_numbers[i] = DirectionNumbers<IntType>(i + 1).numbers();
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
template <
    int maxdims, class IntType, class T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
constexpr unsigned long advance_sequence(int ndims, T *point, IntType *x, unsigned long index)
{
    assert(ndims <= maxdims);
    constexpr auto direction_numbers = detail::get_direction_numbers_array<IntType, maxdims>();
    for (int j = 0; j < ndims; ++j)
    {
        x[j] = x[j] ^ direction_numbers[j][count_trailing_ones(index) + 1];
        point[j] = x[j] / std::pow(static_cast<T>(2), 64);
    }
    return index + 1;
}

template <
    int ndims, int maxdims = ndims, class IntType, class T,
    typename = std::enable_if_t<std::is_floating_point_v<T>>>
constexpr unsigned long advance_sequence(T *point, IntType *x, unsigned long index)
{
    static_assert(ndims <= maxdims);
    constexpr auto direction_numbers = detail::get_direction_numbers_array<IntType, maxdims>();
    detail::compile_time_for<0, ndims>([=](auto I) {
        x[I()] = x[I()] ^ direction_numbers[I()][count_trailing_ones(index) + 1];
        point[I()] = point[I()] / std::pow(static_cast<T>(2), 64);
    });
    return index + 1;
}

} // namespace detail

template <
    int ndims = -1, int maxdims = (ndims > 0 ? ndims : sizeof(coeffs) / sizeof(coeffs[0]) + 1),
    class T = double, class IntType = uint64_t>
class SequenceIterator : public std::iterator<std::input_iterator_tag, const std::vector<T>>
{
    static_assert(ndims == -1 || ndims >= 1);

  public:
    template <bool NotDynamic = (ndims > 0), std::enable_if_t<NotDynamic, int> = 0>
    SequenceIterator() : m_index{1}, m_x(ndims, 0), m_point(ndims, 0)
    {
    }

    template <bool Dynamic = (ndims == -1), std::enable_if_t<Dynamic, int> = 0>
    SequenceIterator(int N) : m_index{1}, m_x(N, 0), m_point(N, 0)
    {
        assert(N >= 1 && "Dimension less than 1 doesn't make sense\n");
    }

    bool operator!=(const SequenceIterator &other) const noexcept
    {
        return (m_point.size() != other.m_point.size() || m_index != other.m_index);
    }

    const std::vector<T> &operator*() const noexcept { return m_point; }

    const std::vector<T> *operator->() const noexcept { return &m_point; }

    SequenceIterator &operator++() noexcept
    {
        if constexpr (ndims == -1)
        {
            m_index = detail::advance_sequence<maxdims>(
                m_point.size(), m_point.data(), m_x.data(), m_index);
        }
        else
        {
            m_index =
                detail::advance_sequence<ndims, maxdims>(m_point.data(), m_x.data(), m_index);
        }

        return *this;
    }

    SequenceIterator &operator++(int) noexcept { return ++(*this); }

  private:
    unsigned long m_index;
    std::vector<IntType> m_x;
    std::vector<T> m_point;
};

} // namespace Sobol

/********************************************************************************
 * Below this point is just tests, which can be run by compiling "testsuite.cpp"
 * Reading the tests may be instructive if you're trying to learn how the code
 * functions.
 *******************************************************************************/
#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Test compile_time_for")
{

    int x[10] = {0};

    Sobol::detail::compile_time_for<0, 10>([&](auto I) {
        for (int i = I(); i < 10; ++i)
        {
            x[i] += 1;
        }
    });

    for (int i = 0; i < 10; ++i)
    {
        REQUIRE(x[i] == i + 1);
    }

} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED

#endif // SOBOL_SEQUENCE_HEADER_INCLUDED
