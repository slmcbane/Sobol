/*
 * Copyright 2022 Sean McBane
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in the
 * Software without restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
 * Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 * The core algorithm for computing Sobol sequences was implemented referencing the
 * code by Stephen Joe and Frances Kuo accessed at
 * <https://web.maths.unsw.edu.au/~fkuo/sobol/>. The accompanying file,
 * 'new-joe-kuo-6.21201' and the machine-generated 'constants.hpp' contain the
 * direction numbers computed by the same authors. The copyright notice from their
 * work is reproduced here as required:
 *
 * -----------------------------------------------------------------------------
 * Licence pertaining to sobol.cc and the accompanying sets of direction numbers
 * -----------------------------------------------------------------------------
 * Copyright (c) 2008, Frances Y. Kuo and Stephen Joe
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *
 *     * Neither the names of the copyright holders nor the names of the
 *       University of New South Wales and the University of Waikato
 *       and its contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef SOBOL_SEQUENCE_HEADER_INCLUDED
#define SOBOL_SEQUENCE_HEADER_INCLUDED

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <memory>
#include <type_traits>
#include <utility>

#include "constants.hpp"

namespace Sobol
{

namespace detail
{

#ifdef __GNUC__

// Whether this works is compiler-dependent. Replace with the version below if it
// won't compile.
constexpr inline int count_trailing_ones(long unsigned i) { return __builtin_ctzl(~i); }

#else

constexpr inline int count_trailing_ones(long unsigned i)
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

#if __cplusplus >= 201703L

template <class T, int P>
constexpr T two_to_the_p()
{
    static_assert(std::is_floating_point<T>::value, "T should be a floating point type");

    if constexpr (P < 64)
    {
        return static_cast<T>(static_cast<uint64_t>(1) << P);
    }
    else
    {
        return two_to_the_p<T, 63>() * two_to_the_p<T, 64 - P>();
    }
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
        DirectionNumbers(unsigned dim)
        : m_numbers{0}
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
std::shared_ptr<std::array<IntType, DirectionNumbers<IntType>::nbits>>
get_direction_numbers(unsigned ndims)
{
    assert(ndims <= sizeof(coeffs) / sizeof(coeffs[0]) + 1);
    using numbers_type = std::array<IntType, DirectionNumbers<IntType>::nbits>;

    std::shared_ptr<numbers_type> direction_numbers(
        new numbers_type[ndims], std::default_delete<numbers_type[]>());
    for (unsigned i = 0; i < ndims; ++i)
    {
        *(direction_numbers.get() + i) = DirectionNumbers<IntType>(i + 1).numbers();
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
IntType advance_sequence(
    unsigned ndims, IntType *x, IntType index,
    const std::array<IntType, DirectionNumbers<IntType>::nbits> *dnums)
{
    for (unsigned j = 0; j < ndims; ++j)
    {
        x[j] = x[j] ^ dnums[j][count_trailing_ones(index - 1)];
    }
    return index + 1;
}

#if __cplusplus >= 201703L

template <class IntType, size_t Dim>
constexpr IntType advance_sequence(std::array<IntType, Dim> &x, IntType index)
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
    static_assert(
        std::is_integral<IntType>::value && std::is_unsigned<IntType>::value,
        "IntType should be an unsigned integer");

  public:
    Sequence(unsigned N)
        : m_x(new IntType[N]),
          m_dnums(detail::get_direction_numbers<IntType>(N)), m_index{1}, m_dim{N}
    {
        std::fill(m_x.get(), m_x.get() + N, 0);
    }

    Sequence clone()
    {
        Sequence new_sequence;
        new_sequence.m_x.reset(new IntType[m_dim]);
        for (unsigned i = 0; i < m_dim; ++i)
        {
            new_sequence.m_x[i] = m_x[i];
        }
        new_sequence.m_dnums = m_dnums;
        new_sequence.m_index = m_index;
        new_sequence.m_dim = m_dim;
        return new_sequence;
    }

    unsigned dimension() const noexcept { return m_dim; }

    void advance() noexcept
    {
        m_index = detail::advance_sequence(dimension(), m_x.get(), m_index, m_dnums.get());
    }

    template <class T>
    void advance(T *dest) noexcept
    {
        static_assert(std::is_floating_point<T>::value, "Point type should be floating point");
        advance();
        get_point(dest);
    }

    template <class T>
    void get_point(T *dest) const noexcept
    {
        static_assert(std::is_floating_point<T>::value, "Point type should be floating point");
        std::transform(
            m_x.get(), m_x.get() + m_dim, dest,
            [](IntType x)
            { return x / std::pow(static_cast<T>(2), DirectionNumbers<IntType>::nbits); });
    }

  private:
    std::unique_ptr<IntType[]> m_x;
    std::shared_ptr<std::array<IntType, DirectionNumbers<IntType>::nbits>> m_dnums;
    IntType m_index;
    unsigned m_dim;

    Sequence() = default;
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
    static_assert(
        Dim >= 1 && std::is_integral<IntType>::value && std::is_unsigned<IntType>::value);

  public:
    constexpr CompileTimeSequence() noexcept : m_x{0}, m_index{1} {}

    constexpr void advance() noexcept { m_index = detail::advance_sequence(m_x, m_index); }

    template <class T>
    constexpr void get_point(T *dest) const noexcept
    {
        static_assert(std::is_floating_point<T>::value);
        for (unsigned i = 0; i < Dim; ++i)
        {
            dest[i] = m_x[i] / detail::two_to_the_p<T, DirectionNumbers<IntType>::nbits>();
        }
    }

    template <class T>
    constexpr void advance(T *dest) noexcept
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
 * The tests also serve as examples for usage.
 *
 * These tests assume the use of the direction numbers included with this
 * repository and are not correct with different sets of direction numbers.
 *******************************************************************************/
#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Single-variable runtime generation - first 10 numbers")
{
    Sobol::Sequence<> seq32(1);
    Sobol::Sequence<uint8_t> seq8(1);
    Sobol::Sequence<uint16_t> seq16(1);
    Sobol::Sequence<uint64_t> seq64(1);

    double val, val2;

    seq32.get_point(&val);
    REQUIRE(val == 0);
    seq8.get_point(&val2);
    REQUIRE(val2 == val);
    seq16.get_point(&val2);
    REQUIRE(val2 == val);
    seq64.get_point(&val2);
    REQUIRE(val2 == val);

    seq32.advance(&val);
    REQUIRE(val == 0.5);
    seq8.advance(&val2);
    REQUIRE(val2 == val);
    seq16.advance(&val2);
    REQUIRE(val2 == val);
    seq64.advance(&val2);
    REQUIRE(val2 == val);

    seq32.advance(&val);
    REQUIRE(val == 0.75);
    seq8.advance(&val2);
    REQUIRE(val2 == val);
    seq16.advance(&val2);
    REQUIRE(val2 == val);
    seq64.advance(&val2);
    REQUIRE(val2 == val);

    seq32.advance(&val);
    REQUIRE(val == 0.25);
    seq8.advance(&val2);
    REQUIRE(val2 == val);
    seq16.advance(&val2);
    REQUIRE(val2 == val);
    seq64.advance(&val2);
    REQUIRE(val2 == val);

    seq32.advance(&val);
    REQUIRE(val == 0.375);
    seq8.advance(&val2);
    REQUIRE(val2 == val);
    seq16.advance(&val2);
    REQUIRE(val2 == val);
    seq64.advance(&val2);
    REQUIRE(val2 == val);

    seq32.advance(&val);
    REQUIRE(val == 0.875);
    seq8.advance(&val2);
    REQUIRE(val2 == val);
    seq16.advance(&val2);
    REQUIRE(val2 == val);
    seq64.advance(&val2);
    REQUIRE(val2 == val);

    seq32.advance(&val);
    REQUIRE(val == 0.625);
    seq8.advance(&val2);
    REQUIRE(val2 == val);
    seq16.advance(&val2);
    REQUIRE(val2 == val);
    seq64.advance(&val2);
    REQUIRE(val2 == val);

    seq32.advance(&val);
    REQUIRE(val == 0.125);
    seq8.advance(&val2);
    REQUIRE(val2 == val);
    seq16.advance(&val2);
    REQUIRE(val2 == val);
    seq64.advance(&val2);
    REQUIRE(val2 == val);

    seq32.advance(&val);
    REQUIRE(val == 0.1875);
    seq8.advance(&val2);
    REQUIRE(val2 == val);
    seq16.advance(&val2);
    REQUIRE(val2 == val);
    seq64.advance(&val2);
    REQUIRE(val2 == val);

    seq32.advance(&val);
    REQUIRE(val == 0.6875);
    seq8.advance(&val2);
    REQUIRE(val2 == val);
    seq16.advance(&val2);
    REQUIRE(val2 == val);
    seq64.advance(&val2);
    REQUIRE(val2 == val);
} // TEST_CASE

#if __cplusplus >= 201703L

constexpr auto get_first_ten_terms()
{
    Sobol::CompileTimeSequence<2, uint8_t> sequence;

    std::array<std::array<float, 2>, 10> arr{0};
    sequence.get_point(arr[0].data());

    for (int i = 0; i < 9; ++i)
    {
        sequence.advance(arr[i + 1].data());
    }
    return arr;
}

constexpr auto first_ten = get_first_ten_terms();

TEST_CASE("Test compile-time generation of the first few points")
{
    REQUIRE(first_ten[0][0] == 0);
    REQUIRE(first_ten[0][1] == 0);
    REQUIRE(first_ten[1][0] == 0.5);
    REQUIRE(first_ten[1][1] == 0.5);
    REQUIRE(first_ten[2][0] == 0.75);
    REQUIRE(first_ten[2][1] == 0.25);
    REQUIRE(first_ten[3][0] == 0.25);
    REQUIRE(first_ten[3][1] == 0.75);
    REQUIRE(first_ten[4][0] == 0.375);
    REQUIRE(first_ten[4][1] == 0.375);
    REQUIRE(first_ten[5][0] == 0.875);
    REQUIRE(first_ten[5][1] == 0.875);
    REQUIRE(first_ten[6][0] == 0.625);
    REQUIRE(first_ten[6][1] == 0.125);
    REQUIRE(first_ten[7][0] == 0.125);
    REQUIRE(first_ten[7][1] == 0.625);
    REQUIRE(first_ten[8][0] == 0.1875);
    REQUIRE(first_ten[8][1] == 0.3125);
    REQUIRE(first_ten[9][0] == 0.6875);
    REQUIRE(first_ten[9][1] == 0.8125);
}

TEST_CASE("Compare a bunch of points between run-time and compile-time generation")
{
    Sobol::Sequence run_time_sequence(20);
    Sobol::CompileTimeSequence<20> compile_time_sequence;

    for (int i = 0; i < 66666; ++i)
    {
        run_time_sequence.advance();
        compile_time_sequence.advance();
    }

    double rtval[20];
    double ctval[20];
    for (int i = 0; i < 1000; ++i)
    {
        run_time_sequence.advance(rtval);
        compile_time_sequence.advance(ctval);

        for (int j = 0; j < 20; ++j)
        {
            REQUIRE(rtval[j] == ctval[j]);
        }
    }
} // TEST_CASE

#endif // 201703L

TEST_CASE("Test a bunch of random, high-dimensional points")
{
    Sobol::Sequence<uint64_t> sequence(20201);

    unsigned i = 1;
    auto point = std::unique_ptr<double[]>(new double[20201]);

    while (i != 114)
    {
        sequence.advance();
        i += 1;
    }

    sequence.get_point(point.get());
    REQUIRE(point[8490] == 0.7578125);
    REQUIRE(point[20200] == 0.8984375);

    while (i != 589)
    {
        sequence.advance();
        i += 1;
    }

    sequence.get_point(point.get());
    REQUIRE(point[603] == 0.1083984375);

    while (i != 618)
    {
        sequence.advance();
        i += 1;
    }

    sequence.get_point(point.get());
    REQUIRE(point[4376] == 0.4033203125);

    while (i != 944)
    {
        sequence.advance();
        i += 1;
    }

    sequence.get_point(point.get());
    REQUIRE(point[15308] == 0.5009765625);

    while (i != 4592)
    {
        sequence.advance();
        i += 1;
    }

    sequence.get_point(point.get());
    REQUIRE(point[17373] == 0.9517822265625);

    while (i != 7788)
    {
        sequence.advance();
        i += 1;
    }

    sequence.get_point(point.get());
    REQUIRE(point[16732] == 0.7423095703125);
} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED

#endif // SOBOL_SEQUENCE_HEADER_INCLUDED
