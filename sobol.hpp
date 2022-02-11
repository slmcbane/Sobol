#ifndef SOBOL_SEQUENCE_HEADER_INCLUDED
#define SOBOL_SEQUENCE_HEADER_INCLUDED

#include <type_traits>
#include <utility>

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

} // namespace detail

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
