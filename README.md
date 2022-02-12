Sobol
=====
Generate Sobol sequences using a modern C++ interface.

The code in this repository is based on [1] and [2] and the corresponding
source code released by the authors [here](https://web.maths.unsw.edu.au/~fkuo/sobol/).
This code is copyright 2008, Frances Y. Kuo and Stephen Joe. I have referenced only
the key bit-twiddling algorithms from their code in mine and otherwise constructed
a new skeleton around the main framework. I wanted an implementation of Sobol sequences
in a streaming fashion - only one allocation is required, when initializing a sequence.

Usage
-----
The main class is `Sobol::Sequence<IntType = uint32_t, FloatType = double>`.
Initialize a sequence using `auto seq = Sequence(dim);` where dim >= 1.
Retrieve the current point using `seq.get_point()` (returns a `unique_ptr`).
The dimension can be accessed using `seq.dimension()`. Retrieving a point does
not advance the sequence; advance the sequence using `seq.advance()` (which also
returns the new point).

Using a 32-bit integer for IntType will only allow sampling 2^32 points (after this
undefined behavior is encountered). To sample more points, use `uint64_t` for the
template parameter. You might wish to use a smaller int type to get faster performance
if you need few points.

`Sobol::CompileTimeSequence<Dim, IntType = uint32_t, FloatType = double>` has the
same interface, except that `get_point` returns a `std::array`, the constructor
has no arguments, and all methods are `constexpr`. This sequence uses direction numbers 
that are generated via compile-time programming and can generate points at compile time.
At run time, its performance was worse than the normal sequence type in my tests.

With the benchmarking library [hayai](https://bruun.co/2012/02/07/easy-cpp-benchmarking)
available, you can compile 'benchmarks.cpp' to run benchmarks for several different use 
cases.

Design choices
--------------
I have made the default to generate sequences using double-precision floats and a
32-bit integer type for direction numbers. In benchmarking I found that 32-bit generation
is more than twice as fast as 64-bit. It is possible to use **any** unsigned integer
type, so you could construct sequences based on 8-bit integers if you need 256 values
or less. Exceeding the limit of the integer type chosen results in undefined behavior.

The main class is `Sobol::Sequence` which initializes direction numbers when it is
constructed; for large dimensions this initialization is fairly heavy (running into
milliseconds), so if you need to construct such a sequence many times, you may wish
to create the initial sequence then copy it for additional instances using
`Sequence::clone` (memory internally is allocated with `unique_ptr` so there is no
copy constructor).
For small dimensions or smaller integer types initialization is insignficant.
Given an initialized
sequence object, iterating over points is fast and allocation-free.

I included the ability to generate the sequence at compile time using the
`Sobol::CompileTimeSequence<Dim>` template. Using this sequence generator is slower
than the main class at run time; however, you could use it to generate a fixed array
of Sobol points at compile time that will then be embedded in the code. This may give
improved performance in some applications.

