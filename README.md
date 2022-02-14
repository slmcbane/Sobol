Sobol
=====
Generate Sobol sequences using a modern C++ interface.

The code in this repository is based on [1] and [2] and the corresponding
source code released by the authors [here](https://web.maths.unsw.edu.au/~fkuo/sobol/).
This code is copyright 2008, Frances Y. Kuo and Stephen Joe; the full copyright notice
is reproduced in 'sobol.hpp' and 'constants.hpp'. I have referenced only
the key bit-twiddling algorithms from their code in mine and otherwise constructed
a new skeleton around the main algorithm. I wanted an implementation of Sobol sequences
in a streaming fashion - only one allocation is required, when initializing a sequence -
and the possibility of computing direction numbers at compile time, or using a different
integer type.

Usage
-----
The main class is `Sobol::Sequence<IntType = uint32_t, FloatType = double>`.
Initialize a sequence using `Sequence(dim);` where dim >= 1.
Retrieve the current point using `seq.get_point(dest)`, which writes the point
into the pointer `dest` (make sure this does not write out of bounds!).
The dimension can be accessed using `seq.dimension()`. Retrieving a point does
not advance the sequence; advance the sequence using `seq.advance()`. To advance
the sequence and retrieve the next point in one operation, use `seq.advance(dest)`.

Using a 32-bit integer for IntType will only allow sampling 2^32 points (after this
undefined behavior is encountered). To sample more points, use `uint64_t` for the
template parameter. You might wish to use a smaller int type to get faster performance
if you need few points, this will only work for a few dimensions before direction
numbers start to overflow the type.

If compiling with C++17 support,
`Sobol::CompileTimeSequence<Dim, IntType = uint32_t, FloatType = double>` has the
same interface, but works at compile time. This sequence uses direction numbers 
that are generated at compile time and can generate points at compile time. Run time
performance is highly dependent on compiler optimizations, but may be extremely fast
for low-dimensional sequences (due to precomputed direction numbers embedded directly
in the sequence generation). I did not observe a consistent trend in benchmarks - it
depends on your compiler and what flags you give it.

With the benchmarking library [hayai](https://bruun.co/2012/02/07/easy-cpp-benchmarking)
available, you can compile 'benchmarks.cpp' to run benchmarks for several different use 
cases.

As a final note, initializing the sequence for large dimensions is relatively slow
since it has to compute the direction numbers. "Slow" is relative - for thousands of
dimensions the initialization may take > 1 ms. If you will be creating new
high-dimensional sequences often, you may wish to create the starting sequence object
one time and use `Sequence::clone()` to create a copy. There is no copy constructor
due to the use of `std::unique_ptr` internally; I wanted it to be clear that copying
the sequence will duplicate allocations. The direction numbers will be shared and not
allocated twice.

Copyright
---------
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in the
Software without restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
