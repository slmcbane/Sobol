#include "sobol.hpp"
#include "hayai.hpp"

BENCHMARK(RunTime, SetUpDimension10, 1000, 1000)
{
    auto seq = Sobol::Sequence<>(10);
}

BENCHMARK(RunTime, SetUpDimension25, 1000, 500)
{
    auto seq = Sobol::Sequence<>(25);
}

BENCHMARK(RunTime, SetUpDimension200, 1000, 100)
{
    auto seq = Sobol::Sequence<>(200);
}

BENCHMARK(RunTime, SetUpDimension10_64bit, 1000, 1000)
{
    auto seq = Sobol::Sequence<uint64_t>(10);
}

BENCHMARK(RunTime, SetUpDimension25_64bit, 1000, 500)
{
    auto seq = Sobol::Sequence<uint64_t>(25);
}

BENCHMARK(RunTime, SetUpDimension200_64bit, 1000, 100)
{
    auto seq = Sobol::Sequence<uint64_t>(200);
}

class SequenceFixture10 : public hayai::Fixture
{
    public:
    virtual void SetUp()
    {
        sequence.reset(new Sobol::Sequence<>(10));
    }

    virtual void TearDown()
    {
    }

    std::unique_ptr<Sobol::Sequence<>> sequence;
};

class SequenceFixture25 : public hayai::Fixture
{
    public:
    virtual void SetUp()
    {
        sequence.reset(new Sobol::Sequence<>(25));
    }

    virtual void TearDown() {} 
    std::unique_ptr<Sobol::Sequence<>> sequence;
};

class SequenceFixture200 : public hayai::Fixture
{
    public:
    virtual void SetUp()
    {
        sequence.reset(new Sobol::Sequence<>(200));
    }

    virtual void TearDown() {} 

    std::unique_ptr<Sobol::Sequence<>> sequence;
};

BENCHMARK_F(SequenceFixture10, SampleDimension10, 1000, 100)
{
    double point[10];
    for (int i = 0; i < 1000; ++i)
    {
        sequence->advance(point);
    }
}

BENCHMARK_F(SequenceFixture25, SampleDimension25, 1000, 100)
{
    double point[25];
    for (int i = 0; i < 1000; ++i)
    {
        sequence->advance(point);
    }
}

BENCHMARK_F(SequenceFixture200, SampleDimension200, 1000, 20)
{
    double point[200];
    for (int i = 0; i < 1000; ++i)
    {
        sequence->advance(point);
    }
}

class SequenceFixture10_64bit : public hayai::Fixture
{
    public:
    virtual void SetUp()
    {
        sequence.reset(new Sobol::Sequence<uint64_t>(10));
    }

    virtual void TearDown()
    {
    }

    std::unique_ptr<Sobol::Sequence<uint64_t>> sequence;
};

class SequenceFixture25_64bit : public hayai::Fixture
{
    public:
    virtual void SetUp()
    {
        sequence.reset(new Sobol::Sequence<uint64_t>(25));
    }

    virtual void TearDown() {} 

    std::unique_ptr<Sobol::Sequence<uint64_t>> sequence;
};

class SequenceFixture200_64bit : public hayai::Fixture
{
    public:
    virtual void SetUp()
    {
        sequence.reset(new Sobol::Sequence<uint64_t>(200));
    }

    virtual void TearDown() {} 

    std::unique_ptr<Sobol::Sequence<uint64_t>> sequence;
};

BENCHMARK_F(SequenceFixture10_64bit, SampleDimension10, 1000, 100)
{
    double point[10];
    for (int i = 0; i < 1000; ++i)
    {
        sequence->advance(point);
    }
}

BENCHMARK_F(SequenceFixture25_64bit, SampleDimension25, 1000, 100)
{
    double point[25];
    for (int i = 0; i < 1000; ++i)
    {
        sequence->advance(point);
    }
}

BENCHMARK_F(SequenceFixture200_64bit, SampleDimension200, 1000, 20)
{
    double point[200];
    for (int i = 0; i < 1000; ++i)
    {
        sequence->advance(point);
    }
}

#if __cplusplus >= 201703L

BENCHMARK(CompileTime, SampleDimension10, 1000, 20)
{
    Sobol::CompileTimeSequence<10> sequence;
    double point[10];

    for (int i = 0; i < 1000; ++i)
    {
        sequence.advance(point);
    }
}
BENCHMARK(CompileTime, SampleDimension25, 1000, 20)
{
    Sobol::CompileTimeSequence<25> sequence;
    double point[25];

    for (int i = 0; i < 1000; ++i)
    {
        sequence.advance(point);
    }
}

#endif

int main(int argc, char *argv[])
{
    hayai::ConsoleOutputter consoleOutputter;

    hayai::Benchmarker::AddOutputter(consoleOutputter);
    hayai::Benchmarker::RunAllTests();

    return 0;
}