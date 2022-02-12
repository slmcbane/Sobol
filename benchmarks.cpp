#include "sobol.hpp"
#include "hayai.hpp"

#include <optional>

BENCHMARK(RunTime, SetUpDimension10, 1000, 1000)
{
    const auto &arr = Sobol::Sequence(10).get_point();
}

BENCHMARK(RunTime, SetUpDimension25, 1000, 500)
{
    const auto &arr = Sobol::Sequence(25).get_point();
}

BENCHMARK(RunTime, SetUpDimension200, 1000, 100)
{
    const auto &arr = Sobol::Sequence(200).get_point();
}

BENCHMARK(RunTime, SetUpDimension10_64bit, 1000, 1000)
{
    const auto &arr = Sobol::Sequence<uint64_t>(10).get_point();
}

BENCHMARK(RunTime, SetUpDimension25_64bit, 1000, 500)
{
    const auto &arr = Sobol::Sequence<uint64_t>(25).get_point();
}

BENCHMARK(RunTime, SetUpDimension200_64bit, 1000, 100)
{
    const auto &arr = Sobol::Sequence<uint64_t>(200).get_point();
}

class SequenceFixture10 : public hayai::Fixture
{
    public:
    virtual void SetUp()
    {
        sequence.emplace(10);
    }

    virtual void TearDown()
    {
    }

    std::optional<Sobol::Sequence<uint32_t>> sequence;
};

class SequenceFixture25 : public hayai::Fixture
{
    public:
    virtual void SetUp()
    {
        sequence.emplace(25);
    }

    virtual void TearDown() {} 

    std::optional<Sobol::Sequence<uint32_t>> sequence;
};

class SequenceFixture200 : public hayai::Fixture
{
    public:
    virtual void SetUp()
    {
        sequence.emplace(200);
    }

    virtual void TearDown() {} 

    std::optional<Sobol::Sequence<uint32_t>> sequence;
};

BENCHMARK_F(SequenceFixture10, SampleDimension10, 1000, 100)
{
    for (int i = 0; i < 1000; ++i)
    {
        const auto &point = sequence->advance();
    }
}

BENCHMARK_F(SequenceFixture25, SampleDimension25, 1000, 100)
{
    for (int i = 0; i < 1000; ++i)
    {
        const auto &point = sequence->advance();
    }
}

BENCHMARK_F(SequenceFixture200, SampleDimension200, 1000, 20)
{
    for (int i = 0; i < 1000; ++i)
    {
        const auto &point = sequence->advance();
    }
}

class SequenceFixture10_64bit : public hayai::Fixture
{
    public:
    virtual void SetUp()
    {
        sequence.emplace(10);
    }

    virtual void TearDown()
    {
    }

    std::optional<Sobol::Sequence<uint64_t>> sequence;
};

class SequenceFixture25_64bit : public hayai::Fixture
{
    public:
    virtual void SetUp()
    {
        sequence.emplace(25);
    }

    virtual void TearDown() {} 

    std::optional<Sobol::Sequence<uint64_t>> sequence;
};

class SequenceFixture200_64bit : public hayai::Fixture
{
    public:
    virtual void SetUp()
    {
        sequence.emplace(200);
    }

    virtual void TearDown() {} 

    std::optional<Sobol::Sequence<uint64_t>> sequence;
};

BENCHMARK_F(SequenceFixture10_64bit, SampleDimension10, 1000, 100)
{
    for (int i = 0; i < 1000; ++i)
    {
        const auto &point = sequence->advance();
    }
}

BENCHMARK_F(SequenceFixture25_64bit, SampleDimension25, 1000, 100)
{
    for (int i = 0; i < 1000; ++i)
    {
        const auto &point = sequence->advance();
    }
}

BENCHMARK_F(SequenceFixture200_64bit, SampleDimension200, 1000, 20)
{
    for (int i = 0; i < 1000; ++i)
    {
        const auto &point = sequence->advance();
    }
}

BENCHMARK(CompileTime, SampleDimension10, 1000, 20)
{
    Sobol::CompileTimeSequence<10> sequence;

    for (int i = 0; i < 1000; ++i)
    {
        const auto &point = sequence.advance();
    }
}
BENCHMARK(CompileTime, SampleDimension25, 1000, 20)
{
    Sobol::CompileTimeSequence<25> sequence;

    for (int i = 0; i < 1000; ++i)
    {
        const auto &point = sequence.advance();
    }
}

int main(int argc, char *argv[])
{
    hayai::ConsoleOutputter consoleOutputter;

    hayai::Benchmarker::AddOutputter(consoleOutputter);
    hayai::Benchmarker::RunAllTests();

    return 0;
}