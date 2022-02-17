#include "sobol.hpp"

#include <cstdio>

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        fprintf(
            stderr, "Expects 2 command-line arguments: number of points to generate, and "
                    "dimensionality\n");
        exit(1);
    }

    int n = atoi(argv[1]);
    int dim = atoi(argv[2]);

    Sobol::RandomizedSequence<uint64_t> seq(dim);
    std::unique_ptr<double[]> point(new double[dim]);

    for (int i = 0; i < n; ++i)
    {
        seq.get_point(point.get());
        seq.advance();

        for (int j = 0; j < dim - 1; ++j)
        {
            printf("%1.14E ", point[j]);
        }
        printf("%1.14E\n", point[dim - 1]);
    }

    return 0;
}

