#include "sobol.hpp"

#include <iostream>

int main()
{
    Sobol::SequenceIterator<-1, 1> iterator(1);

    for (int i = 0; i < 10; ++i)
    {
        std::cout << iterator->operator[](0) << '\n';
        iterator++;
    }

    return 0;
}