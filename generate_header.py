import sys

def process_line(line, degrees, coeffs, dnumbers):
    if not len(line) >= 4:
        print("Input file does not have the expected format!", file = sys.stderr)
        exit(2)
    
    dim = int(line[0])
    if not dim >= 2:
        print("Input file does not have the expected format - dimension < 2", file = sys.stderr)
        exit(2)

    degree = int(line[1])
    if not degree >= 1:
        print("Input file does not have expected format - degree < 1", file = sys.stderr)
        exit(2)
    elif not len(line) == 4 + degree - 1:
        print("Input file does not have expected format - line length is {} but should be {}"
            .format(len(line), 4 + degree - 1),
            file = sys.stderr)
        exit(2)

    degrees.append(degree)
    coeffs.append(int(line[2]))
    dnumbers.append([int(line[i]) for i in range(3, len(line))])

#################################################################################

if not __name__ == "__main__":
    print("Expected to be run as the main function", file = sys.stderr)
    exit(1)

nargs = len(sys.argv)

if not nargs == 2:
    print("Usage: python3 generate_header.py [input file name]")
    exit(0)

s = []
a = []
m_i = []

with open(sys.argv[1], 'r') as infile:
    line = infile.readline().split()
    if not (len(line) == 4 and line[0] == 'd' and line[1] == 's' and line[2] == 'a' and line[3] == 'm_i'):
        print("Input file does not have the expected format!", file = sys.stderr)
        exit(2)

    for line in infile:
        process_line(line.split(), s, a, m_i)

print('''
#ifndef SOBOL_CONSTANTS_HEADER_INCLUDED
#define SOBOL_CONSTANTS_HEADER_INCLUDED

#include <cstdint>

namespace Sobol {
    
constexpr uint8_t degrees[] = {''')

for i in range(0, len(s) - 1):
    print('  {},'.format(s[i]))

print('''{}
}};

constexpr uint32_t coeffs[] = {{'''.format(s[-1]))

for i in range(0, len(a) - 1):
    print('  {},'.format(a[i]))

print('''{}
}};

constexpr uint64_t ms[{}][{}] = {{'''.format(a[-1], len(m_i), max(s)))

for i in range(0, len(m_i) - 1):
    print('  {', end = '')
    for j in range(0, len(m_i[i])-1):
        print('{}, '.format(m_i[i][j]), end = '')
    print('{}}},'.format(m_i[i][-1]))

print('  {', end = '')

for j in range(0, len(m_i[-1])-1):
    print('{}, '.format(m_i[-1][j]), end = '')

print('{}}}'.format(m_i[-1][-1]))

print('''};

} // namespace Sobol

#endif // SOBOL_CONSTANTS_HEADER_INCLUDED''')
