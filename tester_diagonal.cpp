#include <fstream>
#include <iostream>
using namespace std;

template <typename T>
T read_val(ifstream &input)
{
    T val;
    input.read(reinterpret_cast<char *>(&val), sizeof(val));
    return val;
}

int main(int argc, char *argv[])
{
    if (argc != 2) {
        printf("Usage: tester_diagonal <input.dat>.\n");
        return 1;
    }
    
    ifstream input(argv[1], ios::binary);
    int size1 = read_val<int>(input), size2 = read_val<int>(input);
    
    if (size1 != size2) {
        cerr << "Matrix is not diagonal: " << size1 << " x " << size2 << "." << endl;
        return 2;
    }
    
    for (int i = 0; i < size1; i++) {
        for (int j = 0; j < size1; j++) {
            double val = read_val<double>(input);
            
            if (i == j && val != 1) {
                cerr << "Matrix is not diagonal: " << val << " but expected 1 at (" << i << ", " << j << ")." << endl;
                return 3;
            } else if (i != j && val != 0) {
                cerr << "Matrix is not diagonal: " << val << " but expected 0 at (" << i << ", " << j << ")." << endl;
                return 4;
            }
        }
    }
    input.close();
    
    cout << "Matrix is correct." << endl;
}
