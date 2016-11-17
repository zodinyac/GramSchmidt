#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 3) {
        printf("Usage: generator_diagonal <output.dat> <size>.\n");
        return 1;
    }
    
    for (int i = 0; i < strlen(argv[2]); i++) {
        char ch = argv[2][i];
        if (!isdigit(ch)) {
            cerr << "Non digit symbol '" << ch << "' in string '" << argv[2] << "'." << endl;
            return 2;
        }
    }
    
    int size = strtol(argv[2], NULL, 10);
    if (size == 0) {
        cerr << "Error converting '" << argv[2] << "' to int." << endl;
        return 3;
    }
    
    ofstream output(argv[1], ios::binary);
    if (!output) {
        cerr << "Can't open file '" << argv[1] << "' for write." << endl;
        return 4;
    }
    
    output.write(reinterpret_cast<const char *>(&size), sizeof(size));
    output.write(reinterpret_cast<const char *>(&size), sizeof(size));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            double val = (i == j);
            output.write(reinterpret_cast<const char *>(&val), sizeof(val));
        }
    }
    output.close();
}
