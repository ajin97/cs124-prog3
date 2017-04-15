#include <iostream>
#include <fstream>

using namespace std;

const int NUM_ELEMENTS = 100;
const int64_t MAX_VALUE = 1000000000000; // 10^12
int64_t A[NUM_ELEMENTS];

// returns random integer in [1, 10^12]
int64_t rand_number() {
    int64_t low = rand();
    int64_t high = rand();

    int64_t rand_val = low | (high << 32);
    return ((rand_val % MAX_VALUE) + 1);
}

// fills A with NUM_ELEMENTS random integers in [1, 10^12]
void generate_rand_numbers(int write_file) {
    for (int i = 0; i < NUM_ELEMENTS; i++) {
        A[i] = rand_number();
    }

    if (write_file) {
        ofstream fout ("input.txt");
        // write array A to file
        for (int i = 0; i < NUM_ELEMENTS; i++) {
            fout << A[i] << endl;
        }
    }
}

// read NUM_ELEMENTS integers from the given file into array A
void read_numbers_from_file(char* file) {
    ifstream fin (file);
    // read numbers from file
    for (int i = 0; i < NUM_ELEMENTS; i++) {
        fin >> A[i];
    }
}

// main function
int main(int argc, char *argv[]) {
    srand(time(0));
    if (argc != 1 && argc != 2) {
        cout << "Usage: ./kk inputfile" << endl;
        return 0;
    }

    char* input_file = NULL;
    if (argc == 2) {
        input_file = argv[1];
    }

    if (input_file) {
        read_numbers_from_file(input_file);
    } else {
        int write_file = 1;
        generate_rand_numbers(write_file);
    }

    return 0;
}
