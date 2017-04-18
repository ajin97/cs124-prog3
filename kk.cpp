#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

const int MAX_ITER = 25000;
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

// generate random standard solution
int* generate_rand_sol_standard() {
    int* S = (int*) malloc(sizeof(int) * NUM_ELEMENTS);
    for (int i = 0; i < NUM_ELEMENTS; i++) {
        int rand_num = rand() % 2;
        if (rand_num == 0) {
            S[i] = 1;
        } else {
            S[i] = -1;
        }
    }
    return S;
}

// generate random pre-partitioning solution
int* generate_rand_sol_pp() {
    int* P = (int*) malloc(sizeof(int) * NUM_ELEMENTS);
    for (int i = 0; i < NUM_ELEMENTS; i++) {
        P[i] = rand() % NUM_ELEMENTS;
    }
    return P;
}

// generates random neighbor S' of S
int* generate_rand_neighbor_standard(int* S) {
    int* S_prime = (int*) malloc(sizeof(int) * NUM_ELEMENTS);
    for (int k = 0; k < NUM_ELEMENTS; k++) {
        S_prime[k] = S[k];
    }

    int i = rand() % NUM_ELEMENTS;
    int j = rand() % NUM_ELEMENTS;
    while (i == j) {
        j = rand() % NUM_ELEMENTS;
    }

    S_prime[i] *= -1;
    int rand_num = rand() % 2;
    if (rand_num == 0) {
        S_prime[j] *= -1;
    }

    return S_prime;
}

// generates random neighbor P' of P
int* generate_rand_neighbor_pp(int* P) {
    int* P_prime = (int*) malloc(sizeof(int) * NUM_ELEMENTS);
    for (int i = 0; i < NUM_ELEMENTS; i++) {
        P_prime[i] = P[i];
    }
    int rand_i = rand() % NUM_ELEMENTS;
    int rand_j = rand() % NUM_ELEMENTS;
    while (P_prime[rand_i] == rand_j) {
        rand_i = rand() % NUM_ELEMENTS;
        rand_j = rand() % NUM_ELEMENTS;
    }
    P_prime[rand_i] = rand_j;
    return P_prime;
}

// perform KK heuristic algorithm on array arr
int64_t kk_alg(int64_t* arr) {
    int64_t biggest = 0;
    int biggest_idx = -1;
    int64_t second_biggest = 0;
    int second_biggest_idx = -1;
    int flag = 1;

    while (flag) {
        biggest = 0;
        biggest_idx = -1;
        second_biggest = 0;
        second_biggest_idx = -1;

        // find two biggest elements
        for (int i = 0; i < NUM_ELEMENTS; i++) {
            if (arr[i] >= biggest) {
                second_biggest = biggest;
                second_biggest_idx = biggest_idx;
                biggest = arr[i];
                biggest_idx = i;
            } else if (arr[i] > second_biggest) {
                second_biggest = arr[i];
                second_biggest_idx = i;
            }
        }

        arr[biggest_idx] = biggest - second_biggest;
        arr[second_biggest_idx] = 0;

        if (second_biggest == 0) {
            flag = 0;
        }
    }

    return biggest;
}

// calculate residue on standard representation S
int64_t calc_residue_standard(int* S) {
    int64_t residue = 0;
    for (int i = 0; i < NUM_ELEMENTS; i++) {
        residue += A[i] * S[i];
    }
    return llabs(residue);
}

// calculate residue on pre-partition P using KK
int64_t calc_residue_pp(int* P) {
    // derive new sequence A' from A
    int64_t* A_prime = (int64_t*) calloc(NUM_ELEMENTS, sizeof(int64_t));
    for (int i = 0; i < NUM_ELEMENTS; i++) {
        A_prime[P[i]] += A[i];
    }
    // run the KK heuristic algorithm on the result A'
    int64_t res = kk_alg(A_prime);
    free(A_prime);
    return res;
}

int* copy_P(int* P) {
    int* P_copy = (int*) malloc(sizeof(int) * NUM_ELEMENTS);
    for (int i = 0; i < NUM_ELEMENTS; i++) {
        P_copy[i] = P[i];
    }
    return P_copy;
}

/*
// repeatedly generate random solutions to the problem
int64_t repeated_random_pp() {
    int* P = generate_rand_sol_pp();
    for (int i = 0; i < MAX_ITER; i++) {
        int* P_prime = generate_rand_sol_pp();
        if (calc_residue_pp(P_prime) < calc_residue_pp(P)) {
            free(P);
            P = copy_P(P_prime);
        } else {
            free(P_prime);
        }
    }
    int64_t best_residue = calc_residue_pp(P);
    free(P);
    return best_residue;
}
*/

// repeatedly generate random solutions to the problem
int64_t repeated_random_pp() {
    int* P = generate_rand_sol_pp();
    int64_t best_residue = calc_residue_pp(P);
    free(P);
    for (int i = 0; i < MAX_ITER; i++) {
        P = generate_rand_sol_pp();
        int64_t curr_residue = calc_residue_pp(P);
        if (curr_residue < best_residue) {
            best_residue = curr_residue;
        }
        free(P);
    }
    return best_residue;
}

int64_t repeated_random_standard() {
    int* S = generate_rand_sol_standard();
    int64_t best_residue = calc_residue_standard(S);
    free(S);
    for (int i = 0; i < MAX_ITER; i++) {
        S = generate_rand_sol_standard();
        int64_t curr_residue = calc_residue_standard(S);
        if (curr_residue < best_residue) {
            best_residue = curr_residue;
        }
        free(S);
    }
    return best_residue;
}

// hill climbing
int64_t hill_climbing_pp() {
    int* P = generate_rand_sol_pp();
    for (int i = 0; i < MAX_ITER; i++) {
        int* P_prime = generate_rand_neighbor_pp(P);
        if (calc_residue_pp(P_prime) < calc_residue_pp(P)) {
            free(P);
            P = copy_P(P_prime);
        } else {
            free(P_prime);
        }
    }
    int64_t best_residue = calc_residue_pp(P);
    free(P);
    return best_residue;
}

int64_t hill_climbing_standard() {
    int* S = generate_rand_sol_standard();
    for (int i = 0; i < MAX_ITER; i++) {
        int* S_prime = generate_rand_neighbor_standard(S);
        if (calc_residue_standard(S_prime) < calc_residue_standard(S)) {
            free(S);
            S = copy_P(S_prime);
        } else {
            free(S_prime);
        }
    }
    int64_t best_residue = calc_residue_standard(S);
    free(S);
    return best_residue;
}

// get probability for simulated annealing
float get_prob(int64_t res_P_pr, int64_t res_P, int iter) {
    float T_iter = pow(10.0, 10.0) * pow(0.8, floor(iter / 300.0));
    float prob = exp(-1.0 * (res_P_pr - res_P) / T_iter);
    return prob;
}

// simulated annealing
int64_t simulated_annealing_standard() {
    int* S = generate_rand_sol_standard();
    int64_t best_residue = calc_residue_standard(S);
    for (int i = 0; i < MAX_ITER; i++) {
        int* S_prime = generate_rand_neighbor_standard(S);
        int64_t res_S_prime = calc_residue_standard(S_prime);
        int64_t res_S = calc_residue_standard(S);
        if (res_S_prime < res_S) {
            free(S);
            S = copy_P(S_prime);
        } else {
            float prob = get_prob(res_S_prime, res_S, i);
            float random = ((float) rand()) / ((float) RAND_MAX);
            if (random < prob) {
                free(S);
                S = copy_P(S_prime);
            } else {
                free(S_prime);
            }
        }

        int64_t curr_residue = calc_residue_standard(S);
        if (curr_residue < best_residue) {
            best_residue = curr_residue;
        }
    }
    free(S);
    return best_residue;
}

int64_t simulated_annealing_pp() {
    int* P = generate_rand_sol_pp();
    int64_t best_residue = calc_residue_pp(P);
    for (int i = 0; i < MAX_ITER; i++) {
        int* P_prime = generate_rand_neighbor_pp(P);
        int64_t res_P_prime = calc_residue_pp(P_prime);
        int64_t res_P = calc_residue_pp(P);
        if (res_P_prime < res_P) {
            free(P);
            P = copy_P(P_prime);
        } else {
            float prob = get_prob(res_P_prime, res_P, i);
            float random = ((float) rand()) / ((float) RAND_MAX);
            if (random < prob) {
                free(P);
                P = copy_P(P_prime);
            } else {
                free(P_prime);
            }
        }

        int64_t curr_residue = calc_residue_pp(P);
        if (curr_residue < best_residue) {
            best_residue = curr_residue;
        }
    }
    free(P);
    return best_residue;
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

    /*
    if (input_file) {
        read_numbers_from_file(input_file);
    } else {
        int write_file = 1;
        generate_rand_numbers(write_file);
    }
    */

    for (int i = 0; i < 100; i++) {
        generate_rand_numbers(0);
        cout << simulated_annealing_standard() << endl;
    }

    return 0;
}
