#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

const int NUM_ELEMENTS = 100;
// Should be at least 25,000
const int MAX_ITER = 25000;
const int64_t MAX_VALUE = 1000000000000; // 10^12
int64_t A[NUM_ELEMENTS];


// for get_vars
int maxnum, maxidx, nextmax, nextmaxindx;


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

void get_max(int array[]) 
{
    
    maxidx = -1;
    nextmaxindx = -1;
    maxnum = 0;
    nextmax = 0;

    for (int i = 0; i < NUM_ELEMENTS; i++) 
    {
        if (array[i] > maxnum) 
        {
            nextmax = maxnum;
            nextmaxindx = maxidx;
            maxnum = array[i];
            maxidx = i;
        } else if (array[i] > nextmax) 
        {
            nextmax = array[i];
            nextmaxindx = i;
        }
    }
}

int KK(int array[]) 
{
    int residue = 0;
    while (1) 
    {
        get_max(array);
        if (nextmax == 0)
        {
            residue = maxnum;
            break;
        }
        residue = abs(maxnum - nextmax);
        array[nextmaxindx] = 0;
        array[maxidx] = residue;
    }
    return residue;
}

int Sarray[100];

void new_rand_array (){
  
  int num;

  for (int i = 0; i < 100; i++)
  {
    int fnum = rand() % 2 ;

    if (fnum > 0.5)
    {
      num = 1;
    }
    else{
      num = -1;
    }

    Sarray[i] = num;
  }
}

int rand_neighbor_array[100];

// NOT SURE IF THIS IS CORRECT
void new_rand_neighbor_array ()
{
  // generate new random neighbor array to original array
  int i = rand() % 100 ;
  int j = rand() % 100 ;

  // to ensure that i != j
  while (i == j)
  {
    j = rand() % 100 ;
  }

  int chance = rand() % 2 ;
  if (chance == 1)
    rand_neighbor_array[i] = -1 * rand_neighbor_array[i];

  chance = rand() % 2 ;
  if (chance == 1)
    rand_neighbor_array[j] = -1 * rand_neighbor_array[j];

}

int standard_residue(int array[])
{
  int residue = 0;
  for (int i = 0; i < 100; i++)
  {
    residue = residue + (A[i] * array[i]);
  }
  return abs(residue);
}

float T (int iter)
{
  return (pow(10, 10) * pow(0.8, floor(iter/300) ));
}

void standard_solution()
{
    
  int initial_Sarray[100];
  
  int val;

  for (int i = 0; i < 100; i++)
  {
    int fnum = rand() % 2 ;

    if (fnum > 0.5)
    {
      val = 1;
    }
    else{
      val = -1;
    }

    initial_Sarray[i] = val;
  }


  new_rand_array();

  int initial_residue = standard_residue(initial_Sarray);

  printf("Standard Regular Residue = %d\n", initial_residue);

  // Repeated random Residue

  for (int i = 0; i < MAX_ITER; i++)
  {
    new_rand_array();
    int new_residue = standard_residue(Sarray);

    if (new_residue < initial_residue)
    {
      initial_residue = new_residue;
    }
  }

  printf("Standard Repeated random Residue (25,000 iter) = %d\n", abs(initial_residue) );

  // Hill climbing Residue

  // reset initial residue
  initial_residue = abs(standard_residue(initial_Sarray));

  int temp_array[100];

  // init rand_neighbor_array and temp_array
  for (int i = 0; i < 100; i++)
  {
    rand_neighbor_array[i] = initial_Sarray[i];
    temp_array[i] = initial_Sarray[i];
  }

  int neighbor_residue;

  for (int i = 0; i < MAX_ITER; i++)
  {
    new_rand_neighbor_array();
    neighbor_residue = standard_residue(rand_neighbor_array);

    if (neighbor_residue < initial_residue)
    {
      initial_residue = neighbor_residue;
      
      for (int j = 0; j < 100; j++)
      {
        temp_array[j] = rand_neighbor_array[j];
      }
      
    } else {
      
      for (int j = 0; j < 100; j++)
      {
        rand_neighbor_array[j] = temp_array[j]; 
      }
    }

  }

  printf("Standard Hill climbing Residue (25,000 iter) = %d\n", initial_residue );

  // Simulated annealing Reside

  // reset initial residue
  initial_residue = standard_residue(initial_Sarray);

  // init rand_neighbor_array and temp_array
  for (int i = 0; i < 100; i++)
  {
    rand_neighbor_array[i] = initial_Sarray[i];
    temp_array[i] = initial_Sarray[i];
  }

  int temp2_residue;
  temp2_residue = standard_residue(temp_array);
  
  for (int i = 0; i < MAX_ITER; i++)
  {
    new_rand_neighbor_array();
    neighbor_residue = standard_residue(rand_neighbor_array);

    if (neighbor_residue < initial_residue)
    {
      initial_residue = neighbor_residue; // s = s'
      
      for (int j = 0; j < 100; j++)
        {
          temp_array[j] = rand_neighbor_array[j];
        }

    } else {
      // THIS PROBABILITY IS 100% INCORRECT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //bool TrueFalse = (rand() % 100) < 50;
      // (exp((neighbor_residue - initial_residue) / i))

       // printf("neighbor = %d\n", neighbor_residue );
       // printf("initial. = %d\n", initial_residue );
       // printf("sub..... = %d\n", neighbor_residue - initial_residue );
       // printf("cool.... = %f\n", T(i) );
       // printf("no exp.... = %f\n", (neighbor_residue - initial_residue) / T(i) );


       //printf("prob = %f\n", (exp(-(neighbor_residue - initial_residue) / T(i))) );
       float val = (exp(-(neighbor_residue - initial_residue) / T(i)));

      if (val == 1)
      {
        initial_residue = neighbor_residue; // s = s'

        for (int j = 0; j < 100; j++)
        {
          temp_array[j] = rand_neighbor_array[j];
        }
      } else {
        for (int j = 0; j < 100; j++)
        {
          rand_neighbor_array[j] = temp_array[j]; 
        }
      }
    }

    if (initial_residue < temp2_residue)
      temp2_residue = initial_residue;    // S'' = S
  }

  // return initial_residue
  printf("Standard Simulated annealing Residue (25,000 iter) = %d\n", temp2_residue );

}

// main function
int main(int argc, char *argv[]) {

   /* initialize random seed: */
  srand (time(NULL));

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

    standard_solution();


    return 0;
}
