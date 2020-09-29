# Genetic_TSP_Solver
A travelling salesman problem solver with genetic algorithm

### Usage
```
usage: tsp_solver.py [-h] [--output_file OUTPUT_FILE] [-p POPULATION] [-f FITNESS_EVALUATIONS] [-c CONVERGENCE_FACTOR] [-s SELECT_RANK] [-r REPLACE_RANK] [--criteria STOPPING_CRITERIA]
                     [-m MUTATION_FACTOR] [-j JOBS]
                     Target_file

Travelling Salesman Problem solver

positional arguments:
  Target_file           Target file name

optional arguments:
  -h, --help            show this help message and exit
  --output_file OUTPUT_FILE
                        Output file name
  -p POPULATION         How many genes are in the pool?
  -f FITNESS_EVALUATIONS
                        The total number of fitness calculation (The program will exit after calculating beyond this number)
  -c CONVERGENCE_FACTOR
                        How fast the number of randomly generated genes reduces? (Per a generation)
  -s SELECT_RANK        How many parents will be selected? (format=0.XXX)
  -r REPLACE_RANK       How many parents will be replaced? (format=0.XXX)
  --criteria STOPPING_CRITERIA
                        How many generations can pass without improving fitness?
  -m MUTATION_FACTOR    How many factors of the gene will be mutated? (Per a mutation function called)
  -j JOBS               How many threads will calculate fitness?
```