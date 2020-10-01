import sys
import math
import random
import argparse
import optimistic_lock as ol
import threading
from multiprocessing import Process, Queue

# Initial Config & Runtime data
config = {
    "file_name" : None,             # Target tsp file name
    "Output_file" : "solution.csv", # Output file name

    "population" : 50,              # How many genes are in the pool?
    "fitness_evaluations" : 10000,  # The total number of fitness calculation (The program will exit after calculating beyond this number)
    "mutation_factor" : 3,          # How many factors of the gene will be mutated? (Per a mutation function called)
    "convergence_factor" : 10,       # How fast the number of randomly generated genes reduce? (Per a generation)
    "select_ranking" : 0.3,         # How many parents will be selected?
    "replace_ranking" : 0.5,        # How many parents will be replaced?
    "stopping_criteria" : 2000,     # How many generations can pass without improving fitness?
    "total_jobs" : 1,               # How many threads will calculate fitness?

    # Runtime Data
    "dimension" : 0,                # How many nodes are there?
    "Coordinates" : None,
    "Distance_table" : None,
    "table_lock" : None,            # Optimistic sequence lock for distance_table
}

class Coordinate():
    def __init__(self, num, x, y):
        self.num = num
        self.x = x
        self.y = y
    
    def distance(self, other, distance_table, lock):
        #return distance_table[self.num-1][other.num-1]
        seq = lock.read_lock()
        t = None
        while True:
            t = distance_table[self.num-1][other.num-1]
            if lock.validate(seq):
                break
            seq = lock.read_lock()
        if t is None:
            a = abs(self.x - other.x) ** 2
            b = abs(self.y - other.y) ** 2
            distance = math.sqrt(a + b)

            write_seq = lock.write_lock()
            distance_table[self.num - 1][other.num - 1] = distance
            distance_table[other.num - 1][self.num - 1] = distance
            lock.unlock(write_seq)

            return distance
        else:
            return t

class Gene():
    def __init__(self, travel_list=None, dimension=-1, config=None):
        self.config = config
        if travel_list is not None:
            self.gene = travel_list
            self.dimension = len(travel_list)
        else:
            if dimension == -1:
                print("No Specified Dimension when Generating a Gene")
                exit(-1)
            self.dimension = dimension
            #self.__random_generate()
            self.__greedy_generate()
        self.calculated_fitness = self.calc_fitness()

    def get_fitness(self):
        return self.calculated_fitness

    def __repr__(self):
        return str(self.calculated_fitness)
    
    def __lt__(self, other):
        return self.calculated_fitness > other.calculated_fitness
        # Low fitness means higher value

    def __random_generate(self):
        self.gene = list(range(1, self.dimension + 1))
        random.seed()
        random.shuffle(self.gene)
        return

    def __greedy_generate(self):
        random.seed()
        start_point = random.randrange(1, self.dimension + 1)
        self.gene = list()
        self.gene.append(start_point)
        gene_len = 1
        distance_table = self.config["Distance_table"]

        current_point = start_point
        visited = [False for i in range(self.dimension + 1)]
        visited[current_point] = True     # Mark start_point as visited

        while gene_len < self.dimension:
            min_index = None
            min_distance = None
            for i in range(1, self.dimension + 1):
                if visited[i]:
                    continue
                if min_index is None:
                    min_index = i
                    min_distance = distance_table[current_point - 1][i - 1]
                    continue
                current_distance = distance_table[current_point - 1][i - 1]
                if current_distance < min_distance:
                    current_index = i
                    min_distance = current_distance

            #if min_index is None:
            #    break
            current_point = min_index
            self.gene.append(current_point)
            visited[current_point] = True
            gene_len += 1

    
    # Mutate the gene itself
    # Should be called after crossover()
    def mutate(self):
        gene_length = self.dimension
        mutation_factor = config["mutation_factor"]
        for i in range(mutation_factor):
            mutate_point1 = random.randrange(0, gene_length)
            mutate_point2 = random.randrange(0, gene_length)
            # Do point swap mutation
            self.gene[mutate_point1] , self.gene[mutate_point2] = self.gene[mutate_point2], self.gene[mutate_point1]
        return

    def crossover(self, other):
        gene_length = self.dimension
        if gene_length != other.dimension:
            print("Crossover: The length of two genes are different")
            exit(-1)
        random.seed()
        crossover_point = random.randrange(0, gene_length - 1)

        base_list = self.gene[0:crossover_point]
        for i in other.gene:
            if i not in base_list:
                base_list.append(i)
        child = Gene(travel_list=base_list, config=self.config)
        if child.dimension != gene_length:
            print("Crossover: The length of the child is different")
            exit(-1)
        return child

    def calc_fitness(self):
        coordinates = self.config["Coordinates"]
        sum = 0
        a = None
        b = None
        for i in self.gene:
            if a is None:
                a = coordinates[i]
                continue
            b = coordinates[i]
            sum = sum + a.distance(b, self.config['Distance_table'], self.config['table_lock'])
            a = b
        return sum

def threaded_create_random_pool(num, p):
    for i in range(num):
        new_gene = Gene(dimension=p.dimension, config=p.config)
        with p.lock:
            p.pool.append(new_gene)

class Pool():
    def __init__(self, config):
        self.config = config
        self.lock = threading.Lock()
        self.population = config["population"]
        self.dimension = config["dimension"]
        self.pool = list()
        self.create_random_pool()

    # Multi-threaded func
    def create_random_pool(self):
        jobs = self.config["total_jobs"]
        if jobs == 1:
            for i in range(self.population):
                new_gene = Gene(dimension=self.dimension, config=self.config)
                print(new_gene.get_fitness())
                self.pool.append(new_gene)
            return
        population_per_job = int(self.population / jobs)
        remain = self.population - population_per_job * jobs
        job_list = [population_per_job for i in range(jobs)]
        job_list[-1] += remain
        threads = list()
        for i in range(jobs):
            t = threading.Thread(target=threaded_create_random_pool, args=(job_list[i], self))
            threads.append(t)
            t.start()
        for i in threads:
            t.join()

    def replace(self, num):
        for i in range(num):
            self.pool.pop(-1)
        return

    # Multi-threaded func
    def fill_random_gene(self, num):
        jobs = self.config["total_jobs"]
        if jobs == 1:
            for i in range(num):
                new_gene = Gene(dimension=self.dimension, config=self.config)
                self.pool.append(new_gene)
            return
        population_per_job = int(num / jobs)
        remain = num - population_per_job * jobs
        job_list = [population_per_job for i in range(jobs)]
        job_list[-1] += remain
        threads = list()
        for i in range(jobs):
            t = threading.Thread(target=threaded_create_random_pool, args=(job_list[i], self))
            threads.append(t)
            t.start()
        for i in threads:
            t.join()

    def select(self, num):
        if num <= 1:
            return None
        selected_parents = list()
        for i in range(num):
            selected_parents.append(self.pool[i])
        return selected_parents
    
    def generation(self):
        generation = 0
        current_best_gene = None
        population = self.config["population"]
        random_num = population
        stopping_counter = self.config["stopping_criteria"]
        remain_fitness_evaluations = self.config["fitness_evaluations"]

        while(stopping_counter > 0):
            print("Generation: ", generation)

            # Sort the genes of the pool high to low fitness
            self.pool.sort(reverse=True)
            if self.pool[0] == current_best_gene:
                stopping_counter -= 1
            else:
                stopping_counter = self.config["stopping_criteria"]
            current_best_gene = self.pool[0]
            print("Current Best Gene: ", current_best_gene)

            random_num -= self.config["convergence_factor"]
            if random_num < 0:
                random_num = 0
            select_num = int((population - random_num) * self.config["select_ranking"])
            replace_num = int((population - random_num) * self.config["replace_ranking"])

            self.replace(int(replace_num + random_num))
            self.fill_random_gene(int(random_num))
            selected_parents = self.select(select_num)
            if selected_parents is not None:
                self.crossover(selected_parents, select_num, replace_num)
            generation += 1

            remain_fitness_evaluations -= int(replace_num + random_num)
            if remain_fitness_evaluations <= 0:
                break
            
            # Check invariant
            if self.population != len(self.pool):
                print("Population changed")
                exit(-1)
        return current_best_gene

    def crossover(self, selected_parents, select_num, replace_num):
        while replace_num > 0:
            random.seed()
            selectA, selectB = 0, 0
            while selectA == selectB:
                selectA = random.randrange(0, select_num)
                selectB = random.randrange(0, select_num)
            child = selected_parents[selectA].crossover(selected_parents[selectB])
            child.mutate()
            self.pool.append(child)
            replace_num -= 1
        if len(self.pool) != self.config["population"]:
            print("Crossover: Pool population has been changed")
            exit(-1)
        return

def tsp(config):
    try:
        coordinates = read_file(config)
    except Exception as e:
        print("Error when reading a file: EOF error")
        exit(-1)
    config["Coordinates"] = coordinates
    build_distance_table(config)
    result = ga(config)
    print("Result: ", result.gene)
    print_file(result, config)
    return

def fill_distance_table(config):
    table = config['Distance_table']
    dimension = config["dimension"]
    coordinates = config["Coordinates"]
    print("Filling Distance table...")
    progress = 1 / (dimension) * 100
    for i in range(1, dimension + 1):
        progress += progress
        print("Progess: ", round(progress, 2), "%", end="")
        a = coordinates[i]
        for j in range(i, dimension + 1):
            b = coordinates[j]
            distance = a.distance(b, table, config['table_lock'])
            table[i - 1][j - 1] = distance
            table[j - 1][i - 1] = distance
        print("\r", end="")
    print("Done")


def build_distance_table(config):
    dimension = config["dimension"]
    print("Building empty table...")
    table = [[None for col in range(dimension)] for row in range(dimension)]
    print("Done")
    config['Distance_table'] = table
    config["table_lock"] = ol.Optimistic_lock()
    fill_distance_table(config)
    return

def ga(config):
    new_pool = Pool(config)
    best_gene = new_pool.generation()
    return best_gene

#################################################################
def read_file(config):
    print("Reading given file...")
    coordinates = dict()
    with open(config["file_name"], "r") as f:
        name = f.readline()         # NAME
        comment = f.readline()      # COMMENT
        type_ = f.readline()        # TYPE
        dimension = f.readline()    # DIMENSION
        edge_weight_type = f.readline() # EDGE_WEIGHT_TYPE
        temp = ""
        while temp != "NODE_COORD_SECTION\n":
            temp = f.readline()            # loop until find NODE_COORD_SECTION

        dimension = dimension.split()
        dimension = dimension[2]
        config["dimension"] = int(dimension)
        for i in range(int(dimension)):
            string = f.readline()
            string = string.split()
            num = int(string[0])
            x = float(string[1])
            y = float(string[2])
            coordinates[num] = Coordinate(num, x, y)
        if("EOF\n" != f.readline()): # EOF
            raise Exception
    print("Done")
    return coordinates

def print_file(result, config):
    with open(config["Output_file"], "w") as f:
        for i in result.gene:
            f.write(str(i) + "\n")
    return

def define_argparse(config):
    parser = argparse.ArgumentParser(description="Travelling Salesman Problem solver")
    parser.add_argument("Target_file", type=str, help="Target file name")
    parser.add_argument("--output_file", type=str, help="Output file name", default=config["Output_file"])
    parser.add_argument('-p', type=int, help="How many genes are in the pool?", dest="population", default=config["population"])
    parser.add_argument('-f', type=int, help="The total number of fitness calculation (The program will exit after calculating beyond this number)", dest="fitness_evaluations", default=config["fitness_evaluations"])
    parser.add_argument('-c', type=int, help="How fast the number of randomly generated genes reduces? (Per a generation)", dest="convergence_factor", default=config["convergence_factor"])
    parser.add_argument('-s', type=float, help="How many parents will be selected? (format=0.XXX)", dest="select_rank", default=config["select_ranking"])
    parser.add_argument('-r', type=float, help="How many parents will be replaced? (format=0.XXX)", dest="replace_rank", default=config["replace_ranking"])
    parser.add_argument('--criteria', type=int, help="How many generations can pass without improving fitness?", dest="stopping_criteria", default=config["stopping_criteria"])
    parser.add_argument('-m', type=int, help="How many factors of the gene will be mutated? (Per a mutation function called)", dest="mutation_factor", default=config["mutation_factor"])
    parser.add_argument('-j', type=int, help="How many threads will calculate fitness?", dest="jobs", default=config["total_jobs"])
    return parser

def parse_options(config, parser):
    args = parser.parse_args()
    config["file_name"] = args.Target_file
    config["Output_file"] = args.output_file
    config["population"] = args.population
    config["fitness_evaluations"] = args.fitness_evaluations
    config["convergence_factor"] = args.convergence_factor
    config["select_ranking"] = args.select_rank
    config["replace_ranking"] = args.replace_rank
    config["stopping_criteria"] = args.stopping_criteria
    config["mutation_factor"] = args.mutation_factor
    config["total_jobs"] = args.jobs
    return config

if __name__ == "__main__":
    parser = define_argparse(config)
    config = parse_options(config, parser)
    tsp(config)