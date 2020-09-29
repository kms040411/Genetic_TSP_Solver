import sys
import math
import random
import argparse

# Initial Config
config = {
    "debug" : True,                 # Debug mode
    "file_name" : None,             # Target tsp file name

    "population" : 0,              # How many genes are in the pool?
    "fitness_evalutaions" : 10,
    "mutation_factor" : 5,          # How many factors of the gene will be mutated? (Per a mutation function called)
    "convergence_factor" : 1,       # How fast the number of randomly generated genes reduce? (Per a generation)
    "select_ranking" : 0.5,         # How many parents will be selected?
    "replace_ranking" : 0.2,        # How many parents will be replaced?
    "stopping_criteria" : 10,       # How many generations can pass without improving fitness?

    # Runtime Data
    "dimension" : 0,                # How many nodes are there?
    "Coordinates" : None,
    "Distance_table" : None,
    "current_random_genes" : 0,
}

class Coordinate():
    def __init__(self, num, x, y):
        self.num = num
        self.x = x
        self.y = y
    
    def distance(self, other, distance_table):
        if distance_table[self.num-1][other.num-1] != distance_table[other.num-1][self.num-1]:
            print("Coordinate-distance(): table error")
            exit(-1)
        if distance_table[self.num-1][other.num-1] is None:
            a = abs(self.x - other.x) ** 2
            b = abs(self.y - other.y) ** 2
            distance = math.sqrt(a + b)
            distance_table[self.num - 1][other.num - 1] = distance
            distance_table[other.num - 1][self.num - 1] = distance
            return distance
        else:
            return distance_table[self.num-1][other.num-1]
        

class Gene():
    def __init__(self, travel_list=None, dimension=-1, config=None):
        if travel_list is not None:
            self.gene = travel_list
            self.dimension = len(travel_list)
        else:
            if dimension == -1:
                print("No Specified Dimension when Generating a Gene")
                exit(-1)
            self.dimension = dimension
            self.__random_generate()
        self.config = config
        self.calculated_fitness = self.calc_fitness()

    def get_fitness(self):
        return self.calculated_fitness

    def __repr__(self):
        #return str(self.gene)
        return str(self.calculated_fitness)
    
    def __lt__(self, other):
        return self.calculated_fitness > other.calculated_fitness
        # Low fitness means higher value

    def __random_generate(self):
        self.gene = list(range(1, self.dimension + 1))
        random.seed()
        random.shuffle(self.gene)
        return
    
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
        #print(child)
        if child.dimension != gene_length:
            print("Crossover: The length of the child is different")
            exit(-1)
        return child

    def calc_fitness(self):
        coordinates = self.config["coordinates"]
        sum = 0
        a = None
        b = None
        for i in self.gene:
            if a is None:
                a = coordinates[i]
                continue
            b = coordinates[i]
            sum = sum + a.distance(b, self.config['Distance_table'])
            a = b
        #print(sum)
        return sum

class Pool():
    def __init__(self, config):
        self.config = config
        self.population = config["population"]
        self.dimension = config["dimension"]
        self.pool = list()
        self.create_random_pool()

    def create_random_pool(self):
        for i in range(self.population):
            new_gene = Gene(dimension=self.dimension, config=self.config)
            self.pool.append(new_gene)

    def replace(self, num):
        for i in range(num):
            self.pool.pop(-1)
        return

    def fill_random_gene(self, num):
        for i in range(num):
            new_gene = Gene(dimension=self.dimension, config=self.config)
            #print(new_gene)
            self.pool.append(new_gene)
        return

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
        while(True):
            print("Generation: ", generation)

            # Sort the genes of the pool high to low fitness
            self.pool.sort(reverse=True)
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
    build_distance_table(config)
    result = ga(coordinates, config)

def build_distance_table(config):
    dimension = config["dimension"]
    table = [[None for col in range(dimension)] for row in range(dimension)]
    config['Distance_table'] = table
    return

def ga(coordinates, config):
    config["coordinates"] = coordinates
    new_pool = Pool(config)
    new_pool.generation()

#################################################################
def read_file(config):
    coordinates = dict()
    with open(config["file_name"], "r") as f:
        name = f.readline()         # NAME
        comment = f.readline()      # COMMENT
        type_ = f.readline()        # TYPE
        dimension = f.readline()    # DIMENSION
        edge_weight_type = f.readline() # EDGE_WEIGHT_TYPE
        _ = f.readline()            # NODE_COORD_SECTION

        dimension = dimension.split()
        dimension = dimension[2]
        config["dimension"] = int(dimension)
        for i in range(int(dimension)):
            string = f.readline()
            string = string.split()
            num = int(string[0])
            x = float(string[1])
            y = float(string[2])
            #d_print(str(num) + ", " + str(x) + ", " + str(y))
            coordinates[num] = Coordinate(num, x, y)
        if("EOF\n" != f.readline()): # EOF
            raise Exception
    return coordinates

def d_print(string):
    if config["debug"]:
        print(str(string))

def define_argparse():
    parser = argparse.ArgumentParser(description="Travelling Salesman Problem solver")
    parser.add_argument("Target_file", type=str, help="Target file name")
    parser.add_argument('-p', type=int, help="Population", dest="population", default=50)
    parser.add_argument('-f', type=int, help="Fitness evaluations", dest="fitness_evaluations", default=10)
    return parser

def parse_options(config, parser):
    args = parser.parse_args()
    config["file_name"] = args.Target_file
    config["population"] = args.population
    config["fitness_evalutions"] = args.fitness_evaluations
    return config

if __name__ == "__main__":
    parser = define_argparse()
    config = parse_options(config, parser)
    tsp(config)