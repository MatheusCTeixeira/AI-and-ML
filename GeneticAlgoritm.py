# ---------------------------------------------------------------------------- #
#                 Genetic Algorithm Implementation From Scratch                #
#                                                                              #
#                                                                              #
#  Autor: Matheus Cândido Teixeira                                             #
#  Data: 16/12/2019                                                            #
# ---------------------------------------------------------------------------- #

import random
import math

# Lista de cidades.
cities = {
    "A": (0, 0),
    "B": (0, 100),
    "C": (0, 200),
    "D": (0, 300),
    "E": (0, 400),
    "F": (0, 500),
    "G": (0, 600),
    "H": (0, 700),
    "I": (0, 800),
    "J": (0, 900),
    "K": (0, 1000),
    "L": (0, 1100),
    "M": (0, 1200),
    "N": (0, 1300),
    "O": (0, 1400),
}

# Calcula a distância entre duas cidades.
def distance_between(citySource, cityDestiny):
    xA, yA = cities[citySource]
    xB, yB = cities[cityDestiny]

    # Distância euclidiana.
    return math.sqrt((xB-xA)**2 + (yB-yA)**2)


# Gera o genes. Cada cromossomo é o ID de uma cidade.
def generate_initial_genes(n_cromossomes):
    genes = []

    for i in range(n_cromossomes):
        random_number = random.randint(0, len(cities)-1)
        genes.append(list(cities.keys())[random_number])

    return genes

# Calcula a distância total do conjunto de cidades.
def calc_gene_distance(gene):
    distance = 0

    for i in range(0, len(cities) - 1):
        city_origin = gene[i]
        city_source = gene[i+1]
        distance += distance_between(city_origin, city_source)

    return distance

# Cálcula a aptidão do genes.
def fitness(genes):
    value = 1400                                # Valor máximo.
    duplicated_city_penality = 10000            # Penalidade por duplicação.

    # O número de cidades duplicadas.
    num_duplications = len(cities) - len(set(genes))

    # Aplica penalidade por cidades duplicadas.
    value -= duplicated_city_penality * num_duplications

    # Aplica penalidade por distância.
    value -= calc_gene_distance(genes)

    return value


# Heurística 1
# Buscar cidade que minimiza a distância percorrida.
def heuristic_1(gene, selected_cromossome):
    cromossome = gene[selected_cromossome]
    cities_names = filter(lambda city: city != cromossome, list(cities.keys()))

    fits = []
    for city in cities_names:
        gene[selected_cromossome] = city
        fits.append((city, fitness(gene)))

    fits.sort(key= lambda city: city[1], reverse=True)
    # print(fits)
    gene[selected_cromossome] = cromossome

    return fits[0]

class Individual:
    fitness_func = None # Método estático.

    def __init__(self, gene):
        self.initialize(gene)

    def initialize(self, genes):
        # Inicializa o gene do indivíduo com valores aleatórios.
        self.genes = genes
        self.n_genes = len(genes)

    def mutation(self, prob):
        # Determina se o indivíduo vai sofrer mutação ou não. Se o valor aleató-
        # rio de probabilidade for maior que o específicado por parâmetro, então
        # ele deve sofrer mutação.
        n_crom = random.randint(0, int(prob * (len(self.genes)-1)))

        # if not will_mutate:
            # return
        cities_names = list(cities.keys())
        sz = len(cities_names)

        # Faz o swap entre cromossomos.
        for i in range(n_crom):
            idx = random.randint(0, len(self.genes) - 1)
            self.genes[idx] = heuristic_1(self.genes, idx)[0]

    # Em vez de remover apenas a parte inicial do gene, esse algoritmo remove
    # subtitui apenas uma porção dele.
    def replace_genes_portion(self, _from, _to, _by):
        new_gen = self.genes[:_from]            # Seleciona a primeira parte.
        new_gen.extend(_by)                     # A parte a ser subtituida.
        new_gen.extend(self.genes[_to:])        # Adiciona a parte final.

        # Após esse processo, ambas as partes devem ter o mesmo tamanho.
        assert(len(self.genes) == len(new_gen))

        self.genes = new_gen

    # Retorna uma porção do genes.
    def get_genes_portion(self, _from, _to):
        return self.genes[_from: _to]

    # Realiza a duplicação do genes.
    def mitosis(self):
        genes = self.genes.copy()
        other = Individual(genes)

        return other

    # Calcula a aptidão (fitness) do individuo.
    def get_fitness(self):
        return Individual.fitness_func(self.genes)

class Population:
    def __init__(self,                            \
            n_individuals,                        \
            n_cromossomos,                        \
            max_generations,                      \
            mutation_prob,                        \
            generator_func=generate_initial_genes,\
            fitness_func=fitness):

        self.population = \
            [Individual(generator_func(n_cromossomos)) \
                for i in range(n_cromossomos)]


        self.n_individuals  = n_individuals    # O número de indivíduos.
        self.n_cromossomos  = n_cromossomos    # O número de cromossomos.
        self.generation     = 0                # Geração atual.
        self.max_generation = max_generations  # Quantidade máxima de gerações.
        self.mutation_prob  = mutation_prob    # Probabilidade de mutação.
        self.generator_func = generator_func   # Função para gerar genes.
        self.fitness_func   = fitness_func     # Função para calcular aptidão.

    # Faz a seleção natural, isto é, após o crossover deve haver mais indivíduos
    # que a população suporta. Portanto, o excedente é eliminado, ou seja, os
    # menos aptos.
    def do_selection(self):
        self.population.sort(
            key=lambda individual: individual.get_fitness(),
            reverse=True)

        self.population = self.population[:self.n_individuals]

    def _random_idxs(self, n, min_value, max_value):
        result = [random.randint(min_value, max_value) for x in range(n)]
        result.sort()

        return result

    def do_crossover(self):
        # Seleciona o indivíduos para o crossover.
        self.population.sort(
            key=lambda individual: individual.get_fitness(),
            reverse=True)

        self.descendents = [] # A nova geração de indivíduos.

        # Os 10 melhores reproduzem.
        for i in range(0, min(10, len(self.population) - 1)):
            father = self.population[i]
            mother = self.population[i+1]

            # Escolhe o índice de começo e fim para realizar o crossover.
            low_crom_idx, high_crom_idx =\
                self._random_idxs(2, 0, self.n_cromossomos)

            # Extrai a porção do gene a ser trocada.
            father_dna = father.get_genes_portion(low_crom_idx, high_crom_idx)
            mother_dna = father.get_genes_portion(low_crom_idx, high_crom_idx)

            # Gera o primeiro descente.
            descendent_1 = father.mitosis()
            # Substitui pela porção do DNA do outro ancestral.
            descendent_1.\
                replace_genes_portion(low_crom_idx, high_crom_idx, mother_dna)
            self.descendents.append(descendent_1)

            # Mesmo processo para o segundo.
            descendent_2 = mother.mitosis()
            descendent_2.\
                replace_genes_portion(low_crom_idx, high_crom_idx, father_dna)
            self.descendents.append(descendent_2)

        # Adiciona a nova geração à população.
        # self.population.extend(descendents)

        # Note: Neste ponto há len(descende) a mais. Na fase de seleção os mais
        # aptos deve ser selecionados para a próxima geração.

    # Realiza a mutação. Note que a mutação é individual. Em um ciclo, alguns
    # indivíduos podem ou não sofre mutação.
    def do_mutation(self):
        for individual in self.descendents:
            individual.mutation(self.mutation_prob)

        self.population.extend(self.descendents)


    def next_generation(self):
        self.generation += 1 # Incrementa o ciclo/geração.

        self.do_crossover() # Realiza o crossover.
        self.do_selection() # Faz a "seleção natural".
        self.do_mutation()  # Faz a mutação.

    # Retorna o maior valor de aptidão/fitness contido na população.
    def highest_fitness(self):
        fitness_score_list = [ind.get_fitness() for ind in self.population]
        max_fitness = max(fitness_score_list)

        return max_fitness

    # Calcula o valor médio da aptidão/fitness da popução.
    def average_fitness(self):
        fitness_score_list = [ind.get_fitness() for ind in self.population]
        total_fitness = sum(fitness_score_list)
        average = total_fitness / self.n_individuals

        return average

    # Exibe todos os indivíduos da população.
    def show_all_individual(self):
        self.population.sort(key=lambda ind: ind.get_fitness())
        for individual in self.population:
            fitness = individual.get_fitness()
            distance = calc_gene_distance(individual.genes)

            print("%s -> f: %-2f d:%-2f"%
                  (str(individual.genes), fitness, distance))

    # Itera até atingir algum critério de parada.
    def until_perfection(self, highest_fit, avg_fit):
        h_value, hv_tol = highest_fit
        a_value, av_tol = avg_fit

        absolute_error = lambda a, b: abs(a - b)
        last_highest_value = self.highest_fitness()

        while absolute_error(self.highest_fitness(), h_value) > hv_tol and \
              absolute_error(self.average_fitness(), a_value) > av_tol and \
              self.generation < self.max_generation:

            if self.generation % 500 == 0:
                highest_fitness = self.highest_fitness()
                print("[%-5d] hg: %4.2f, av: %4.2f, prob: %4.2f"%\
                    (self.generation, \
                    highest_fitness, \
                    self.average_fitness(),\
                    self.mutation_prob))

            # Mais uma geração.
            self.next_generation()

        if self.generation == self.max_generation:
            print("Any solution was found.")
        elif absolute_error(self.highest_fitness(), h_value) < hv_tol:
            print("Highest value reached.")
        elif absolute_error(self.average_fitness(), a_value) < av_tol:
            print("Average fitness reached.")


def main():
    n_cromossomos = len(cities)
    sz_population = 200
    max_fitness = 1400

    Individual.fitness_func = fitness
    population = Population(sz_population, n_cromossomos, 5000, 0.40)
    population.until_perfection((0, 50.0), (max_fitness, 1E2))
    population.show_all_individual()

    return 0

def test():
    n_cromossomes = len(cities)

    genes = generate_initial_genes(n_cromossomes)

    print("genes: %s"%genes)
    print("distance: %f"%calc_gene_distance(genes))
    print("fitness: %f"%fitness(genes))
    print("duplications: %d"%(len(genes)-len(set(genes))))

    enhence = heuristic_1(genes, 1)
    genes[1] = enhence[0]
    print("fitness: %f"%fitness(genes))

if __name__ == "__main__":
    main()
