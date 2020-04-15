#include "../include/genetic_algorithm.hpp"

namespace ai {

genetic_algorithm::genetic_algorithm(
    std::size_t population_size,
    double probability_crossover,
    double probability_mutation)
    :
    _population_size{population_size},
    _probability_crossover{probability_crossover},
    _probability_mutation{probability_mutation}
{

}

genetic_algorithm::~genetic_algorithm()
{

}

void genetic_algorithm::set_initial_population(
    std::vector<individual_ptr> population)
{
    _population = population;
}

void genetic_algorithm::train(
    const genetic_algorithm::stop_condition_t& condition)
{
    do {
        _population = roulette_wheel_selection();

        auto new_genaration = crossover_population(_population);

        mutate_population(new_genaration);

        _population.insert(
            _population.end(),
            new_genaration.begin(),
            new_genaration.end());

    } while(!condition.should_stop(get_population_state()));
}

std::vector<genetic_algorithm::individual_ptr> genetic_algorithm::
    roulette_wheel_selection()
{
    auto population_fitness = get_population_state()._population_fitness;

    auto sum = std::accumulate(
        population_fitness.begin(),
        population_fitness.end(), 0);

    auto normalized_population_fitness =
        std::for_each(
            population_fitness.begin(),
            population_fitness.end(),
            [sum](double v){
                return v/sum;
            });

    auto random_selected = get_random_vector((int)_population_size/2);

    return select_individuals(population_fitness, random_selected);
}

void genetic_algorithm::mutate_population(
    const std::vector<individual_ptr>& population)
{
    for (auto& ind: population)
        ind->mutate();
}

std::vector<genetic_algorithm::individual_ptr> genetic_algorithm::
    crossover_population(
        std::vector<individual_ptr> parents)
{
    std::vector<individual_ptr> new_generation;

    for (std::size_t i = 0; i < _population.size(); ++i) {
        auto children = _population[i]->crossover(*_population[i+1]);
        new_generation.insert(
            new_generation.end(),
            children.begin(),
            children.end());
    }

    return new_generation;
}

std::vector<genetic_algorithm::individual_ptr> genetic_algorithm::
    select_individuals(
        std::vector<double> pop_fitness, std::vector<double> probabilities)
{
    std::vector<individual_ptr> selected_individuals{};

    for (auto prob: probabilities) {
        std::size_t idx = 0;
        auto accum = pop_fitness[0];
        while (accum < prob) {
            accum += pop_fitness[idx+1];
            ++idx;
        };
        selected_individuals.push_back(_population[idx]);
    }

    return selected_individuals;
}

genetic_algorithm::fitness_data_t genetic_algorithm::
    get_population_state() const
{
    fitness_data_t result{};

    for (auto ind: _population) {
        result._population_fitness.push_back(ind->fitness());
    }

    auto accum = std::accumulate(
        result._population_fitness.begin(),
        result._population_fitness.end(), 0);

    result._avg_fitness = accum / _population_size;

    return result;
}


}