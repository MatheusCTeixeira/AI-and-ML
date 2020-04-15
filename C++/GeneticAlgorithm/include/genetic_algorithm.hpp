#pragma once

#include <algorithm>
#include <vector>
#include <memory>

#include "individual.hpp"
#include "random.hpp"

namespace ai {

class genetic_algorithm {
public:
    struct fitness_data_t;
    struct stop_condition_t;

    using individual_ptr = std::shared_ptr<individual>;

    genetic_algorithm(
        std::size_t population_size,
        double probability_crossover,
        double probability_mutation);

    ~genetic_algorithm();

    void set_initial_population(std::vector<individual_ptr> population);
    void train(const stop_condition_t& condition);

protected:
    std::vector<individual_ptr> roulette_wheel_selection();

    void mutate_population(
        const std::vector<individual_ptr>& population);

    std::vector<individual_ptr> crossover_population(
        std::vector<individual_ptr> parents);

    std::vector<individual_ptr> select_individuals(
        std::vector<double> pop_fitness, std::vector<double> prop);

private:
    fitness_data_t get_population_state() const;

    std::vector<individual_ptr> _population;

    std::size_t _population_size;
    std::size_t _generation;

    double _probability_crossover;
    double _probability_mutation;

    std::size_t _max_generation;
};

struct genetic_algorithm::fitness_data_t {
    std::vector<double> _population_fitness;
    double _avg_fitness;
};

struct genetic_algorithm::stop_condition_t {

    virtual bool should_stop(const fitness_data_t&) const = 0;
};


}