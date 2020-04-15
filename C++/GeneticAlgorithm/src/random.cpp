#include "../include/random.hpp"

namespace ai {

double get_random_number(double a, double b)
{
    std::uniform_real_distribution<> dis(a, b);

    return dis(gen);
}

std::vector<double> get_random_vector(
    std::size_t n,
    double a,
    double b)
{
    std::vector<double> output{};

    for (std::size_t i = 0; i < n; ++i)
        output.push_back(get_random_number(a, b));

    return output;
}

std::size_t get_random_int(std::size_t a, std::size_t b)
{
    std::uniform_int_distribution<> dis(a, b);

    return dis(gen);
}

}