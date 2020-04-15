#pragma once

#include <vector>
#include <random>

namespace ai {
    static std::random_device rd;
    static std::mt19937 gen(rd());

    double get_random_number(double a = 0.0, double b = 1.0);

    std::vector<double> get_random_vector(
        std::size_t n,
        double a = 0.0,
        double b = 1.0);

    std::size_t get_random_int(std::size_t a, std::size_t b);


}