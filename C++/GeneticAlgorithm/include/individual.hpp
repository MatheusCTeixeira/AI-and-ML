#pragma once

#include <algorithm>
#include <set>
#include <random>
#include <memory>
#include <vector>

namespace ai {

class individual
{
public:
    using individual_ptr = std::shared_ptr<individual>;

    static std::vector<char> generate_individual(
        std::set<char> *alleles,
        std::size_t locus_number,
        bool allow_gene_duplication = false);

    virtual ~individual();

    virtual std::vector<individual_ptr> crossover(
        const individual& other) const = 0;

    virtual void mutate() = 0;

    virtual double fitness() const = 0;

protected:
    individual(
        std::vector<char> genes,
        double pc,
        double pm);

private:
    std::vector<char> _genes;
    std::set<char> *_alleles;

    double _pc;
    double _pm;
};


}