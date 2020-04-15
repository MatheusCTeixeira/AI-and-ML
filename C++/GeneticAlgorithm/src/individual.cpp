#include "../include/individual.hpp"

namespace ai {

std::vector<char> individual::generate_individual(
    std::set<char> *alleles,
    std::size_t locus_number,
    bool allow_gene_duplication)
{
    auto n_of_alleles = alleles->size();

    std::vector<char> genes{};

    if (!allow_gene_duplication && alleles->size() < locus_number)
        throw std::invalid_argument{"Few alleles."};

    auto b = alleles->begin();
    while (genes.size() < locus_number) {
        if (b == alleles->end()) b = alleles->begin();
        genes.push_back(*(b++));
    }

    std::random_shuffle(genes.begin(), genes.end());

    return genes;
}


}