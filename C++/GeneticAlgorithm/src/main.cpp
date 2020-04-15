#include <iostream>
#include <algorithm>
#include "../include/individual.hpp"

using namespace std;


int main()
{

    std::set<char> alleles{{'A', 'B', 'C', 'D', 'E'}};
    auto genes = ai::individual::generate_individual(&alleles, 6, true);
    genes = ai::individual::generate_individual(&alleles, 6, true);
    genes = ai::individual::generate_individual(&alleles, 6, true);

    std::string genes_str{genes.begin(), genes.end()};

    cout << genes_str << endl;
    return 0;
}