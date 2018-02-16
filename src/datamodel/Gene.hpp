#ifndef Gene_hpp
#define Gene_hpp

#include <array>
#include "parameters.hpp"

namespace varmodel {

enum GeneSource {
    SOURCE_POOL = 0,
    SOURCE_RECOMBINATION = 1,
    SOURCE_MUTATION = 2,
    SOURCE_IMMIGRATION = 3
};

struct Gene {
    Gene(uint64_t id) : id(id) { }
    
    uint64_t const id;
    
    uint64_t source;
    uint64_t deactivationClass;
    uint64_t actBiasClass;
    bool is_functional;
    std::array<uint64_t, N_LOCI> alleles;
};

} // namespace varmodel

#endif // #ifndef Gene_hpp
