#ifndef Gene_hpp
#define Gene_hpp

#include <array>
#include "parameters.hpp"

namespace varmodel {

struct Gene {
    Gene(uint64_t id) : id(id) { }
    
    uint64_t const id;
    
    std::array<uint64_t, N_LOCI> alleles;
};

} // namespace varmodel

#endif // #ifndef Gene_hpp
