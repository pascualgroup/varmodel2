#ifndef Strain_hpp
#define Strain_hpp

#include <stdint.h>
#include <array>
#include "parameters.hpp"

namespace varmodel {

struct Gene;

struct Strain {
    Strain(uint64_t id) : id(id) { }
    
    uint64_t const id;
    std::array<Gene *, N_GENES_PER_STRAIN> genes_sorted;
};

} // namespace varmodel

#endif // #ifndef Strain_hpp
