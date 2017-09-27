#ifndef Strain_hpp
#define Strain_hpp

#include <stdint.h>
#include <vector>

namespace varmodel {

struct Gene;

struct Strain {
    Strain(uint64_t id) : id(id) { }
    
    uint64_t const id;
    std::vector<Gene *> genes_sorted;
};

} // namespace varmodel

#endif // #ifndef Strain_hpp
