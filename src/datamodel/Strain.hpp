#ifndef Strain_hpp
#define Strain_hpp

#include <stdint.h>
#include <unordered_set>

namespace varmodel {

struct Gene;

struct Strain {
    Strain(uint64_t id) : id(id) { }
    
    uint64_t const id;
    std::unordered_set<Gene *> genes;
};

} // namespace varmodel

#endif // #ifndef Strain_hpp
