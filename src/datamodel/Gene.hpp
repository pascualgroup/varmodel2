#ifndef Gene_hpp
#define Gene_hpp

#include <vector>

namespace varmodel {

struct Gene {
    Gene(uint64_t id) : id(id) { }
    
    uint64_t const id;
    
    std::vector<uint64_t> alleles;
};

} // namespace varmodel

#endif // #ifndef Gene_hpp
