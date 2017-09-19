#ifndef Strain_hpp
#define Strain_hpp

namespace varmodel {

struct Strain {
    Strain(uint64_t id) : id(id) { }
    
    uint64_t const id;
    
    std::vector<Gene *> genes;
};

} // namespace varmodel

#endif // #ifndef Strain_hpp
