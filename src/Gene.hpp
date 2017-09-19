#ifndef Gene_hpp
#define Gene_hpp

namespace varmodel {

struct Gene {
    Strain(uint64_t id) : id(id) { }
    
    uint64_t const id;
    
    double transmissibility;
    double immunity_loss_rate;
};

} // namespace varmodel

#endif // #ifndef Gene_hpp
