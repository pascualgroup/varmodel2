#ifndef Infection_hpp
#define Infection_hpp

#include <stdint.h>
#include <sqlite3.h>

namespace varmodel {

struct Strain;
struct Host;

struct Infection { 
    Infection(uint64_t id) : id(id) { }
    
    uint64_t const id;
    
    Strain * strain;
    Host * host;
    
    uint64_t gene_index;
    bool active;
    double transition_time;
};

} // namespace varmodel

#endif // #ifndef Infection_hpp
