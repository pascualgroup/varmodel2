#ifndef Infection_hpp
#define Infection_hpp

#include <stdint.h>
#include <sqlite3.h>

#include <vector>

namespace varmodel {

struct Strain;
struct Host;
struct Gene;

struct Infection { 
    Infection(uint64_t id) : id(id) { }
    
    uint64_t const id;
    
    Strain * strain;
    Host * host;
    
    int64_t expression_index;
    bool active;
    
    double transition_time;
    double mutation_time;
    double recombination_time;
    double clearance_time;
    
    std::vector<Gene *> expression_order;
};

} // namespace varmodel

#endif // #ifndef Infection_hpp
