#ifndef Infection_hpp
#define Infection_hpp

#include <stdint.h>
#include <sqlite3.h>

#include "parameters.hpp"
#include <array>

namespace varmodel {

struct Strain;
struct Host;
struct Gene;

struct Infection { 
    Infection(uint64_t id) : id(id) { }
    
    uint64_t const id;
    
    Strain * strain;
    Host * host;
    uint64_t hostInfection_id;
    int64_t expression_index;
    
    double transition_time;
    double mutation_time;
    double recombination_time;
    double infected_time;
    double clearance_time;
    
    std::array<Gene *, N_GENES_PER_STRAIN> expression_order;
};

} // namespace varmodel

#endif // #ifndef Infection_hpp
