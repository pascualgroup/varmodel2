#ifndef Population_hpp
#define Population_hpp

#include <sqlite3.h>
#include "Host.hpp"
#include "IndexedSet.hpp"

namespace varmodel {

/*** Population type declaration ***/

struct Population {
    Population(uint64_t id) : id(id) { }
    
    uint64_t const id;
    
    // Scalar fields
    uint64_t ind;
    uint64_t transmission_count;
    double next_biting_time;
    double next_immigration_time;
    double next_IRS_rate_change_time;
    uint64_t current_IRS_id = 0;
    uint64_t within_IRS_id = 0;
    double IRS_biting_rate = 1;
    double IRS_immigration_rate_factor = 1;
    uint64_t MDA_id = 0;
    double next_MDA_time;
    bool MDA_effective_period = false;
    double MDA_immigration_rate_factor = 1;

    // Collections of objects
    IndexedSet<Host> hosts;
};

} // namespace varmodel

#endif // #ifndef Population_hpp
