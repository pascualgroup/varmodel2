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
    uint64_t order;
    uint64_t transmission_count;
    double next_biting_time;
    double next_immigration_time;
    
    // Collections of objects
    IndexedSet<Host> hosts;
};

} // namespace varmodel

#endif // #ifndef Population_hpp
