#ifndef Population_hpp
#define Population_hpp

#include <sqlite3.h>
#include <vector>
#include "Host.hpp"
#include "IndexedMap.hpp"
#include "dbtypes.hpp"

namespace varmodel {

/*** Forward type definitions ***/

struct BitingEvent;
struct ImmigrationEvent;

/*** Population type declaration ***/

struct Population {
    uint64_t const id;
    
    // Scalar fields
    uint64_t index;
    uint64_t transmission_count;
    
    // References to objects
    BitingEvent * biting_event;
    ImmigrationEvent * immigration_event;
    
    // Collections of objects
    IndexedMap<Host> hosts;
};

} // namespace varmodel

#endif // #ifndef Population_hpp
