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

DB_TYPE struct Population {
    uint64_t const id;
    
    // Scalar fields
    DB_FIELD uint64_t index;
    DB_FIELD uint64_t transmission_count;
    
    // References to objects
    DB_REF BitingEvent * biting_event;
    DB_REF ImmigrationEvent * immigration_event;
    
    // Collections of objects
    DB_REFLIST IndexedMap<Host> hosts;
    
    // Constructor
    Population(uint64_t id);
};


/*** PopulationManager type declaration ***/

struct PopulationManager {
    uint64_t next_id;
    IndexedMap<Population> populations;
    
    // Constructors
    PopulationManager();
    PopulationManager(uint64_t next_id);
    
    // Object management
    Population * create();
    Population * create(uint64_t id);
    Population * population_for_id(uint64_t id);
    
    // Database management
    void load_from_db(sqlite3 * db);
    void resolve_relationships_from_db(sqlite3 * db);
    void resolve_population_hosts_table(sqlite3 * db);
    void write_to_db(sqlite3 * db);
    void write_population_table(sqlite3 * db);
    void write_population_hosts_table(sqlite3 * db);
};

} // namespace varmodel

#endif // #ifndef Population_hpp
