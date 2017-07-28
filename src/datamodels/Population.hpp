#ifndef __Population_hpp__
#define __Population_hpp__

#include <sqlite3.h>
#include <vector>
#include "Host.hpp"
#include "IndexedMap.hpp"

namespace varmodel {

struct Population {
    const int64_t id;
    
    // Scalar fields 
    int64_t transmission_count;
    
    // One-to-many relationships
    IndexedMap<Host> hosts;
    
    // Constructor
    Population(int64_t id, int64_t transmission_count);
};

struct PopulationManager {
    int64_t next_id;
    IndexedMap<Population> populations;
    
    // Constructors
    PopulationManager();
    PopulationManager(int64_t next_id);
    
    // Object management
    Population * create();
    Population * create(int64_t id, int64_t transmission_count);
    Population * population_for_id(int64_t id);
    
    // Database management
    void load_from_db(sqlite3 * db);
    void resolve_relationships_from_db(sqlite3 * db);
    void resolve_population_hosts_table(sqlite3 * db);
    void write_to_db(sqlite3 * db);
    void write_population_table(sqlite3 * db);
    void write_population_hosts_table(sqlite3 * db);
};

} // namespace varmodel

#endif // #ifndef __Population_hpp__
