#ifndef Population_hpp
#define Population_hpp

#include <sqlite3.h>
#include <vector>
#include "Host.hpp"
#include "IndexedMap.hpp"

namespace varmodel {

struct Population {
    uint64_t const id;
    
    // Scalar fields
    double const biting_rate;
    uint64_t transmission_count;
    
    // One-to-many relationships
    IndexedMap<Host> hosts;
    
    // Constructor
    Population(uint64_t id, double biting_rate, uint64_t transmission_count);
};

struct PopulationManager {
    uint64_t next_id;
    IndexedMap<Population> populations;
    
    // Constructors
    PopulationManager();
    PopulationManager(uint64_t next_id);
    
    // Object management
    Population * create(double biting_rate, uint64_t transmission_count);
    Population * create(uint64_t id, double biting_rate, uint64_t transmission_count);
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
