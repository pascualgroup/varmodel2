#ifndef PopulationManager_hpp
#define PopulationManager_hpp

#include "Population.hpp"
#include "IndexedMap.hpp"

namespace varmodel {

struct PopulationManager {
    uint64_t next_id;
    IndexedMap<Population> collection;
    
    // Object management
    Population * create();
    Population * create(uint64_t id);
    Population * object_for_id(uint64_t id);
    
    // Checkpointing methods
    void load_from_checkpoint(sqlite3 * db, char const * const table_name);
    void resolve_relationships(sqlite3 * db, BitingEventManager * biting_event_manager, ImmigrationEventManager * immigration_event_manager);
    void save_to_checkpoint(sqlite3 * db, char const * const table_name, char const * const hosts_table_name);
}

} // namespace varmodel

#endif // #define PopulationManager_hpp
