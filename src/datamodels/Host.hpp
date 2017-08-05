#ifndef Host_hpp
#define Host_hpp

#include <sqlite3.h>
#include "Infection.hpp"
#include "IndexedMap.hpp"

namespace varmodel {

struct Host {
    uint64_t const id;
    
    // One-to-many relationships
    IndexedMap<Infection> infections;
    
    Host(uint64_t id);
};

struct HostManager {
    uint64_t next_id;
    IndexedMap<Host> hosts;
    
    HostManager();
    HostManager(uint64_t next_id);
    
    Host * create();
    Host * create(uint64_t id);
    Host * host_for_id(uint64_t id);
    
    void load_from_db(sqlite3 * db);
    void resolve_relationships_from_db(sqlite3 * db);
    
    void write_to_db(sqlite3 * db);
    void write_host_table(sqlite3 * db);
};

} // namespace varmodel

#endif // #ifndef Host_hpp
