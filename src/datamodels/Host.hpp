#ifndef __Host_hpp__
#define __Host_hpp__

#include <sqlite3.h>
#include "IndexedMap.hpp"

namespace varmodel {

struct Host { 
    const int64_t id;
    
    Host(int64_t id);
};

struct HostManager {
    int64_t next_id;
    IndexedMap<Host> hosts;
    
    HostManager();
    HostManager(int64_t next_id);
    
    Host * create();
    Host * create(int64_t id);
    Host * host_for_id(int64_t id);
    
    void load_from_db(sqlite3 * db);
    void resolve_relationships_from_db(sqlite3 * db);
    
    void write_to_db(sqlite3 * db);
    void write_host_table(sqlite3 * db);
};

} // namespace varmodel

#endif // #ifndef __Host_hpp__
