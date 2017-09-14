#ifndef Immunity_hpp
#define Immunity_hpp

#include <sqlite3.h>
#include "IndexedMap.hpp"

struct Host;

namespace varmodel {

struct Immunity { 
    uint64_t const id;
    Host * host;
    
    Immunity(uint64_t id);
};

struct ImmunityManager {
    uint64_t next_id;
    IndexedMap<Immunity> immunities;
    
    ImmunityManager();
    ImmunityManager(uint64_t next_id);
    
    ImmunityManager * create();
    ImmunityManager * create(uint64_t id);
    ImmunityManager * immunity_for_id(uint64_t id);
    
    void load_from_db(sqlite3 * db);
    void write_to_db(sqlite3 * db);
};

} // namespace varmodel }

#endif // #define Immunity_hpp
