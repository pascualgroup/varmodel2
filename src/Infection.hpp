#ifndef Infection_hpp
#define Infection_hpp

#include <sqlite3.h>
#include "IndexedMap.hpp"

struct Host;

namespace varmodel {

struct Infection { 
    uint64_t const id;
    Host * host;
    
    Infection(uint64_t id);
};

struct InfectionManager {
    uint64_t next_id;
    IndexedMap<Infection> infections;
    
    InfectionManager();
    InfectionManager(uint64_t next_id);
    
    InfectionManager * create();
    InfectionManager * create(uint64_t id);
    InfectionManager * infection_for_id(uint64_t id);
    
    void load_from_db(sqlite3 * db);
    void write_to_db(sqlite3 * db);
};

} // namespace varmodel

#endif // #ifndef Infection_hpp
