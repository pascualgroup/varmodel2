#ifndef Host_hpp
#define Host_hpp

#include <sqlite3.h>
#include "Infection.hpp"
#include "IndexedMap.hpp"

namespace varmodel {

struct Host {
    Host(uint64_t id) : id(id) { }
    
    uint64_t const id;
    
    double birth_time;
    double death_time;
    
    // One-to-many relationships
    IndexedMap<Infection> infections;
};

} // namespace varmodel

#endif // #ifndef Host_hpp
