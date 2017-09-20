#ifndef Host_hpp
#define Host_hpp

#include <stdint.h>
#include <sqlite3.h>
#include <unordered_set>

namespace varmodel {

struct Population;
struct Infection;
struct Immunity;

struct Host {
    Host(uint64_t id) : id(id) { }
    
    uint64_t const id;
    
    Population * population;
    
    double birth_time;
    double death_time;
    
    // One-to-many relationships
    std::unordered_set<Infection *> infections;
    std::unordered_set<Immunity *> immunities;
};

} // namespace varmodel

#endif // #ifndef Host_hpp
