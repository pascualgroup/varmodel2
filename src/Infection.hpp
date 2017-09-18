#ifndef Infection_hpp
#define Infection_hpp

#include <stdint.h>
#include <sqlite3.h>

namespace varmodel {

struct Host;

struct Infection { 
    Infection(uint64_t id) : id(id) { }
    
    uint64_t const id;
    Host * host;
};

} // namespace varmodel

#endif // #ifndef Infection_hpp
