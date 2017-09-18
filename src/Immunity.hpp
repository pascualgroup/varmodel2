#ifndef Immunity_hpp
#define Immunity_hpp

#include <stdint.h>
#include <sqlite3.h>

namespace varmodel {

struct Host;

struct Immunity {
    Immunity(uint64_t id) : id(id) { }
     
    uint64_t const id;
    Host * host;
};

} // namespace varmodel }

#endif // #define Immunity_hpp
