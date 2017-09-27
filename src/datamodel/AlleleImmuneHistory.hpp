#ifndef AlleleImmuneHistory_hpp
#define AlleleImmuneHistory_hpp

#include <stdint.h>
#include <sqlite3.h>
#include <vector>

namespace varmodel {

struct LocusImmunity;

struct AlleleImmuneHistory {
    AlleleImmuneHistory(uint64_t id) : id(id) { }
     
    uint64_t const id;
    std::vector<LocusImmunity *> immunity_by_locus;
};

} // namespace varmodel }

#endif // #define AlleleImmuneHistory_hpp
