#ifndef ImmuneHistory_hpp
#define ImmuneHistory_hpp

#include <stdint.h>
#include <sqlite3.h>
#include <array>
#include "parameters.hpp"

namespace varmodel {

struct LocusImmunity;

struct ImmuneHistory {
    ImmuneHistory(uint64_t id) : id(id) { }
     
    uint64_t const id;
    std::array<LocusImmunity *, N_LOCI> immunity_by_locus;
};

} // namespace varmodel }

#endif // #define ImmuneHistory_hpp
