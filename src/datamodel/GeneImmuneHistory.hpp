#ifndef GeneImmuneHistory_hpp
#define GeneImmuneHistory_hpp

#include <stdint.h>
#include <sqlite3.h>
#include "IndexedSet.hpp"

namespace varmodel {

struct Gene;

struct GeneImmuneHistory {
    GeneImmuneHistory(uint64_t id) : id(id) { }
     
    uint64_t const id;
    IndexedSet<Gene> immune_genes;
};

} // namespace varmodel }

#endif // #define GeneImmuneHistory_hpp
