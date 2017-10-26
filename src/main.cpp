#include "varmodel.hpp"
#include <sqlite3.h>
#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include "parameters.hpp"
#include "error_handling.hpp"

int main(int argc, char const * argv[]) {
    program_name = argv[0];
    
    // Load seed from parameter if provided
    // If not provided, parameter will be used
    bool override_seed;
    uint64_t random_seed;
    if(argc > 1) {
        override_seed = true;
        
        char const * seed_str = argv[1];
        errno = 0;
        long seed_long = strtol(seed_str, NULL, 0);
        assert(!errno);
        // Make sure it's unsigned and signed 64-bit safe for the hell of it
        assert(seed_long > 0);
        assert(seed_long < std::numeric_limits<int64_t>::max());
        random_seed = seed_long; 
    }
    else {
        override_seed = false;
        random_seed = 0;
    }
    
    sqlite3_config(SQLITE_CONFIG_LOG, handle_sqlite_error, NULL);
    register_signal_handler();
    
    varmodel::run(override_seed, random_seed);
}
