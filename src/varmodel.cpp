#include "varmodel.hpp"

#include <iostream>

namespace varmodel {

void run() {
    host_manager = new HostManager();
    population_manager = new PopulationManager();
    
    std::cerr << CHECKPOINT_LOAD_FILENAME << std::endl;
}

} // namespace varmodel

