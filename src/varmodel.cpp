#include "varmodel.hpp"
#include "state.hpp"

namespace varmodel {

void initialize() {
    host_manager = new HostManager();
    population_manager = new PopulationManager();
    event_manager = new EventManager();
}

void run() {
}

} // namespace varmodel

