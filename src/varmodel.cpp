#include "varmodel.hpp"
#include "state.hpp"

namespace varmodel {

void initialize() {
    host_manager = new HostManager();
    population_manager = new PopulationManager();
    event_manager = new EventManager();
    
    checkpoint_event = event_manager->create_checkpoint_event();
    extern BitingRateUpdateEvent * biting_rate_update_event;
}

void load_checkpoint() {
}

void save_checkpoint() {
}

void run() {
}

void update_biting_rate() {
}

} // namespace varmodel

