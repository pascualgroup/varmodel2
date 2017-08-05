#include "varmodel.hpp"
#include "state.hpp"

namespace varmodel {

void initialize_populations();
void initialize_population(uint64_t index);
void initialize_hosts(Population * pop, uint64_t n_hosts);

void initialize() {
    host_manager = new HostManager();
    population_manager = new PopulationManager();
    event_manager = new EventManager();
    
    checkpoint_event = event_manager->create_checkpoint_event();
    
    initialize_populations();
}

void initialize_populations() {
    for(uint64_t i = 0; i < N_POPULATIONS; i++) {
        initialize_population(i);
    }
}

void initialize_population(uint64_t index) {
    double biting_rate = BITING_RATE_MEAN;
    double transmission_count = 0;
    
    Population * pop = population_manager->create(biting_rate, transmission_count);
    initialize_hosts(pop, N_HOSTS[index]);
}

void initialize_hosts(Population * pop, uint64_t n_hosts) {
    for(uint64_t i = 0; i < n_hosts; i++) {
        pop->hosts.add(host_manager->create());
    } 
}

void load_checkpoint() {
}

void save_checkpoint() {
}

void run() {
}

void do_update_biting_rate_event() {
}

void do_biting_event(Population * population) {
}

void do_immigration_event(Population * population) {
}

void do_death_event(Host * host) {
}

void do_transition_event(Infection * infection) {
}

void do_clearance_event(Infection * infection) {
}

} // namespace varmodel

