#include "varmodel.hpp"
#include "state.hpp"
#include "random.hpp"

#include "Event.hpp"

namespace varmodel {

#pragma mark \
*** Helper function declarations ***

void initialize_populations();
void initialize_population(uint64_t index);
void initialize_population_events(Population * pop);
void initialize_population_hosts(Population * pop);

double draw_biting_time(Population * pop);
double draw_immigration_time(Population * pop);


#pragma mark \
*** Initialization function implementations ***

void initialize() {
    current_time = 0.0;
    
    // Don't use the "new" keyword anywhere except here: these object managers
    // handle allocating/freeing memory for all other objects.
    rng = new rng_t(RANDOM_SEED);
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
    Population * pop = population_manager->create();
    pop->index = index;
    pop->transmission_count = 0;
    
    BitingEvent * biting_event = event_manager->create_biting_event(
        pop, draw_biting_time(pop)
    );
    pop->biting_event = biting_event;
    biting_event->population = pop;
    
    ImmigrationEvent * immigration_event = event_manager->create_immigration_event(
        pop, draw_immigration_time(pop)
    );
    pop->immigration_event = immigration_event;
    initialize_population_hosts(pop);
}

void initialize_population_hosts(Population * pop) {
    for(uint64_t i = 0; i < N_HOSTS[pop->index]; i++) {
        pop->hosts.add(host_manager->create());
    } 
}


#pragma mark \
*** Checkpoint function implementations ***

void load_checkpoint() {
}

void save_checkpoint() {
}

void run() {
}


#pragma mark \
*** Event callback implementations ***

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


#pragma mark \
*** Random draw helper function implementations *** 

double draw_biting_time(Population * pop) {
    double biting_rate = BITING_RATE_MEAN[pop->index] * (
        1.0 + BITING_RATE_RELATIVE_AMPLITUDE[pop->index] * cos(
            2 * M_PI * ((current_time / T_YEAR) - BITING_RATE_PEAK_PHASE[pop->index])
        )
    );
    
    return current_time + draw_exponential(biting_rate);
}

double draw_immigration_time(Population * pop) {
    return current_time + draw_exponential(IMMIGRATION_RATE[pop->index]);
}

} // namespace varmodel

