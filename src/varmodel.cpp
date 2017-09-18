#include "varmodel.hpp"
#include "PopulationManager.hpp"
#include "HostManager.hpp"
#include "InfectionManager.hpp"
#include "ImmunityManager.hpp"
#include "random.hpp"
#include "parameters.hpp"
#include "EventQueue.hpp"

namespace varmodel {

#pragma mark \
*** Simulation state ***

rng_t * rng;

double current_time;
double next_checkpoint_time;

PopulationManager * population_manager;
HostManager * host_manager;
InfectionManager * infection_manager;
ImmunityManager * immunity_manager;

#pragma mark \
*** Event queues for different types of event ***

bool compare_biting_time(Population * p1, Population * p2) {
    if(p1->next_biting_time == p2->next_biting_time) {
        return p1->id < p2->id;
    }
    return p1->next_biting_time < p2->next_biting_time;
}

typedef EventQueue<Population, compare_biting_time> BitingQueue; 
BitingQueue * biting_queue;

bool compare_immigration_time(Population * p1, Population * p2) {
    if(p1->next_immigration_time == p2->next_immigration_time) {
        return p1->id < p2->id;
    }
    return p1->next_immigration_time < p2->next_immigration_time;
}

typedef EventQueue<Population, compare_immigration_time> ImmigrationQueue;
ImmigrationQueue * immigration_queue;

bool compare_death_time(Host * h1, Host * h2)
{
    if(h1->death_time == h2->death_time) {
        return h1->id < h2->id;
    }
    return h1->death_time < h2->death_time;
}

typedef EventQueue<Host, compare_death_time> DeathQueue;
DeathQueue * death_queue;

#pragma mark \
*** Helper function declarations ***

void do_next_event(); 
void do_checkpoint_event();
void do_biting_event();
void do_immigration_event();
void do_death_event();

void initialize_event_queues();
void initialize_biting_queue();
void initialize_immigration_queue();
void initialize_death_queue();

void initialize_populations();
void initialize_population(uint64_t index);
void initialize_population_events(Population * pop);
void initialize_population_hosts(Population * pop);

double draw_biting_time(Population * pop);
double draw_immigration_time(Population * pop);

double draw_exponential(double lambda);

#pragma mark \
*** Main simulation loop ***

enum class EventType {
    NONE,
    CHECKPOINT,
    BITING,
    IMMIGRATION,
    DEATH
};

void run() {
//    Population * pop = population_manager->object_for_id(0);
//    printf("Population %llu\n", pop->id);
//    for(Host * host : pop->hosts.as_vector()) {
//        printf("Host %llu\n", host->id);
//    }
//    
//    while(current_time < T_END) { 
//        do_next_event();
//    }
}

void do_next_event() {
    double next_event_time = std::numeric_limits<double>::infinity();
    EventType next_event_type = EventType::NONE;
    
    // Find the event type with the lowest-timed event (simple linear search)
    
    if(next_checkpoint_time < next_event_time) {
        next_event_time = next_checkpoint_time;
        next_event_type = EventType::CHECKPOINT;
    }
    
    if(biting_queue->size() > 0 && biting_queue->head()->next_biting_time < next_event_time) {
        next_event_time = biting_queue->head()->next_biting_time;
        next_event_type = EventType::BITING; 
    }
    
    if(immigration_queue->size() > 0 && immigration_queue->head()->next_immigration_time < next_event_time) {
        next_event_time = biting_queue->head()->next_immigration_time;
        next_event_type = EventType::IMMIGRATION;
    }
    
    if(death_queue->size() > 0 && death_queue->head()->death_time < next_event_time) {
        next_event_time = death_queue->head()->death_time;
        next_event_type = EventType::DEATH;
    }
    
    switch(next_event_type) {
        case EventType::NONE:
            break;
        case EventType::CHECKPOINT:
            do_checkpoint_event();
            break;
        case EventType::BITING:
            do_biting_event();
            break;
        case EventType::IMMIGRATION:
            do_immigration_event();
            break;
        case EventType::DEATH:
            do_death_event();
            break;
    }
}

void do_checkpoint_event() {
    printf("do_checkpoint_event()\n");
    save_checkpoint();
    next_checkpoint_time += CHECKPOINT_SAVE_PERIOD;
}

void do_biting_event() {
    printf("do_biting_event()\n");
    assert(biting_queue->size() > 0);
    Population * pop = biting_queue->head();
    pop->next_biting_time += draw_biting_time(pop);
}

void do_immigration_event() {
    printf("do_immigration_event()\n");
}

void do_death_event() {
    printf("do_death_event()\n");
}

#pragma mark \
*** Initialization function implementations ***

void initialize() {
    current_time = 0.0;
    next_checkpoint_time = CHECKPOINT_SAVE_PERIOD;
    
    // Don't use the "new" keyword anywhere except here to create object managers.
    // These object managers handle allocating/freeing memory for all other objects.
    rng = new rng_t(RANDOM_SEED);
    
    population_manager = new PopulationManager();
    host_manager = new HostManager();
    infection_manager = new InfectionManager();
    immunity_manager = new ImmunityManager();
    
    initialize_populations();
    initialize_event_queues();
}

void initialize_populations() {
    for(uint64_t i = 0; i < N_POPULATIONS; i++) {
        initialize_population(i);
    }
}

void initialize_event_queues() {
    initialize_biting_queue();
    initialize_immigration_queue();
    initialize_death_queue();
}

void initialize_biting_queue() {
    biting_queue = new BitingQueue();
}

void initialize_immigration_queue() {
    immigration_queue = new ImmigrationQueue();
}

void initialize_death_queue() {
    death_queue = new DeathQueue();
}

void initialize_population(uint64_t index) {
    Population * pop = population_manager->create();
    pop->index = index;
    pop->transmission_count = 0;
    
    initialize_population_hosts(pop);
}

void initialize_population_hosts(Population * pop) {
    for(uint64_t i = 0; i < N_HOSTS[pop->index]; i++) {
        pop->hosts.add(host_manager->create());
    } 
}


#pragma mark \
*** Checkpoint function implementations ***

void save_checkpoint() {
    sqlite3 * db;
    int result;
    
    result = sqlite3_open(CHECKPOINT_SAVE_FILENAME.c_str(), &db);
    assert(result == SQLITE_OK);
    
    result = sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, NULL);
    assert(result == SQLITE_OK);
    
    population_manager->save_to_checkpoint(db);
    host_manager->save_to_checkpoint(db);
    infection_manager->save_to_checkpoint(db);
    immunity_manager->save_to_checkpoint(db);
    
    result = sqlite3_exec(db, "END TRANSACTION", NULL, NULL, NULL);
    assert(result == SQLITE_OK);
    
    result = sqlite3_close(db);
    assert(result == SQLITE_OK);
}

void load_checkpoint() {
}


#pragma mark \
*** Random draw helper function implementations *** 

double draw_exponential(double lambda) {
    return std::exponential_distribution<>(lambda)(*rng); 
}

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

