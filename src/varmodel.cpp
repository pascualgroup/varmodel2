#include "varmodel.hpp"

#include <unistd.h>

#include "StrainManager.hpp"
#include "GeneManager.hpp"
#include "PopulationManager.hpp"
#include "HostManager.hpp"
#include "InfectionManager.hpp"
#include "ImmunityManager.hpp"
#include "random.hpp"
#include "util.hpp"
#include "parameters.hpp"
#include "EventQueue.hpp"

namespace varmodel {

#pragma mark \
*** Enum parameters ***

enum class SelectionMode {
    SPECIFIC_IMMUNITY,
    GENERALIZED_IMMUNITY,
    NEUTRALITY
};
SelectionMode SELECTION_MODE_ENUM;

#pragma mark \
*** Simulation state ***

rng_t * rng;

double current_time;
double next_checkpoint_time;

StrainManager * strain_manager;
GeneManager * gene_manager;

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

void initialize_gene_pool();

void initialize_populations();
void initialize_population(uint64_t index);
void initialize_population_events(Population * pop);
void initialize_population_hosts(Population * pop);

void initialize_event_queues();
void initialize_biting_queue();
void initialize_immigration_queue();
void initialize_death_queue();

Host * create_new_host(Population * pop);
Gene * get_gene_with_alleles(std::vector<uint64_t> const & alleles);

void do_next_event(); 
void do_checkpoint_event();
void do_biting_event();
void do_immigration_event();
void do_death_event();

double draw_biting_time(Population * pop);
double draw_immigration_time(Population * pop);

double draw_exponential(double lambda);
uint64_t draw_uniform_index(uint64_t size);
double draw_uniform_real(double min, double max);

#pragma mark \
*** Initialization function implementations ***

void validate_and_load_parameters() {
    assert(SAMPLE_DB_FILENAME == "" || !file_exists(SAMPLE_DB_FILENAME));
    assert(HOST_SAMPLING_PERIOD > 0.0);
    assert(!SAMPLE_TRANSMISSION_EVENTS || TRANSMISSION_EVENT_SAMPLING_SKIP > 0);
    assert(T_END >= 0.0);
    assert(T_BURNIN >= 0.0);
    assert(T_BURNIN < T_END);
    
    assert(!SAVE_TO_CHECKPOINT || !file_exists(CHECKPOINT_SAVE_FILENAME));
    assert(!SAVE_TO_CHECKPOINT || CHECKPOINT_SAVE_PERIOD > 0.0);
    
    assert(!LOAD_FROM_CHECKPOINT || file_exists(CHECKPOINT_LOAD_FILENAME));
    
    if(SELECTION_MODE == "SPECIFIC_IMMUNITY") {
        SELECTION_MODE_ENUM = SelectionMode::SPECIFIC_IMMUNITY;
    }
    else if(SELECTION_MODE == "GENERALIZED_IMMUNITY") {
        SELECTION_MODE_ENUM = SelectionMode::GENERALIZED_IMMUNITY;
    }
    else if(SELECTION_MODE == "NEUTRALITY") {
        SELECTION_MODE_ENUM = SelectionMode::NEUTRALITY;
    }
    else {
        assert(false);
    }
    
    assert(RANDOM_SEED > 0);
    assert(T_YEAR > 0.0);
    assert(T_END > 0.0);
    
    assert(N_GENES_IN_POOL >= 1);
    assert(N_GENES_PER_STRAIN >= 1);
    assert(N_LOCI >= 1);
    assert(N_ALLELES.size() == N_LOCI);
    for(auto value : N_ALLELES) {
        assert(value >= 1);
    }
    
    assert(GENE_TRANSMISSIBILITY >= 0.0 && GENE_TRANSMISSIBILITY <= 1.0);
    assert(IMMUNITY_LOSS_RATE >= 0.0);
    
    assert(P_MUTATION >= 0.0 && P_MUTATION <= 1.0);
    assert(P_STRAIN_RECOMBINATION >= 0.0 && P_STRAIN_RECOMBINATION <= 1.0);
    assert(T_LIVER_STAGE >= 0.0);
    
    if(USE_HOST_LIFETIME_PDF) {
        assert(add(HOST_LIFETIME_PDF) > 0.0);
        for(auto value : HOST_LIFETIME_PDF) {
            assert(value >= 0.0);
        }
    }
    else {
        assert(MEAN_HOST_LIFETIME > 0.0);
        assert(MAX_HOST_LIFETIME > MEAN_HOST_LIFETIME);
    }
    
    assert(N_POPULATIONS >= 1);
    for(auto value : N_HOSTS) {
        assert(value >= 1);
    }
    
    assert(BITING_RATE_MEAN.size() == N_POPULATIONS);
    for(auto value : BITING_RATE_MEAN) {
        assert(value >= 0.0);
    }
    assert(BITING_RATE_RELATIVE_AMPLITUDE.size() == N_POPULATIONS);
    for(auto value : BITING_RATE_RELATIVE_AMPLITUDE) {
        assert(value >= 0.0 && value <= 1.0);
    }
    assert(BITING_RATE_PEAK_PHASE.size() == N_POPULATIONS);
    for(auto value : BITING_RATE_PEAK_PHASE) {
        assert(value >= 0.0 && value <= 1.0);
    }
    
    assert(IMMIGRATION_RATE.size() == N_POPULATIONS);
    for(auto value : IMMIGRATION_RATE) {
        assert(value >= 0.0);
    }
}

void initialize() {
    current_time = 0.0;
    next_checkpoint_time = CHECKPOINT_SAVE_PERIOD;
    
    // Don't use the "new" keyword anywhere except here to create object managers.
    // These object managers handle allocating/freeing memory for all other objects.
    rng = new rng_t(RANDOM_SEED);
    
    strain_manager = new StrainManager();
    gene_manager = new GeneManager();
    population_manager = new PopulationManager();
    host_manager = new HostManager();
    infection_manager = new InfectionManager();
    immunity_manager = new ImmunityManager();
    
    initialize_gene_pool();
    initialize_populations();
    initialize_event_queues();
}

void initialize_gene_pool() {
    // Create gene pool
    for(uint64_t i = 0; i < N_GENES_IN_POOL; i++) {
        // Draw random alleles different from all other existing genes
        std::vector<uint64_t> alleles(N_LOCI);
        do {
            for(uint64_t j = 0; j < N_LOCI; j++) {
                alleles[j] = draw_uniform_index(N_ALLELES[j]);
            }
        } while(get_gene_with_alleles(alleles) != NULL);
        
        Gene * gene = gene_manager->create();
        gene->alleles = alleles;
    }
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
    
    initialize_population_hosts(pop);
}

void initialize_population_hosts(Population * pop) {
    for(uint64_t i = 0; i < N_HOSTS[pop->index]; i++) {
        create_new_host(pop);
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
    
    for(Host * host : host_manager->objects()) {
        death_queue->add(host);
    }
}


#pragma mark \
*** Object management functions ***

Host * create_new_host(Population * pop) {
    Host * host = host_manager->create();
    host->population = pop;
    
    double lifetime;
    if(USE_HOST_LIFETIME_PDF) {
        assert(false); // TODO: implement
    }
    else {
        lifetime = std::min(
            draw_exponential(1.0 / MEAN_HOST_LIFETIME),
            MAX_HOST_LIFETIME
        );
    }
    host->birth_time = current_time;
    host->death_time = host->birth_time + lifetime;
    
    pop->hosts.add(host);
    
    return host;
}

Gene * get_gene_with_alleles(std::vector<uint64_t> const & alleles) {
    for(Gene * gene : gene_manager->objects()) {
        assert(gene->alleles.size() == N_LOCI);
        if(gene->alleles == alleles) {
            return gene;
        }
    }
    return nullptr;
}

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
    while(current_time < T_END) { 
        do_next_event();
    }
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
    
    printf("next_event_time: %f\n", next_event_time);
    printf("next_event_type: %d\n", next_event_type);
    
    // Execute the next event
    assert(current_time <= next_event_time);
    current_time = next_event_time;
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
    assert(death_queue->size() > 0);
    
    Host * host = death_queue->head();
    Population * pop = host->population;
    printf("Removing host id %llu from population %llu\n", host->id, pop->id);
    
    pop->hosts.remove(host);
    death_queue->remove(host);
    host_manager->destroy(host);
    
    host = create_new_host(pop);
    printf("Created new host id %llu in population %llu\n", host->id, pop->id);
    death_queue->add(host);
}

#pragma mark \
*** Checkpoint function implementations ***

void save_checkpoint() {
    std::string old_checkpoint_filename = CHECKPOINT_SAVE_FILENAME + "-old";
    if(file_exists(CHECKPOINT_SAVE_FILENAME)) {
        assert(!rename(CHECKPOINT_SAVE_FILENAME.c_str(), old_checkpoint_filename.c_str()));
    }
    
    sqlite3 * db;
    assert(!file_exists(CHECKPOINT_SAVE_FILENAME));
    sqlite3_open(CHECKPOINT_SAVE_FILENAME.c_str(), &db);
    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, NULL);
    
    strain_manager->save_to_checkpoint(db);
    gene_manager->save_to_checkpoint(db);
    population_manager->save_to_checkpoint(db);
    host_manager->save_to_checkpoint(db);
    infection_manager->save_to_checkpoint(db);
    immunity_manager->save_to_checkpoint(db);
    
    sqlite3_exec(db, "END TRANSACTION", NULL, NULL, NULL);
    int result = sqlite3_close(db);
    assert(result == SQLITE_OK);
    
    if(file_exists(old_checkpoint_filename)) {
        assert(!unlink(old_checkpoint_filename.c_str()));
    }
}

void load_checkpoint() {
}


#pragma mark \
*** Random draw helper function implementations *** 

double draw_exponential(double lambda) {
    return std::exponential_distribution<>(lambda)(*rng); 
}

uint64_t draw_uniform_index(uint64_t size) {
    return std::uniform_int_distribution<uint64_t>(0, size - 1)(*rng);
}

double draw_uniform_real(double min, double max) {
    return std::uniform_real_distribution<double>(min, max)(*rng);
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

