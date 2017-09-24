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

#pragma mark \
*** Simulation state ***

rng_t * rng;

double current_time;
double next_checkpoint_time;

StrainManager * strain_manager;
std::unordered_map<std::vector<Gene *>, Strain *, HashVector<Gene *>> genes_strain_map;

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

bool compare_death_time(Host * h1, Host * h2) {
    if(h1->death_time == h2->death_time) {
        return h1->id < h2->id;
    }
    return h1->death_time < h2->death_time;
}
typedef EventQueue<Host, compare_death_time> DeathQueue;
DeathQueue * death_queue;

bool compare_transition_time(Infection * o1, Infection * o2) {
    if(o1->transition_time == o2->transition_time) {
        return o1->id < o2->id;
    }
    return o1->transition_time < o1->transition_time;
}
typedef EventQueue<Infection, compare_transition_time> TransitionQueue;
TransitionQueue * transition_queue;

bool compare_mutation_time(Infection * o1, Infection * o2) {
    if(o1->mutation_time == o2->mutation_time) {
        return o1->id < o2->id;
    }
    return o1->mutation_time < o1->mutation_time;
}
typedef EventQueue<Infection, compare_mutation_time> MutationQueue;
MutationQueue * mutation_queue;

bool compare_recombination_time(Infection * o1, Infection * o2) {
    if(o1->recombination_time == o2->recombination_time) {
        return o1->id < o2->id;
    }
    return o1->recombination_time < o1->recombination_time;
}
typedef EventQueue<Infection, compare_recombination_time> RecombinationQueue;
RecombinationQueue * recombination_queue;

bool compare_clearance_time(Infection * o1, Infection * o2) {
    if(o1->clearance_time == o2->clearance_time) {
        return o1->id < o2->id;
    }
    return o1->clearance_time < o1->clearance_time;
}
typedef EventQueue<Infection, compare_clearance_time> ClearanceQueue;
ClearanceQueue * clearance_queue;

#pragma mark \
*** Helper function declarations ***

void initialize_gene_pool();

void initialize_populations();
void initialize_population(uint64_t order);
void initialize_population_events(Population * pop);
void initialize_population_hosts(Population * pop);
void initialize_population_infections(Population * pop);

void destroy_host(Host * host);
Host * create_new_host(Population * pop, bool initial);
Gene * get_gene_with_alleles(std::vector<uint64_t> const & alleles);
Strain * generate_random_strain();
Strain * create_strain(std::vector<Gene *> const & genes);
Gene * draw_random_gene();
Strain * get_strain_with_genes(std::vector<Gene *> genes);

void infect_host(Host * host, Strain * strain);

double draw_transition_time(Infection * infection);
double draw_mutation_time(Infection * infection);
double draw_recombination_time(Infection * infection);

void do_next_event(); 

void do_checkpoint_event();

void do_biting_event();
void do_immigration_event();

void do_death_event();

void do_transition_event();
void do_mutation_event();
void do_recombination_event();
void do_clearance_event(); 

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
    
    biting_queue = new BitingQueue();
    immigration_queue = new ImmigrationQueue();
    death_queue = new DeathQueue();
    transition_queue = new TransitionQueue();
    mutation_queue = new MutationQueue();
    recombination_queue = new RecombinationQueue();
    clearance_queue = new ClearanceQueue();
    
    initialize_gene_pool();
    initialize_populations();
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

void initialize_population(uint64_t order) {
    Population * pop = population_manager->create();
    pop->order = order;
    pop->transmission_count = 0;
    pop->next_biting_time = draw_biting_time(pop);
    pop->next_immigration_time = draw_immigration_time(pop);
    
    biting_queue->add(pop);
    immigration_queue->add(pop);
    
    initialize_population_hosts(pop);
    initialize_population_infections(pop);
}

void initialize_population_hosts(Population * pop) {
    for(uint64_t i = 0; i < N_HOSTS[pop->order]; i++) {
        create_new_host(pop, true);
    }
}

void initialize_population_infections(Population * pop) {
    for(uint64_t i = 0; i < N_INITIAL_INFECTIONS[pop->order]; i++) {
        Host * host = pop->hosts.object_at_index(draw_uniform_index(pop->hosts.size())); 
        Strain * strain = generate_random_strain();
        infect_host(host, strain);
    }
}


#pragma mark \
*** Object management functions ***

Host * create_new_host(Population * pop, bool initial) {
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
    if(initial) {
        host->birth_time = -draw_uniform_real(0, lifetime);
    }
    else {
        host->birth_time = current_time;
    }
    host->death_time = host->birth_time + lifetime;
    death_queue->add(host);
    
    pop->hosts.add(host);
    
    return host;
}

void destroy_host(Host * host) {
    Population * pop = host->population;
    printf("Removing host id %llu from population %llu\n", host->id, pop->id);
    
    for(Infection * infection : host->infections) {
        infection_manager->destroy(infection);
    }
    
    for(Immunity * immunity : host->immunities) {
        immunity_manager->destroy(immunity);
    }
    
    pop->hosts.remove(host);
    death_queue->remove(host);
    host_manager->destroy(host);
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

Strain * generate_random_strain() {
    std::vector<Gene *> genes(N_GENES_PER_STRAIN);
    for(uint64_t i = 0; i < N_GENES_PER_STRAIN; i++) {
        genes[i] = draw_random_gene();
    }
    return get_strain_with_genes(genes);
}

Strain * get_strain_with_genes(std::vector<Gene *> genes) {
    std::sort(genes.begin(), genes.end());
    auto itr = genes_strain_map.find(genes);
    if(itr == genes_strain_map.end()) {
        Strain * strain = strain_manager->create();
        for(Gene * gene : genes) {
            strain->genes.add(gene);
        }
        genes_strain_map[genes] = strain;
        return strain;
    }
    return itr->second;
}

Gene * draw_random_gene() {
    return gene_manager->objects()[draw_uniform_index(gene_manager->size())];
}

void infect_host(Host * host, Strain * strain) {
    Infection * infection = infection_manager->create();
    infection->host = host;
    infection->strain = strain;
    
    for(Gene * gene : strain->genes.as_vector()) {
        infection->expression_order.push_back(gene);
    }
    std::shuffle(
        infection->expression_order.begin(), infection->expression_order.end(), *rng
    );
    
    infection->active = false;
    if(T_LIVER_STAGE == 0.0) {
        infection->expression_index = 0;
        infection->transition_time = draw_transition_time(infection);
    }
    else {
        infection->expression_index = -1;
        infection->transition_time = current_time + T_LIVER_STAGE;
    }
    
    infection->mutation_time = draw_mutation_time(infection);
    infection->recombination_time = draw_recombination_time(infection);
    
    host->infections.insert(infection);
}

double draw_transition_time(Infection * infection) {
    if(infection->active) {
        return current_time + draw_exponential(DEACTIVATION_RATE);
    }
    else {
        return current_time + draw_exponential(ACTIVATION_RATE);
    }
}

double draw_mutation_time(Infection * infection) {
    
}

double draw_recombination_time(Infection * infection) {
}

double draw_clearance_time(Infection * infection) {
    return current_time + draw_exponential(CLEARANCE_RATE);
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
    
    printf("biting event pop: %llu\n", pop->id);
    
    // TODO: biting activity
    
    // Update biting event time
    pop->next_biting_time += draw_biting_time(pop);
    biting_queue->update(pop);
}

void do_immigration_event() {
    printf("do_immigration_event()\n");
    assert(immigration_queue->size() > 0);
    Population * pop = immigration_queue->head();
    
    printf("immigration event pop: %llu\n", pop->id);
    
    // TODO: immigration activity
    
    // Update immigration event time
    pop->next_immigration_time += draw_immigration_time(pop);
    immigration_queue->update(pop);
}

void do_death_event() {
    printf("do_death_event()\n");
    assert(death_queue->size() > 0);
    
    Host * host = death_queue->head();
    Population * pop = host->population;
    destroy_host(host);
    Host * new_host = create_new_host(pop, false);
    printf("Created new host id %llu in population %llu\n", new_host->id, pop->id);
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
    double biting_rate = BITING_RATE_MEAN[pop->order] * (
        1.0 + BITING_RATE_RELATIVE_AMPLITUDE[pop->order] * cos(
            2 * M_PI * ((current_time / T_YEAR) - BITING_RATE_PEAK_PHASE[pop->order])
        )
    );
    
    return current_time + draw_exponential(biting_rate);
}

double draw_immigration_time(Population * pop) {
    return current_time + draw_exponential(IMMIGRATION_RATE[pop->order]);
}

} // namespace varmodel

