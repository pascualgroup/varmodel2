#include "varmodel.hpp"

#include <unistd.h>

#include "StrainManager.hpp"
#include "GeneManager.hpp"
#include "PopulationManager.hpp"
#include "HostManager.hpp"
#include "InfectionManager.hpp"
#include "GeneImmuneHistoryManager.hpp"
#include "AlleleImmuneHistoryManager.hpp"
#include "LocusImmunityManager.hpp"
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
GeneImmuneHistoryManager * gene_immune_history_manager;
AlleleImmuneHistoryManager * allele_immune_history_manager;
LocusImmunityManager * locus_immunity_manager;

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

bool compare_immunity_loss_time(Host * o1, Host * o2) {
    if(o1->immunity_loss_time == o2->immunity_loss_time) {
        return o1->id < o2->id;
    }
    return o1->immunity_loss_time < o1->immunity_loss_time;
}
typedef EventQueue<Host, compare_immunity_loss_time> ImmunityLossQueue;
ImmunityLossQueue * immunity_loss_queue;

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
Host * create_host(Population * pop, bool initial);
Gene * get_gene_with_alleles(std::vector<uint64_t> const & alleles);
Strain * generate_random_strain();
Strain * create_strain(std::vector<Gene *> const & genes);
Gene * draw_random_gene();
Strain * get_strain_with_genes(std::vector<Gene *> genes);

Gene * get_current_gene(Infection * infection);

uint64_t get_immune_allele_count(Host * host);

void gain_immunity(Host * host, Gene * gene, bool is_clearing);
void gain_allele_specific_immunity(Host * host, Gene * gene);
double draw_allele_specific_immunity_loss_time(Host * host);
void gain_gene_specific_immunity(Host * host, Gene * gene);
double draw_gene_specific_immunity_loss_time(Host * host);
void gain_general_immunity(Host * host, Gene * gene);
double draw_general_immunity_loss_time(Host * host);

void lose_random_allele_specific_immunity(Host * host);
void lose_gene_specific_immunity(Host * host);
void lose_general_immunity(Host * host);

void update_host_infection_times(Host * host);
void update_infection_times(Infection * infection);

double get_specific_immunity_level(Host * host, Gene * gene);
uint64_t get_active_infection_count(Host * host);
void infect_host(Host * host, Strain * strain);

void perform_infection_transition(Infection * infection);
void clear_infection(Infection * infection);

double draw_transition_time(Infection * infection);
double draw_activation_time(Infection * infection);
double draw_deactivation_time(Infection * infection);

double draw_mutation_time(Infection * infection);
double draw_recombination_time(Infection * infection);
double draw_clearance_time(Infection * infection);

void do_next_event(); 

void do_checkpoint_event();

void do_biting_event();
void do_immigration_event();

void do_immunity_loss_event();
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
    gene_immune_history_manager = new GeneImmuneHistoryManager();
    allele_immune_history_manager = new AlleleImmuneHistoryManager();
    locus_immunity_manager = new LocusImmunityManager();
    
    biting_queue = new BitingQueue();
    immigration_queue = new ImmigrationQueue();
    
    immunity_loss_queue = new ImmunityLossQueue();
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
        create_host(pop, true);
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

Host * create_host(Population * pop, bool initial) {
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
    
    host->immunity_loss_time = std::numeric_limits<double>::infinity();
    
    if(SELECTION_MODE == ALLELE_SPECIFIC_IMMUNITY) {
        AlleleImmuneHistory * allele_immune_history = allele_immune_history_manager->create();
        for(uint64_t i = 0; i < N_LOCI; i++) {
            allele_immune_history->immunity_by_locus.push_back(locus_immunity_manager->create());
        }
        host->allele_immune_history = allele_immune_history;
    }
    else {
        host->allele_immune_history = NULL;
    }
    
    if(SELECTION_MODE == GENE_SPECIFIC_IMMUNITY) {
        host->gene_immune_history = gene_immune_history_manager->create();
    }
    else {
        host->gene_immune_history = NULL;
    }
    
    pop->hosts.add(host);
    
    return host;
}

void destroy_host(Host * host) {
    Population * pop = host->population;
    printf("Removing host id %llu from population %llu\n", host->id, pop->id);
    
    for(Infection * infection : host->infections) {
        infection_manager->destroy(infection);
    }
    
    if(SELECTION_MODE == ALLELE_SPECIFIC_IMMUNITY) {
        for(LocusImmunity * immunity : host->allele_immune_history->immunity_by_locus) {
            locus_immunity_manager->destroy(immunity);
        }
        allele_immune_history_manager->destroy(host->allele_immune_history);
    }
    else if(SELECTION_MODE == GENE_SPECIFIC_IMMUNITY) {
        gene_immune_history_manager->destroy(host->gene_immune_history);
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
    std::sort(genes.begin(), genes.end(),
        [](Gene * o1, Gene * o2) {
            return o1->id < o2->id;
        }
    );
    auto itr = genes_strain_map.find(genes);
    if(itr == genes_strain_map.end()) {
        Strain * strain = strain_manager->create();
        strain->genes_sorted = genes;
        genes_strain_map[genes] = strain;
        return strain;
    }
    return itr->second;
}

Gene * draw_random_gene() {
    return gene_manager->objects()[draw_uniform_index(gene_manager->size())];
}

Gene * get_current_gene(Infection * infection) {
    return infection->expression_order[infection->expression_index];
}

void gain_immunity(Host * host, Gene * gene, bool is_clearing) {
    switch(SELECTION_MODE) {
        case ALLELE_SPECIFIC_IMMUNITY:
            gain_allele_specific_immunity(host, gene);
            break;
        case GENE_SPECIFIC_IMMUNITY:
            gain_gene_specific_immunity(host, gene);
            break;
        case GENERAL_IMMUNITY:
            if(is_clearing) {
                gain_general_immunity(host, gene);
            }
            break;
        case NEUTRALITY:
            break;
    }
}

void gain_allele_specific_immunity(Host * host, Gene * gene) {
    AlleleImmuneHistory * immune_history = host->allele_immune_history;
    for(uint64_t i = 0; i < N_LOCI; i++) {
        LocusImmunity * immunity = immune_history->immunity_by_locus[i];
        auto itr = immunity->immunity_level_by_allele.find(gene->alleles[i]);
        if(itr == immunity->immunity_level_by_allele.end()) {
            immunity->immunity_level_by_allele[gene->alleles[i]] = 1;
        }
        else {
            immunity->immunity_level_by_allele[gene->alleles[i]]++;
        }
    }
    host->immunity_loss_time = draw_allele_specific_immunity_loss_time(host);
}

double draw_allele_specific_immunity_loss_time(Host * host) {
    uint64_t immune_allele_count = get_immune_allele_count(host);
    if(immune_allele_count == 0) {
        return std::numeric_limits<double>::infinity();
    }
    return current_time + draw_exponential(IMMUNITY_LOSS_RATE * immune_allele_count);
}

void gain_gene_specific_immunity(Host * host, Gene * gene) {
}

void gain_general_immunity(Host * host, Gene * gene) {
}

void lose_random_allele_specific_immunity(Host * host) {
    AlleleImmuneHistory * immune_history = host->allele_immune_history;
    uint64_t cur_index = 0;
    uint64_t index = draw_uniform_index(get_immune_allele_count(host));
    bool found = false;
    for(uint64_t i = 0; i < N_LOCI && !found; i++) {
        LocusImmunity * immunity = immune_history->immunity_by_locus[i];
        for(auto kv : immunity->immunity_level_by_allele) {
            if(index == cur_index) {
                assert(kv.second > 0);
                if(kv.second == 1) {
                    immunity->immunity_level_by_allele.erase(kv.first);
                }
                else {
                    immunity->immunity_level_by_allele[kv.first]--;
                }
                found = true;
            }
            if(found) {
                break;
            }
            cur_index++;
        }
    }
    assert(found);
    
    host->immunity_loss_time = draw_allele_specific_immunity_loss_time(host);
    immunity_loss_queue->update(host);
}

uint64_t get_immune_allele_count(Host * host) {
    uint64_t immune_allele_count = 0;
    AlleleImmuneHistory * immune_history = host->allele_immune_history;
    for(uint64_t i = 0; i < N_LOCI; i++) {
        LocusImmunity * immunity = immune_history->immunity_by_locus[i];
        immune_allele_count += immunity->immunity_level_by_allele.size();
    }
    return immune_allele_count;
}

void lose_gene_specific_immunity(Host * host) {
}

void lose_general_immunity(Host * host) {
}

double get_specific_immunity_level(Host * host, Gene * gene) {
    uint64_t immunity_count = 0;
    for(uint64_t i = 0; i < N_LOCI; i++) {
        auto immunity_level_by_allele = host->allele_immune_history->immunity_by_locus[i]->immunity_level_by_allele; 
        auto itr = immunity_level_by_allele.find(gene->alleles[i]);
        if(itr != immunity_level_by_allele.end()) {
            immunity_count += std::min(itr->second, 1ULL);
        }
    }
    return immunity_count / (double)N_LOCI;
}

uint64_t get_active_infection_count(Host * host) {
    uint64_t count = 0;
    for(Infection * infection : host->infections) {
        if(infection->active) {
            count += 1;
        }
    }
    return count;
}

void infect_host(Host * host, Strain * strain) {
    Infection * infection = infection_manager->create();
    infection->host = host;
    infection->strain = strain;
    
    for(Gene * gene : strain->genes_sorted) {
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
    mutation_queue->add(infection);
    
    infection->recombination_time = draw_recombination_time(infection);
    recombination_queue->add(infection);
    
    infection->clearance_time = draw_clearance_time(infection);
    clearance_queue->add(infection);
    
    transition_queue->add(infection);
    
    
    host->infections.insert(infection);
}

void perform_infection_transition(Infection * infection) {
    Host * host = infection->host;
    if(infection->expression_index == -1) {
        assert(!infection->active);
        infection->expression_index = 0;
    }
    else if(infection->active) {
        assert(
            infection->expression_index >= 0 &&
            infection->expression_index < N_GENES_PER_STRAIN - 1
        );
        Gene * gene = get_current_gene(infection);
        gain_immunity(host, gene, false);
        infection->expression_index++;
        infection->active = false;
    }
    else {
        infection->active = true;
    }
    
    update_host_infection_times(infection->host);
}

void clear_infection(Infection * infection) {
    transition_queue->remove(infection);
    mutation_queue->remove(infection);
    recombination_queue->remove(infection);
    clearance_queue->remove(infection);
    
    Host * host = infection->host;
    host->infections.erase(infection);
    
    infection_manager->destroy(infection);
}

void update_host_infection_times(Host * host) {
    for(Infection * infection : host->infections) {
        update_infection_times(infection);
    }
}

void update_infection_times(Infection * infection) {
    if(infection->expression_index >= 0) {
        infection->transition_time = draw_transition_time(infection);
        transition_queue->update(infection);
    }
}

double draw_transition_time(Infection * infection) {
    if(infection->active) {
        return draw_deactivation_time(infection);
    }
    else {
        return draw_activation_time(infection);
    }
}

double draw_activation_time(Infection * infection) {
    if(get_active_infection_count(infection->host) == 0) {
        return current_time;
    }
    else {
        return current_time + draw_exponential(ACTIVATION_RATE);
    }
}

double draw_deactivation_time(Infection * infection) {
    assert(infection->active);
    Host * host = infection->host;
    double rate;
    if(SELECTION_MODE == ALLELE_SPECIFIC_IMMUNITY) {
        Gene * active_gene = get_current_gene(infection);
        double immunity_level =  get_specific_immunity_level(host, active_gene);
        assert(immunity_level >= 0.0 && immunity_level <= 1.0);
        rate = DEACTIVATION_RATE_IMMUNE * DEACTIVATION_RATE_NOT_IMMUNE / (
            DEACTIVATION_RATE_IMMUNE * (1.0 - immunity_level) +
            DEACTIVATION_RATE_NOT_IMMUNE * immunity_level
        );
    }
    else {
        rate = DEACTIVATION_RATE_NOT_IMMUNE;
    }
    return current_time + draw_exponential(rate);
}

double draw_mutation_time(Infection * infection) {
    return current_time + draw_exponential(1.0);
}

double draw_recombination_time(Infection * infection) {
    return current_time + draw_exponential(1.0);
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
    IMMUNITY_LOSS,
    DEATH,
    TRANSITION,
    MUTATION,
    RECOMBINATION,
    CLEARANCE
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
    // This is the least elegant part of the reimplementation, but it should
    // be worth it for simplicity and reduced memory usage.
    // If this section is a speed bottleneck, re-evaluate the choice.
    // It shouldn't come even close.
    if(next_checkpoint_time <= next_event_time) {
        next_event_time = next_checkpoint_time;
        next_event_type = EventType::CHECKPOINT;
    }
    if(
        biting_queue->size() > 0 &&
        biting_queue->head()->next_biting_time < next_event_time
    ) {
        next_event_time = biting_queue->head()->next_biting_time;
        next_event_type = EventType::BITING; 
    }
    if(
        immigration_queue->size() > 0 &&
        immigration_queue->head()->next_immigration_time < next_event_time
    ) {
        next_event_time = biting_queue->head()->next_immigration_time;
        next_event_type = EventType::IMMIGRATION;
    }
    if(
        immunity_loss_queue->size() > 0 &&
        immunity_loss_queue->head()->immunity_loss_time < next_event_time
    ) {
        next_event_time = immunity_loss_queue->head()->immunity_loss_time;
        next_event_type = EventType::IMMUNITY_LOSS; 
    }
    if(
        death_queue->size() > 0 &&
        death_queue->head()->death_time < next_event_time
    ) {
        next_event_time = death_queue->head()->death_time;
        next_event_type = EventType::DEATH;
    }
    if(
        transition_queue->size() > 0 &&
        transition_queue->head()->transition_time < next_event_time
    ) {
        next_event_time = transition_queue->head()->transition_time;
        next_event_type = EventType::TRANSITION;
    }
    if(
        mutation_queue->size() > 0 &&
        mutation_queue->head()->mutation_time < next_event_time
    ) {
        next_event_time = death_queue->head()->death_time;
        next_event_type = EventType::MUTATION;
    }
    if(
        recombination_queue->size() > 0 &&
        recombination_queue->head()->recombination_time < next_event_time
    ) {
        next_event_time = recombination_queue->head()->recombination_time;
        next_event_type = EventType::RECOMBINATION;
    }
    if(
        clearance_queue->size() > 0 &&
        clearance_queue->head()->clearance_time < next_event_time
    ) {
        next_event_time = clearance_queue->head()->clearance_time;
        next_event_type = EventType::CLEARANCE;
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
        case EventType::IMMUNITY_LOSS:
            do_immunity_loss_event();
            break;
        case EventType::DEATH:
            do_death_event();
            break;
        case EventType::TRANSITION:
            do_transition_event();
            break;
        case EventType::MUTATION:
            do_mutation_event();
            break;
        case EventType::RECOMBINATION:
            do_recombination_event();
            break;
        case EventType::CLEARANCE:
            do_clearance_event();
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

void do_immunity_loss_event() {
    printf("do_immunity_loss_event()\n");
    assert(immunity_loss_queue->size() > 0);
    
    Host * host = immunity_loss_queue->head();
    
    switch(SELECTION_MODE) {
        case ALLELE_SPECIFIC_IMMUNITY:
            lose_random_allele_specific_immunity(host);
            break;
        case GENE_SPECIFIC_IMMUNITY:
            lose_gene_specific_immunity(host);
            break;
        case GENERAL_IMMUNITY:
            lose_general_immunity(host);
            break;
        case NEUTRALITY:
            assert(false);
            break;
    }
}

void do_death_event() {
    printf("do_death_event()\n");
    assert(death_queue->size() > 0);
    
    Host * host = death_queue->head();
    Population * pop = host->population;
    destroy_host(host);
    Host * new_host = create_host(pop, false);
    printf("Created new host id %llu in population %llu\n", new_host->id, pop->id);
}

void do_transition_event() {
    printf("do_transition_event()\n");
    assert(transition_queue->size() > 0);
    
    Infection * infection = transition_queue->head();
    if(
        infection->expression_index == N_GENES_PER_STRAIN - 1 &&
        infection->active
    ) {
        clear_infection(infection);
    }
    else {
        perform_infection_transition(infection);
    }
}

void do_mutation_event() {
}

void do_recombination_event() {
}

void do_clearance_event() {
    Infection * infection = clearance_queue->head();
    clear_infection(infection);
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
    gene_immune_history_manager->save_to_checkpoint(db);
    allele_immune_history_manager->save_to_checkpoint(db);
    locus_immunity_manager->save_to_checkpoint(db);
    
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

