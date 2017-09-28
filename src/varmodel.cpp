#include "varmodel.hpp"

#include <unistd.h>

#include "StrainManager.hpp"
#include "GeneManager.hpp"
#include "PopulationManager.hpp"
#include "HostManager.hpp"
#include "InfectionManager.hpp"
#include "AlleleImmuneHistoryManager.hpp"
#include "LocusImmunityManager.hpp"
#include "random.hpp"
#include "util.hpp"
#include "parameters.hpp"
#include "EventQueue.hpp"

namespace varmodel {

#pragma mark \
*** SIMULATION STATE ***

double INF = std::numeric_limits<double>::infinity();

rng_t rng(RANDOM_SEED);

double now = 0.0;
double next_verification_time = VERIFICATION_ON ? 0.0 : INF;
double next_checkpoint_time = SAVE_TO_CHECKPOINT ? 0.0 : INF;

StrainManager strain_manager;
std::unordered_map<
    std::array<Gene *, N_GENES_PER_STRAIN>,
    Strain *,
    HashArray<Gene *, N_GENES_PER_STRAIN>
> genes_strain_map;

GeneManager gene_manager;

PopulationManager population_manager;
HostManager host_manager;
InfectionManager infection_manager;
AlleleImmuneHistoryManager immune_history_manager;
LocusImmunityManager locus_immunity_manager;

#pragma mark \
*** EVENT QUEUES ***

double get_biting_time(Population * p) { return p->next_biting_time; }
EventQueue<Population, get_biting_time> biting_queue; 

double get_immigration_time(Population * p) { return p->next_immigration_time; }
EventQueue<Population, get_immigration_time> immigration_queue;

double get_immunity_loss_time(Host * h) { return h->immunity_loss_time; }
EventQueue<Host, get_immunity_loss_time> immunity_loss_queue;

double get_death_time(Host * host) { return host->death_time; }
EventQueue<Host, get_death_time> death_queue;

double get_transition_time(Infection * infection) { return infection->transition_time; }
EventQueue<Infection, get_transition_time> transition_queue;

double get_mutation_time(Infection * infection) { return infection->mutation_time; }
EventQueue<Infection, get_mutation_time> mutation_queue;

double get_recombination_time(Infection * infection) { return infection->recombination_time; }
EventQueue<Infection, get_recombination_time> recombination_queue;

double get_clearance_time(Infection * infection) { return infection->clearance_time; }
EventQueue<Infection, get_clearance_time> clearance_queue;

#pragma mark \
*** Helper function declarations ***

void initialize_gene_pool();

void initialize_populations();
void initialize_population(uint64_t order);
void initialize_population_events(Population * pop);
void initialize_population_hosts(Population * pop);
void initialize_population_infections(Population * pop);

void destroy_host(Host * host);
Host * create_host(Population * pop);
Gene * get_gene_with_alleles(std::array<uint64_t, N_LOCI> const & alleles);
Strain * generate_random_strain();
Strain * create_strain(std::vector<Gene *> const & genes);
Gene * draw_random_gene();
Strain * get_strain_with_genes(std::array<Gene *, N_GENES_PER_STRAIN> genes);

void destroy_infection(Infection * infection);

Gene * get_current_gene(Infection * infection);

uint64_t get_immune_allele_count(Host * host);

void gain_immunity(Host * host, Gene * gene);
void update_immunity_loss_time(Host * host);

void lose_random_immunity(Host * host);

void update_host_infection_times(Host * host);
void update_infection_times(Infection * infection);

double get_specific_immunity_level(Host * host, Gene * gene);
uint64_t get_active_infection_count(Host * host);
void infect_host(Host * host, Strain * strain);

void perform_infection_transition(Infection * infection);
void clear_infection(Infection * infection);

void update_transition_time(Infection * infection, bool initial);
double draw_activation_time(Infection * infection);
double draw_deactivation_time(Infection * infection);

void update_mutation_time(Infection * infection, bool initial);
void update_recombination_time(Infection * infection, bool initial);
void update_clearance_time(Infection * infection, bool initial);

void do_next_event(); 

void do_verification_event();
void do_checkpoint_event();

void do_biting_event();
void do_immigration_event();

void do_immunity_loss_event();
void do_death_event();

void do_transition_event();
void do_mutation_event();
void do_recombination_event();
void do_clearance_event(); 

void update_biting_time(Population * pop, bool initial);
void update_immigration_time(Population * pop, bool initial);

double draw_exponential_after_now(double lambda);
double draw_exponential(double lambda);
uint64_t draw_uniform_index(uint64_t size);
double draw_uniform_real(double min, double max);

#pragma mark \
*** Printing/debugging helpers ***

uint64_t same_time_count = 0;
uint64_t call_depth = 0;

#define ENTER() { \
    call_depth++; \
    if(PRINT_FUNCTION_TRACE) { \
        for(uint64_t i = 0; i < call_depth; i++) { \
            printf(" "); \
        } \
        printf("ENTER: %s\n", __func__); \
    } \
}

#define LEAVE() { \
    if(PRINT_FUNCTION_TRACE) { \
        for(uint64_t i = 0; i < call_depth; i++) { \
            printf(" "); \
        } \
        printf("LEAVE: %s\n", __func__); \
    } \
    call_depth--; \
}

#define PRINT_DEBUG(level, ...) { \
    if(level >= PRINT_DEBUG_LEVEL) { \
        for(uint64_t i = 0; i <= call_depth; i++) { \
            printf(" "); \
        } \
        printf(__VA_ARGS__); \
    } \
}

#pragma mark \
*** Initialization function implementations ***

void validate_and_load_parameters() {
    ENTER();
    
    assert(SAMPLE_DB_FILENAME == "" || !file_exists(SAMPLE_DB_FILENAME));
    assert(HOST_SAMPLING_PERIOD > 0.0);
    assert(!SAMPLE_TRANSMISSION_EVENTS || TRANSMISSION_EVENT_SAMPLING_SKIP > 0);
    assert(T_END >= 0.0);
    assert(T_BURNIN >= 0.0);
    assert(T_BURNIN < T_END);
    
    assert(!SAVE_TO_CHECKPOINT || !file_exists(CHECKPOINT_SAVE_FILENAME));
    assert(!SAVE_TO_CHECKPOINT || CHECKPOINT_SAVE_PERIOD > 0.0);
    
    assert(!LOAD_FROM_CHECKPOINT || file_exists(CHECKPOINT_LOAD_FILENAME));
    
    assert(VERIFICATION_PERIOD >= 0.0);
    
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
        assert(accumulate(HOST_LIFETIME_PDF) > 0.0);
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
    
    LEAVE();
}

void initialize() {
    ENTER();
    
    initialize_gene_pool();
    initialize_populations();
    
    LEAVE();
}

void initialize_gene_pool() {
    ENTER();
    
    // Create gene pool
    for(uint64_t i = 0; i < N_GENES_IN_POOL; i++) {
        // Draw random alleles different from all other existing genes
        std::array<uint64_t, N_LOCI> alleles;
        do {
            for(uint64_t j = 0; j < N_LOCI; j++) {
                alleles[j] = draw_uniform_index(N_ALLELES[j]);
            }
        } while(get_gene_with_alleles(alleles) != NULL);
        
        Gene * gene = gene_manager.create();
        gene->alleles = alleles;
    }
    
    LEAVE();
}

void initialize_populations() {
    ENTER();
    
    for(uint64_t i = 0; i < N_POPULATIONS; i++) {
        initialize_population(i);
    }
    
    LEAVE();
}

void initialize_population(uint64_t order) {
    ENTER();
    
    Population * pop = population_manager.create();
    pop->order = order;
    pop->transmission_count = 0;
    update_biting_time(pop, true);
    update_immigration_time(pop, true);
    
    initialize_population_hosts(pop);
    initialize_population_infections(pop);
    
    LEAVE();
}

void initialize_population_hosts(Population * pop) {
    ENTER();
    
    for(uint64_t i = 0; i < N_HOSTS[pop->order]; i++) {
        create_host(pop);
    }
    
    LEAVE();
}

void initialize_population_infections(Population * pop) {
    ENTER();
    
    for(uint64_t i = 0; i < N_INITIAL_INFECTIONS[pop->order]; i++) {
        Host * host = pop->hosts.object_at_index(draw_uniform_index(pop->hosts.size())); 
        Strain * strain = generate_random_strain();
        infect_host(host, strain);
    }
    
    LEAVE();
}


#pragma mark \
*** Object management functions ***

Host * create_host(Population * pop) {
    ENTER();
    
    Host * host = host_manager.create();
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
    host->birth_time = now;
    host->death_time = host->birth_time + lifetime;
    death_queue.add(host);
    
    host->completed_infection_count = 0;
    host->immunity_loss_time = INF;
    immunity_loss_queue.add(host);
    
    if(SELECTION_MODE == SPECIFIC_IMMUNITY) {
        AlleleImmuneHistory * immune_history = immune_history_manager.create();
        for(uint64_t i = 0; i < N_LOCI; i++) {
            immune_history->immunity_by_locus[i] = locus_immunity_manager.create();
        }
        host->immune_history = immune_history;
    }
    else {
        host->immune_history = NULL;
    }
    
    pop->hosts.add(host);
    
    LEAVE();
    
    return host;
}

void destroy_host(Host * host) {
    ENTER();
    
    Population * pop = host->population;
    PRINT_DEBUG(5, "Removing host id %llu from population %llu\n", host->id, pop->id);
    
    for(Infection * infection : host->infections) {
        destroy_infection(infection);
    }
    
    if(SELECTION_MODE == SPECIFIC_IMMUNITY) {
        for(LocusImmunity * immunity : host->immune_history->immunity_by_locus) {
            locus_immunity_manager.destroy(immunity);
        }
        immune_history_manager.destroy(host->immune_history);
    }
    
    pop->hosts.remove(host);
    death_queue.remove(host);
    immunity_loss_queue.remove(host);
    host_manager.destroy(host);
    
    LEAVE();
}

void destroy_infection(Infection * infection) {
    ENTER();
    transition_queue.remove(infection);
    mutation_queue.remove(infection);
    recombination_queue.remove(infection);
    clearance_queue.remove(infection);
    infection_manager.destroy(infection);
    LEAVE();
}

Gene * get_gene_with_alleles(std::array<uint64_t, N_LOCI> const & alleles) {
    ENTER();
    
    for(Gene * gene : gene_manager.objects()) {
        assert(gene->alleles.size() == N_LOCI);
        if(gene->alleles == alleles) {
            LEAVE();
            return gene;
        }
    }
    
    LEAVE();
    return nullptr;
}

Strain * generate_random_strain() {
    ENTER();
    
    std::array<Gene *, N_GENES_PER_STRAIN> genes;
    for(uint64_t i = 0; i < N_GENES_PER_STRAIN; i++) {
        genes[i] = draw_random_gene();
    }
    
    Strain * strain = get_strain_with_genes(genes);
    LEAVE();
    return strain;
}

Strain * get_strain_with_genes(std::array<Gene *, N_GENES_PER_STRAIN> genes) {
    ENTER();
    Strain * strain;
    std::sort(genes.begin(), genes.end(),
        [](Gene * o1, Gene * o2) {
            return o1->id < o2->id;
        }
    );
    auto itr = genes_strain_map.find(genes);
    if(itr == genes_strain_map.end()) {
        strain = strain_manager.create();
        strain->genes_sorted = genes;
        genes_strain_map[genes] = strain;
    }
    else {
        strain = itr->second;
    }
    LEAVE();
    return strain;
}

Gene * draw_random_gene() {
    ENTER();
    Gene * gene = gene_manager.objects()[draw_uniform_index(gene_manager.size())];
    LEAVE();
    return gene;
}

Gene * get_current_gene(Infection * infection) {
    ENTER();
    Gene * gene = infection->expression_order[infection->expression_index];
    LEAVE();
    return gene;
}

void gain_immunity(Host * host, Gene * gene) {
    ENTER();
    AlleleImmuneHistory * immune_history = host->immune_history;
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
    update_immunity_loss_time(host);
    LEAVE();
}

void update_immunity_loss_time(Host * host) {
    ENTER();
    uint64_t immune_allele_count = get_immune_allele_count(host);
    host->immunity_loss_time = draw_exponential_after_now(IMMUNITY_LOSS_RATE * immune_allele_count);
    immunity_loss_queue.update(host);
    LEAVE();
}

void gain_general_immunity(Host * host, Gene * gene) {
    ENTER();
    LEAVE();
}

void lose_random_immunity(Host * host) {
    ENTER();
    
    AlleleImmuneHistory * immune_history = host->immune_history;
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
    update_immunity_loss_time(host);
    LEAVE();
}

uint64_t get_immune_allele_count(Host * host) {
    ENTER();
    uint64_t immune_allele_count = 0;
    AlleleImmuneHistory * immune_history = host->immune_history;
    for(uint64_t i = 0; i < N_LOCI; i++) {
        LocusImmunity * immunity = immune_history->immunity_by_locus[i];
        immune_allele_count += immunity->immunity_level_by_allele.size();
    }
    LEAVE();
    return immune_allele_count;
}

double get_specific_immunity_level(Host * host, Gene * gene) {
    ENTER();
    
    uint64_t immunity_count = 0;
    for(uint64_t i = 0; i < N_LOCI; i++) {
        auto immunity_level_by_allele = host->immune_history->immunity_by_locus[i]->immunity_level_by_allele; 
        auto itr = immunity_level_by_allele.find(gene->alleles[i]);
        if(itr != immunity_level_by_allele.end()) {
            immunity_count += std::min(itr->second, 1ULL);
        }
    }
    
    LEAVE();
    return immunity_count / (double)N_LOCI;
}

uint64_t get_active_infection_count(Host * host) {
    ENTER();
    uint64_t count = 0;
    for(Infection * infection : host->infections) {
        if(infection->active) {
            count += 1;
        }
    }
    
    LEAVE();
    return count;
}

void infect_host(Host * host, Strain * strain) {
    ENTER();
    
    Infection * infection = infection_manager.create();
    infection->host = host;
    infection->strain = strain;
    
    for(uint64_t i = 0; i < strain->genes_sorted.size(); i++) {
        infection->expression_order[i] = strain->genes_sorted[i];
    }
    std::shuffle(
        infection->expression_order.begin(), infection->expression_order.end(), rng
    );
    
    infection->active = false;
    if(T_LIVER_STAGE == 0.0) {
        infection->expression_index = 0;
        update_transition_time(infection, true);
    }
    else {
        infection->expression_index = -1;
        infection->transition_time = now + T_LIVER_STAGE;
        transition_queue.add(infection);
    }
    
    update_mutation_time(infection, true);
    update_recombination_time(infection, true);
    
    if(SELECTION_MODE == GENERAL_IMMUNITY) {
        update_clearance_time(infection, true);
    }
    
    host->infections.insert(infection);
    
    LEAVE();
}

void perform_infection_transition(Infection * infection) {
    ENTER();
    
    Host * host = infection->host;
    if(infection->expression_index == -1) {
        assert(!infection->active);
        infection->expression_index = 0;
    }
    else if(infection->active) {
        assert(
            infection->expression_index >= 0 &&
            infection->expression_index < N_GENES_PER_STRAIN
        );
        Gene * gene = get_current_gene(infection);
        
        if(SELECTION_MODE == SPECIFIC_IMMUNITY) {
            gain_immunity(host, gene);
        }
        
        infection->expression_index++;
        infection->active = false;
    }
    else {
        infection->active = true;
    }
    
    if(infection->expression_index == N_GENES_PER_STRAIN) {
        clear_infection(infection);
    }
    else {
        update_host_infection_times(infection->host);
    }
    
    LEAVE();
}

void clear_infection(Infection * infection) {
    ENTER();
    
    Host * host = infection->host;
    host->infections.erase(infection);
    
    destroy_infection(infection);
    host->completed_infection_count++;
    
    if(SELECTION_MODE == GENERAL_IMMUNITY) {
        update_host_infection_times(host);
    }
    
    LEAVE();
}

void update_host_infection_times(Host * host) {
    ENTER();
    for(Infection * infection : host->infections) {
        update_infection_times(infection);
    }
    LEAVE();
}

void update_infection_times(Infection * infection) {
    ENTER();
    if(infection->expression_index >= 0) {
        update_transition_time(infection, false);
    }
    LEAVE();
}

void update_transition_time(Infection * infection, bool initial) {
    ENTER();
    assert(infection->expression_index > -1 && infection->expression_index < N_GENES_PER_STRAIN);
    if(infection->active) {
        infection->transition_time = draw_deactivation_time(infection);
    }
    else {
        infection->transition_time = draw_activation_time(infection);
    }
    if(initial) {
        transition_queue.add(infection);
    }
    else {
        transition_queue.update(infection);
    }
    LEAVE();
}

double draw_activation_time(Infection * infection) {
    ENTER();
    double time;
    if(get_active_infection_count(infection->host) == 0) {
        time = now;
    }
    else {
        time = draw_exponential_after_now(ACTIVATION_RATE);
    }
    LEAVE();
    return time;
}

double draw_deactivation_time(Infection * infection) {
    ENTER();
    assert(infection->active);
    Host * host = infection->host;
    double rate;
    if(SELECTION_MODE == SPECIFIC_IMMUNITY) {
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
    double time = draw_exponential_after_now(rate); 
    LEAVE();
    return time;
}

void update_mutation_time(Infection * infection, bool initial) {
    ENTER();
    infection->mutation_time = draw_exponential_after_now(0.0);
    if(initial) {
        mutation_queue.add(infection);
    }
    else {
        mutation_queue.update(infection);
    }
    LEAVE();
}

void update_recombination_time(Infection * infection, bool initial) {
    ENTER();
    infection->recombination_time = draw_exponential_after_now(0.0);
    if(initial) {
        recombination_queue.add(infection);
    }
    else {
        recombination_queue.update(infection);
    }
    LEAVE();
}

void update_clearance_time(Infection * infection, bool initial) {
    ENTER();
    
    assert(SELECTION_MODE == GENERAL_IMMUNITY);
    
    // TODO: replace with Qixin's function based on vector of parameters
    double rate = 0.0;
    if(infection->active) {
        rate = CLEARANCE_RATE;
    }
    
    infection->clearance_time = draw_exponential_after_now(rate);
    if(initial) {
        clearance_queue.add(infection);
    }
    else {
        clearance_queue.update(infection);
    }
    LEAVE();
}


#pragma mark \
*** Main simulation loop ***

enum class EventType {
    NONE,
    VERIFICATION,
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
    ENTER();
    while(now < T_END) { 
        do_next_event();
    }
    LEAVE();
}

void do_next_event() {
    ENTER();
    double next_event_time = INF;
    EventType next_event_type = EventType::NONE;
    
    // Find the event type with the lowest-timed event (simple linear search)
    // This is the least elegant part of the reimplementation, but it should
    // be worth it for simplicity and reduced memory usage.
    // If this section is a speed bottleneck, re-evaluate the choice.
    // It shouldn't come even close.
    if(VERIFICATION_ON && next_verification_time < next_event_time) {
        next_event_time = next_verification_time;
        next_event_type = EventType::VERIFICATION;
    }
    if(SAVE_TO_CHECKPOINT && next_checkpoint_time < next_event_time) {
        next_event_time = next_checkpoint_time;
        next_event_type = EventType::CHECKPOINT;
    }
    if(biting_queue.next_time() < next_event_time) {
        next_event_time = biting_queue.next_time();
        next_event_type = EventType::BITING; 
    }
    if(immigration_queue.next_time() < next_event_time) {
        next_event_time = immigration_queue.next_time();
        next_event_type = EventType::IMMIGRATION;
    }
    if(immunity_loss_queue.next_time() < next_event_time) {
        next_event_time = immunity_loss_queue.next_time();
        next_event_type = EventType::IMMUNITY_LOSS; 
    }
    if(death_queue.next_time() < next_event_time) {
        next_event_time = death_queue.next_time();
        next_event_type = EventType::DEATH;
    }
    if(transition_queue.next_time() < next_event_time) {
        next_event_time = transition_queue.next_time();
        next_event_type = EventType::TRANSITION;
    }
    if(mutation_queue.next_time() < next_event_time) {
        next_event_time = mutation_queue.next_time();
        next_event_type = EventType::MUTATION;
    }
    if(recombination_queue.next_time() < next_event_time) {
        next_event_time = recombination_queue.next_time();
        next_event_type = EventType::RECOMBINATION;
    }
    if(
        SELECTION_MODE == GENERAL_IMMUNITY &&
        clearance_queue.next_time() < next_event_time
    ) {
        next_event_time = clearance_queue.head()->clearance_time;
        next_event_type = EventType::CLEARANCE;
    }
    
    PRINT_DEBUG(0, "next_event_time: %f\n", next_event_time);
    PRINT_DEBUG(0, "next_event_type: %d\n", next_event_type);
    
    // Execute the next event
    assert(next_event_time >= now);
    
    // If we get the same time a whole bunch of times in a row, we've introduced
    // a bug where the same event is being reused over and over.
    // We make the limit 2 * sum(N_INITIAL_INFECTIONS) to allow a bunch of events at
    // t = T_LIVER_STAGE to be OK (along with an instantaneous activation)
    if(next_event_time == now) {
        same_time_count++;
    }
    else {
        same_time_count = 0;
    }
    assert(same_time_count <= 2 * accumulate(N_INITIAL_INFECTIONS));
    
    now = next_event_time;
    switch(next_event_type) {
        case EventType::NONE:
            break;
        case EventType::VERIFICATION:
            do_verification_event();
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
    LEAVE();
}

void do_verification_event() {
    ENTER();
    verify_simulation_state();
    next_verification_time += VERIFICATION_PERIOD;
    LEAVE();
}

void do_checkpoint_event() {
    ENTER();
    save_checkpoint();
    next_checkpoint_time += CHECKPOINT_SAVE_PERIOD;
    LEAVE();
}

void do_biting_event() {
    ENTER();
    assert(biting_queue.size() > 0);
    Population * pop = biting_queue.head();
    
    PRINT_DEBUG(1, "biting event pop: %llu\n", pop->id);
    
    // TODO: biting activity
    
    // Update biting event time
    update_biting_time(pop, false);
    LEAVE();
}

void do_immigration_event() {
    ENTER();
    assert(immigration_queue.size() > 0);
    Population * pop = immigration_queue.head();
    
    PRINT_DEBUG(1, "immigration event pop: %llu\n", pop->id);
    
    // TODO: immigration activity
    
    // Update immigration event time
    update_immigration_time(pop, false);
    LEAVE();
}

void do_immunity_loss_event() {
    ENTER();
    assert(immunity_loss_queue.size() > 0);
    
    Host * host = immunity_loss_queue.head();
    
    switch(SELECTION_MODE) {
        case SPECIFIC_IMMUNITY:
            lose_random_immunity(host);
            break;
        case GENERAL_IMMUNITY:
        case NEUTRALITY:
            assert(false);
            break;
    }
    LEAVE();
}

void do_death_event() {
    ENTER();
    assert(death_queue.size() > 0);
    
    Host * host = death_queue.head();
    Population * pop = host->population;
    destroy_host(host);
    Host * new_host = create_host(pop);
    PRINT_DEBUG(1, "Created new host id %llu in population %llu\n", new_host->id, pop->id);
    LEAVE();
}

void do_transition_event() {
    ENTER();
    assert(transition_queue.size() > 0);
    
    Infection * infection = transition_queue.head();
    perform_infection_transition(infection);
    LEAVE();
}

void do_mutation_event() {
    ENTER();
    LEAVE();
}

void do_recombination_event() {
    ENTER();
    LEAVE();
}

void do_clearance_event() {
    ENTER();
    assert(SELECTION_MODE == GENERAL_IMMUNITY);
    Infection * infection = clearance_queue.head();
    clear_infection(infection);
    LEAVE();
}

#pragma mark \
*** Verification function implementations ***

void verify_simulation_state() {
    ENTER();
    
    assert(biting_queue.verify_heap());
    assert(immigration_queue.verify_heap());
    assert(immunity_loss_queue.verify_heap());
    assert(death_queue.verify_heap());
    assert(transition_queue.verify_heap());
    assert(mutation_queue.verify_heap());
    assert(recombination_queue.verify_heap());
    assert(clearance_queue.verify_heap());
    
    LEAVE();
}

#pragma mark \
*** Checkpoint function implementations ***

void save_checkpoint() {
    ENTER();
    std::string old_checkpoint_filename = CHECKPOINT_SAVE_FILENAME + "-old";
    if(file_exists(CHECKPOINT_SAVE_FILENAME)) {
        assert(!rename(CHECKPOINT_SAVE_FILENAME.c_str(), old_checkpoint_filename.c_str()));
    }
    
    sqlite3 * db;
    assert(!file_exists(CHECKPOINT_SAVE_FILENAME));
    sqlite3_open(CHECKPOINT_SAVE_FILENAME.c_str(), &db);
    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, NULL);
    
    strain_manager.save_to_checkpoint(db);
    gene_manager.save_to_checkpoint(db);
    population_manager.save_to_checkpoint(db);
    host_manager.save_to_checkpoint(db);
    infection_manager.save_to_checkpoint(db);
    immune_history_manager.save_to_checkpoint(db);
    locus_immunity_manager.save_to_checkpoint(db);
    
    sqlite3_exec(db, "END TRANSACTION", NULL, NULL, NULL);
    int result = sqlite3_close(db);
    assert(result == SQLITE_OK);
    
    if(file_exists(old_checkpoint_filename)) {
        assert(!unlink(old_checkpoint_filename.c_str()));
    }
    LEAVE();
}

void load_checkpoint() {
    ENTER();
    LEAVE();
}


#pragma mark \
*** Random draw helper function implementations *** 

double draw_exponential_after_now(double lambda) {
    ENTER();
    double time;
    if(lambda == 0.0) {
        time = INF;
    }
    else {
        time = now + draw_exponential(lambda);
    }
    LEAVE();
    return time;
}

double draw_exponential(double lambda) {
    ENTER();
    LEAVE();
    return std::exponential_distribution<>(lambda)(rng); 
}

uint64_t draw_uniform_index(uint64_t size) {
    ENTER();
    LEAVE();
    return std::uniform_int_distribution<uint64_t>(0, size - 1)(rng);
}

double draw_uniform_real(double min, double max) {
    ENTER();
    LEAVE();
    return std::uniform_real_distribution<double>(min, max)(rng);
}

void update_biting_time(Population * pop, bool initial) {
    ENTER();
    double biting_rate = BITING_RATE_MEAN[pop->order] * (
        1.0 + BITING_RATE_RELATIVE_AMPLITUDE[pop->order] * cos(
            2 * M_PI * ((now / T_YEAR) - BITING_RATE_PEAK_PHASE[pop->order])
        )
    );
    pop->next_biting_time = draw_exponential_after_now(biting_rate);
    if(initial) {
        biting_queue.add(pop);
    }
    else {
        biting_queue.update(pop);
    }
    LEAVE();
}

void update_immigration_time(Population * pop, bool initial) {
    ENTER();
    pop->next_immigration_time = draw_exponential_after_now(IMMIGRATION_RATE[pop->order]);
    if(initial) {
        immigration_queue.add(pop);
    }
    else {
        immigration_queue.update(pop);
    }
    LEAVE();
}

} // namespace varmodel

