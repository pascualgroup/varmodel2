#include "varmodel.hpp"

#include "StrainManager.hpp"
#include "GeneManager.hpp"
#include "PopulationManager.hpp"
#include "HostManager.hpp"
#include "InfectionManager.hpp"
#include "ImmuneHistoryManager.hpp"
#include "LocusImmunityManager.hpp"
#include "random.hpp"
#include "util.hpp"
#include "parameters.hpp"
#include "EventQueue.hpp"

#include <unistd.h>
#include <sstream>
#include <algorithm>
#include <bitset>

namespace varmodel {

#pragma mark \
*** SIMULATION STATE ***

double INF = std::numeric_limits<double>::infinity();

std::mt19937_64 rng(RANDOM_SEED);

double now = 0.0;
double next_sampling_time = T_BURNIN;
double next_verification_time = VERIFICATION_ON ? 0.0 : INF;
double next_checkpoint_time = SAVE_TO_CHECKPOINT ? 0.0 : INF;
double next_info_time = 0.0;

uint64_t n_infections_cumulative = 0;

std::array<uint64_t, N_LOCI> n_alleles = N_ALLELES_INITIAL;

StrainManager strain_manager;
std::unordered_map<
    std::array<Gene *, N_GENES_PER_STRAIN>,
    Strain *,
    HashArray<Gene *, N_GENES_PER_STRAIN>
> genes_strain_map;

GeneManager gene_manager;
std::unordered_map<
    std::array<uint64_t, N_LOCI>,
    Gene *,
    HashArray<uint64_t, N_LOCI>
> alleles_genes_map;

PopulationManager population_manager;
HostManager host_manager;
InfectionManager infection_manager;
ImmuneHistoryManager immune_history_manager;
LocusImmunityManager locus_immunity_manager;

sqlite3 * sample_db;

sqlite3_stmt * summary_stmt;
sqlite3_stmt * summary_alleles_stmt;

sqlite3_stmt * host_stmt;
sqlite3_stmt * strain_stmt;
sqlite3_stmt * gene_stmt;
sqlite3_stmt * allele_stmt;

sqlite3_stmt * sampled_host_stmt;
sqlite3_stmt * sampled_inf_stmt;
sqlite3_stmt * sampled_strain_stmt;
sqlite3_stmt * sampled_gene_stmt;
sqlite3_stmt * sampled_allele_stmt;

#pragma mark \
*** EVENT QUEUES ***

double get_biting_time(Population * p) { return p->next_biting_time; }
EventQueue<Population, get_biting_time> biting_queue; 

double get_immigration_time(Population * p) { return p->next_immigration_time; }
EventQueue<Population, get_immigration_time> immigration_queue;

double get_next_immunity_loss_time(Host * h) { return h->next_immunity_loss_time; }
EventQueue<Host, get_next_immunity_loss_time> immunity_loss_queue;

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

void initialize(bool override_seed, uint64_t random_seed);
void clean_up();

void validate_and_load_parameters();

void write_summary();
void sample_hosts();
void write_host(Host * host);
void write_strain(Strain * strain, sqlite3_stmt * s_stmt, sqlite3_stmt * g_stmt, sqlite3_stmt * a_stmt);
void write_gene(Gene * gene, sqlite3_stmt * g_stmt, sqlite3_stmt * a_stmt);
void write_sampled_host(Host * host);
void write_sampled_infection(Host * host, Infection * infection);

void load_checkpoint(bool should_load_rng_state);
void save_checkpoint();

void save_global_state_to_checkpoint(sqlite3 * db);
void load_global_state_from_checkpoint(sqlite3 * db, bool should_load_rng_state);
void initialize_event_queues_from_state();

std::string get_rng_as_string();
void set_rng_from_string(std::string const & rng_str);

void initialize_sample_db();
void finalize_sample_db();

void initialize_gene_pool();

void initialize_populations();
void initialize_population(uint64_t order);
void initialize_population_events(Population * pop);
void initialize_population_hosts(Population * pop);
void initialize_population_infections(Population * pop);

Host * create_host(Population * pop);
void destroy_host(Host * host);

Gene * get_gene_with_alleles(std::array<uint64_t, N_LOCI> const & alleles);
Gene * get_or_create_gene(std::array<uint64_t, N_LOCI> const & alleles, GeneSource source, bool is_functional);
Strain * generate_random_strain(uint64_t n_new_genes, GeneSource new_gene_source);
Strain * create_strain(std::vector<Gene *> const & genes);
Gene * draw_random_gene();
std::array<uint64_t, N_LOCI> recombine_alleles(std::array<uint64_t, N_LOCI> const & a1, std::array<uint64_t, N_LOCI> const & a2, uint64_t breakpoint);
bool contains_different_genes(Strain * strain);
Strain * get_strain_with_genes(std::array<Gene *, N_GENES_PER_STRAIN> genes);

Strain * recombine_strains(Strain * s1, Strain * s2);
Strain * mutate_strain(Strain * strain);
Gene * mutate_gene(Gene * gene, GeneSource source);
void recombine_infection(Infection * infection);
double get_gene_similarity(Gene * gene1, Gene * gene2, uint64_t breakpoint);

void destroy_infection(Infection * infection);

Gene * get_current_gene(Infection * infection);

uint64_t get_immune_allele_count(Host * host);

void gain_immunity(Host * host, Gene * gene);
void update_next_immunity_loss_time(Host * host);

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

bool do_next_event(); 

void do_verification_event();
void do_sampling_event();
void do_checkpoint_event();

void do_biting_event();
Host * draw_random_source_host(Population * pop);
Host * draw_random_destination_host(Population * src_pop);
void transmit(Host * src_host, Host * dst_host);
double get_transmission_probability(Infection * infection);

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
uint64_t draw_uniform_index_except(uint64_t size, uint64_t except_index);
double draw_uniform_real(double min, double max);
std::vector<uint64_t> draw_uniform_indices_without_replacement(uint64_t n, uint64_t k);
bool draw_bernoulli(double p);

#pragma mark \
*** Printing/debugging helpers ***

uint64_t same_time_count = 0;
uint64_t call_depth = 0;

#define BEGIN() { \
    call_depth++; \
    if(PRINT_FUNCTION_TRACE) { \
        for(uint64_t i = 0; i < call_depth; i++) { \
            printf(" "); \
        } \
        printf("BEGIN: %s\n", __func__); \
    } \
}

#define RETURN(x) { \
    if(PRINT_FUNCTION_TRACE) { \
        for(uint64_t i = 0; i < call_depth; i++) { \
            printf(" "); \
        } \
        printf("RETURN: %s\n", __func__); \
    } \
    call_depth--; \
    return x; \
}

#define PRINT_DEBUG(level, ...) { \
    if(level <= PRINT_DEBUG_LEVEL) { \
        for(uint64_t i = 0; i <= call_depth; i++) { \
            printf(" "); \
        } \
        printf(__VA_ARGS__); \
        printf("\n"); \
    } \
}

#pragma mark \
*** Initialization function implementations ***

void validate_and_load_parameters() {
    BEGIN();
    
    assert(RANDOM_SEED > 0);
    
    assert(T_YEAR > 0.0);
    assert(T_END >= 0.0);
    assert(T_BURNIN >= 0.0);
    assert(T_BURNIN <= T_END);
    
    assert(SAMPLE_DB_FILENAME == "" || !file_exists(SAMPLE_DB_FILENAME));
    assert(HOST_SAMPLING_PERIOD > 0.0);
    
    assert(!SAVE_TO_CHECKPOINT || !file_exists(CHECKPOINT_SAVE_FILENAME));
    assert(!SAVE_TO_CHECKPOINT || CHECKPOINT_SAVE_PERIOD > 0.0);
    
    assert(!LOAD_FROM_CHECKPOINT || file_exists(CHECKPOINT_LOAD_FILENAME));
    
    assert(!VERIFICATION_ON || VERIFICATION_PERIOD > 0.0);
    
    assert(N_GENES_INITIAL >= 1);
    assert(N_GENES_PER_STRAIN >= 2);
    assert(N_LOCI >= 1);
    assert(N_ALLELES_INITIAL.size() == N_LOCI);
    for(auto value : N_ALLELES_INITIAL) {
        assert(value >= 1);
    }
    
    assert(GENE_TRANSMISSIBILITY >= 0.0 && GENE_TRANSMISSIBILITY <= 1.0);
    assert(IMMUNITY_LOSS_RATE >= 0.0);
    
    assert(MUTATION_RATE >= 0.0);
    assert(T_LIVER_STAGE >= 0.0);
    
    assert(P_ECTOPIC_RECOMBINATION_IS_CONVERSION >= 0.0 && P_ECTOPIC_RECOMBINATION_IS_CONVERSION <= 1.0);
    
    assert(MEAN_HOST_LIFETIME > 0.0);
    assert(MAX_HOST_LIFETIME > 0.0);
    
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
    
    assert(!IMMIGRATION_ON || IMMIGRATION_RATE.size() == N_POPULATIONS);
    if(IMMIGRATION_ON) {
        for(auto value : IMMIGRATION_RATE) {
            assert(value >= 0.0);
        }
        assert(P_IMMIGRATION_INCLUDES_NEW_GENES >= 0.0 && P_IMMIGRATION_INCLUDES_NEW_GENES <= 1.0);
        assert(N_IMMIGRATION_NEW_GENES <= N_GENES_PER_STRAIN);
    }
    
    if(SELECTION_MODE == GENERAL_IMMUNITY) {
        assert(GENERAL_IMMUNITY_PARAMS[0] > 0.0);
        assert(GENERAL_IMMUNITY_PARAMS[1] > 0.0);
        assert(GENERAL_IMMUNITY_PARAMS[2] > 0.0);
        assert(GENERAL_IMMUNITY_PARAMS[3] > 0.0);
    }
    
    RETURN();
}

void initialize(bool override_seed, uint64_t random_seed) {
    BEGIN();
    
    validate_and_load_parameters();
    initialize_sample_db();
    
    if(override_seed) {
        rng.seed(random_seed);
    }
    
    if(LOAD_FROM_CHECKPOINT) {
        load_checkpoint(!override_seed);
    }
    else {
        initialize_gene_pool();
        initialize_populations();
    }
    
    RETURN();
}

void clean_up() {
    BEGIN();
    
    finalize_sample_db();
    
    RETURN();
}

void initialize_gene_pool() {
    BEGIN();
    
    // Create gene pool
    for(uint64_t i = 0; i < N_GENES_INITIAL; i++) {
        // Draw random alleles different from all other existing genes
        std::array<uint64_t, N_LOCI> alleles;
        do {
            for(uint64_t j = 0; j < N_LOCI; j++) {
                alleles[j] = draw_uniform_index(n_alleles[j]);
            }
        } while(get_gene_with_alleles(alleles) != NULL);
        get_or_create_gene(alleles, SOURCE_POOL, true);
    }
    
    RETURN();
}

void initialize_populations() {
    BEGIN();
    
    for(uint64_t i = 0; i < N_POPULATIONS; i++) {
        initialize_population(i);
    }
    
    RETURN();
}

void initialize_population(uint64_t index) {
    BEGIN();
    
    Population * pop = population_manager.create();
    pop->ind = index;
    pop->transmission_count = 0;
    update_biting_time(pop, true);
    
    if(IMMIGRATION_ON) {
        update_immigration_time(pop, true);
    }
    
    initialize_population_hosts(pop);
    initialize_population_infections(pop);
    
    RETURN();
}

void initialize_population_hosts(Population * pop) {
    BEGIN();
    
    for(uint64_t i = 0; i < N_HOSTS[pop->ind]; i++) {
        create_host(pop);
    }
    
    RETURN();
}

void initialize_population_infections(Population * pop) {
    BEGIN();
    
    for(uint64_t i = 0; i < N_INITIAL_INFECTIONS[pop->ind]; i++) {
        Host * host = pop->hosts.object_at_index(draw_uniform_index(pop->hosts.size())); 
        Strain * strain = generate_random_strain(0, SOURCE_POOL);
        infect_host(host, strain);
    }
    
    RETURN();
}

#pragma mark \
*** Sample database output ***

void initialize_sample_db() {
    sqlite3_open(SAMPLE_DB_FILENAME.c_str(), &sample_db);
    
    sqlite3_exec(sample_db, "BEGIN TRANSACTION", NULL, NULL, NULL);
    
    sqlite3_exec(sample_db,
        "CREATE TABLE IF NOT EXISTS hosts (id INTEGER PRIMARY KEY, population_id INTEGER, birth_time REAL, death_time REAL);",
        NULL, NULL, NULL
    );
    
    sqlite3_exec(sample_db,
        "CREATE TABLE IF NOT EXISTS strains (id INTEGER, ind INTEGER, gene_id INTEGER, PRIMARY KEY (id, ind));",
        NULL, NULL, NULL
    );
    
    sqlite3_exec(sample_db,
        "CREATE TABLE IF NOT EXISTS genes (id INTEGER PRIMARY KEY, source INTEGER, is_functional INTEGER);",
        NULL, NULL, NULL
    );
    
    sqlite3_exec(sample_db,
        "CREATE TABLE IF NOT EXISTS alleles (gene_id INTEGER, locus INTEGER, allele INTEGER, PRIMARY KEY (gene_id, locus));",
        NULL, NULL, NULL
    );
    
    sqlite3_exec(sample_db,
        "CREATE TABLE IF NOT EXISTS summary ("
            "time REAL, n_infections INTEGER, n_infected INTEGER, n_infections_cumulative INTEGER, n_circulating_strains INTEGER, n_circulating_genes INTEGER"
        ");",
        NULL, NULL, NULL
    );
    
    sqlite3_exec(sample_db,
        "CREATE TABLE IF NOT EXISTS summary_alleles ("
            "time REAL, locus INTEGER, n_circulating_alleles INTEGER"
        ");",
        NULL, NULL, NULL
    );
    
    sqlite3_exec(sample_db,
        "CREATE TABLE IF NOT EXISTS sampled_hosts ("
        "time REAL, id INTEGER, population_id INTEGER, birth_time REAL, death_time REAL, "
        "PRIMARY KEY (time, id)"
        ");",
        NULL, NULL, NULL
    );
    
    sqlite3_exec(sample_db,
        "CREATE TABLE IF NOT EXISTS sampled_infections (time REAL, host_id INTEGER, infection_id INTEGER, strain_id INTEGER, PRIMARY KEY (time, host_id, infection_id));",
        NULL, NULL, NULL
    );
    
    sqlite3_exec(sample_db,
        "CREATE TABLE IF NOT EXISTS sampled_strains (id INTEGER, ind INTEGER, gene_id INTEGER, PRIMARY KEY (id, ind));",
        NULL, NULL, NULL
    );
    
    sqlite3_exec(sample_db,
        "CREATE TABLE IF NOT EXISTS sampled_genes (id INTEGER PRIMARY KEY, source INTEGER, is_functional INTEGER);",
        NULL, NULL, NULL
    );
    
    sqlite3_exec(sample_db,
        "CREATE TABLE IF NOT EXISTS sampled_alleles (gene_id INTEGER, locus INTEGER, allele INTEGER, PRIMARY KEY (gene_id, locus));",
        NULL, NULL, NULL
    );
    
    sqlite3_prepare_v2(sample_db, "INSERT INTO summary VALUES (?,?,?,?,?,?);", -1, &summary_stmt, NULL);
    sqlite3_prepare_v2(sample_db, "INSERT INTO summary_alleles VALUES (?,?,?);", -1, &summary_alleles_stmt, NULL);
    
    sqlite3_prepare_v2(sample_db, "INSERT INTO hosts VALUES (?,?,?,?);", -1, &host_stmt, NULL);
    sqlite3_prepare_v2(sample_db, "INSERT INTO strains VALUES (?,?,?);", -1, &strain_stmt, NULL);
    sqlite3_prepare_v2(sample_db, "INSERT INTO genes VALUES (?,?,?);", -1, &gene_stmt, NULL);
    sqlite3_prepare_v2(sample_db, "INSERT INTO alleles VALUES (?,?,?);", -1, &allele_stmt, NULL);
    sqlite3_prepare_v2(sample_db, "INSERT INTO sampled_hosts VALUES (?,?,?,?,?);", -1, &sampled_host_stmt, NULL);
    sqlite3_prepare_v2(sample_db, "INSERT INTO sampled_infections VALUES (?,?,?,?);", -1, &sampled_inf_stmt, NULL);
    sqlite3_prepare_v2(sample_db, "INSERT OR IGNORE INTO sampled_strains VALUES (?,?,?);", -1, &sampled_strain_stmt, NULL);
    sqlite3_prepare_v2(sample_db, "INSERT OR IGNORE INTO sampled_genes VALUES (?,?,?);", -1, &sampled_gene_stmt, NULL);
    sqlite3_prepare_v2(sample_db, "INSERT OR IGNORE INTO sampled_alleles VALUES (?,?,?);", -1, &sampled_allele_stmt, NULL);
    
    sqlite3_exec(sample_db, "COMMIT", NULL, NULL, NULL);
    sqlite3_exec(sample_db, "BEGIN TRANSACTION", NULL, NULL, NULL);
}

void finalize_sample_db() {
    sqlite3_exec(sample_db, "COMMIT", NULL, NULL, NULL);
    
    sqlite3_finalize(summary_stmt);
    sqlite3_finalize(summary_alleles_stmt);
    
    sqlite3_finalize(host_stmt);
    sqlite3_finalize(sampled_host_stmt);
    sqlite3_finalize(strain_stmt);
    sqlite3_finalize(gene_stmt);
    sqlite3_finalize(sampled_gene_stmt);
    sqlite3_finalize(sampled_strain_stmt);
    sqlite3_finalize(allele_stmt);
    sqlite3_finalize(sampled_allele_stmt);
    sqlite3_finalize(sampled_inf_stmt);
    
    sqlite3_close(sample_db);
}

#pragma mark \
*** Object management functions ***

Host * create_host(Population * pop) {
    BEGIN();
    
    Host * host = host_manager.create();
    host->population = pop;
    
    double lifetime = std::min(
        draw_exponential(1.0 / MEAN_HOST_LIFETIME),
        MAX_HOST_LIFETIME
    );
    host->birth_time = now;
    host->death_time = host->birth_time + lifetime;
    death_queue.add(host);
    
    host->completed_infection_count = 0;
    host->next_immunity_loss_time = INF;
    immunity_loss_queue.add(host);
    
    if(SELECTION_MODE == SPECIFIC_IMMUNITY) {
        ImmuneHistory * immune_history = immune_history_manager.create();
        for(uint64_t i = 0; i < N_LOCI; i++) {
            immune_history->immunity_by_locus[i] = locus_immunity_manager.create();
        }
        host->immune_history = immune_history;
    }
    else {
        host->immune_history = NULL;
    }
    
    pop->hosts.add(host);
    
    if(OUTPUT_HOSTS) {
        write_host(host);
    }
    
    RETURN(host);
}

void destroy_host(Host * host) {
    BEGIN();
    
    Population * pop = host->population;
    PRINT_DEBUG(5, "Removing host id %llu from population %llu", host->id, pop->id);
    
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
    
    RETURN();
}

void destroy_infection(Infection * infection) {
    BEGIN();
    transition_queue.remove(infection);
    mutation_queue.remove(infection);
    recombination_queue.remove(infection);
    clearance_queue.remove(infection);
    infection_manager.destroy(infection);
    RETURN();
}

Gene * get_gene_with_alleles(std::array<uint64_t, N_LOCI> const & alleles) {
    BEGIN();
    
    auto itr = alleles_genes_map.find(alleles);
    if(itr == alleles_genes_map.end()) {
        RETURN(nullptr);
    }
    RETURN(itr->second);
}

Gene * get_or_create_gene(std::array<uint64_t, N_LOCI> const & alleles, GeneSource source, bool is_functional) {
    BEGIN();
    Gene * gene = get_gene_with_alleles(alleles);
    if(gene == nullptr) {
        gene = gene_manager.create();
        gene->alleles = alleles;
        gene->source = source;
        gene->is_functional = is_functional;
    
        if(OUTPUT_STRAINS) {
            write_gene(gene, gene_stmt, allele_stmt);
        }
    }
    RETURN(gene);
}

Strain * generate_random_strain(uint64_t n_new_genes, GeneSource new_gene_source) {
    BEGIN();
    
    std::array<Gene *, N_GENES_PER_STRAIN> genes;
    for(uint64_t i = 0; i < n_new_genes; i++) {
        genes[i] = mutate_gene(draw_random_gene(), new_gene_source);
    }
    for(uint64_t i = n_new_genes; i < N_GENES_PER_STRAIN; i++) {
        genes[i] = draw_random_gene();
    }
    
    Strain * strain = get_strain_with_genes(genes);
    RETURN(strain);
}

Strain * get_strain_with_genes(std::array<Gene *, N_GENES_PER_STRAIN> genes) {
    BEGIN();
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
        
        if(OUTPUT_STRAINS) {
            write_strain(strain, strain_stmt, NULL, NULL);
        }
    }
    else {
        strain = itr->second;
    }
    RETURN(strain);
}

Strain * recombine_strains(Strain * s1, Strain * s2) {
    BEGIN();
    
    // Choose random subset of all genes from strains s1 and s2
    // [0, N_GENES_PER_STRAIN) mapped to s1
    // [N_GENES_PER_STRAIN, 2 * N_GENES_PER_STRAIN) mapped to s2
    std::bitset<2 * N_GENES_PER_STRAIN> used;
    std::array<Gene *, N_GENES_PER_STRAIN> daughter_genes;
    for(uint64_t i = 0; i < N_GENES_PER_STRAIN; i++) {
        // Choose random unused index across strains
        uint64_t src_index;
        do {
            src_index = draw_uniform_index(2 * N_GENES_PER_STRAIN);
        } while(used[src_index]);
        used[src_index] = true;
        if(src_index < N_GENES_PER_STRAIN) {
            daughter_genes[i] = s1->genes_sorted[src_index]; 
        }
        else {
            daughter_genes[i] = s2->genes_sorted[src_index - N_GENES_PER_STRAIN];
        }
    }
    RETURN(get_strain_with_genes(daughter_genes));
}

double get_gene_similarity(Gene * gene1, Gene * gene2, uint64_t breakpoint) {
    BEGIN();
    
    double p_div = 0;
    double child_div = 0;
    double rho = 0.8; //recombination tolerance;
    double avg_mutation = 5; //average number of mutations per epitope
    for(uint64_t i = 0; i < N_LOCI; i++) {
        if(gene1->alleles[i] != gene2->alleles[i]) {
            p_div += 1;
            if(i < breakpoint) {
                child_div += 1;
            }
        }
    }
    double rho_power = child_div * avg_mutation * (p_div - child_div)* avg_mutation / (p_div * avg_mutation - 1);
    double surv_prob = pow(rho, rho_power);
    
    RETURN(surv_prob);
}

Strain * mutate_strain(Strain * strain) {
    BEGIN();
    uint64_t index = draw_uniform_index(N_GENES_PER_STRAIN);
    std::array<Gene *, N_GENES_PER_STRAIN> genes = strain->genes_sorted;
    genes[index] = mutate_gene(genes[index], SOURCE_MUTATION);
    RETURN(get_strain_with_genes(genes));
}

Gene * mutate_gene(Gene * gene, GeneSource source) {
    BEGIN();
    auto alleles = gene->alleles;
    uint64_t locus = draw_uniform_index(N_LOCI); 
    alleles[locus] = n_alleles[locus];
    n_alleles[locus]++;
    RETURN(get_or_create_gene(alleles, source, gene->is_functional));
}

void recombine_infection(Infection * infection) {
    BEGIN();
    
    if(N_GENES_PER_STRAIN == 1) {
        RETURN();
    }
    
    Gene * new_gene_1;
    Gene * new_gene_2;
    
    // Choose random expression order indices
    uint64_t exp_index_1;
    uint64_t exp_index_2;
    do {
        exp_index_1 = draw_uniform_index(N_GENES_PER_STRAIN);
        exp_index_2 = draw_uniform_index_except(N_GENES_PER_STRAIN, exp_index_1);
    } while(exp_index_1 == exp_index_2);
    
    Gene * src_gene_1 = infection->expression_order[exp_index_1];
    Gene * src_gene_2 = infection->expression_order[exp_index_2];
    if(src_gene_1 == src_gene_2) {
        RETURN();
    }
    
    bool is_conversion = draw_bernoulli(P_ECTOPIC_RECOMBINATION_IS_CONVERSION);
    uint64_t breakpoint = draw_uniform_index(N_LOCI);
    
    // If breakpoint == 0, very little to do
    if(breakpoint == 0) {
        if(is_conversion) {
            new_gene_1 = src_gene_1;
            new_gene_2 = src_gene_1;
        }
        else {
            new_gene_1 = src_gene_1;
            new_gene_2 = src_gene_2;
        }
    }
    // Otherwise, recombine the genes to produce one (under conversion) or two (if normal recombination) new genes
    else {
        double similarity = get_gene_similarity(src_gene_1, src_gene_2, breakpoint);
        
        auto rec_alleles_1 = recombine_alleles(src_gene_1->alleles, src_gene_2->alleles, breakpoint);
        Gene * rec_gene_1 = get_or_create_gene(rec_alleles_1, SOURCE_RECOMBINATION, draw_bernoulli(similarity));
        new_gene_1 = rec_gene_1->is_functional ? rec_gene_1 : src_gene_1;
        
        // Under conversion, the second gene remains unchanged
        if(is_conversion) {
            new_gene_2 = src_gene_2;
        }
        // Under normal combination, both genes have material swapped
        else {
            auto rec_alleles_2 = recombine_alleles(src_gene_2->alleles, src_gene_1->alleles, breakpoint);
            Gene * rec_gene_2 = get_or_create_gene(rec_alleles_2, SOURCE_RECOMBINATION, draw_bernoulli(similarity));
            new_gene_2 = rec_gene_2->is_functional ? rec_gene_2 : src_gene_2;
        }
    }
    
    // If nothing has changed, nothing to do
    if(new_gene_1 == src_gene_1 && new_gene_2 == src_gene_2) {
        RETURN();
    }
    
    // Update expression order and strain
    std::array<Gene *, N_GENES_PER_STRAIN> new_genes = infection->strain->genes_sorted;
    if(new_gene_1 != src_gene_1) {
        infection->expression_order[exp_index_1] = new_gene_1;
        replace_first(new_genes, src_gene_1, new_gene_1);
    }
    if(new_gene_2 != src_gene_2) {
        infection->expression_order[exp_index_2] = new_gene_2;
        replace_first(new_genes, src_gene_2, new_gene_2);
    }
    infection->strain = get_strain_with_genes(new_genes);
    
    RETURN();
}

std::array<uint64_t, N_LOCI> recombine_alleles(
    std::array<uint64_t, N_LOCI> const & a1, std::array<uint64_t, N_LOCI> const & a2, uint64_t breakpoint
) {
    BEGIN();
    std::array<uint64_t, N_LOCI> arc;
    for(uint64_t i = 0; i < N_LOCI; i++) {
        if(i < breakpoint) {
            arc[i] = a1[i];
        }
        else {
            arc[i] = a2[i];
        }
    }
    RETURN(arc);
}

bool contains_different_genes(Strain * strain) {
    BEGIN();
    for(uint64_t i = 1; i < N_GENES_PER_STRAIN; i++) {
        if(strain->genes_sorted[i] != strain->genes_sorted[0]) {
            RETURN(true);
        }
    }
    RETURN(false);
}

Gene * draw_random_gene() {
    BEGIN();
    Gene * gene = gene_manager.objects()[draw_uniform_index(gene_manager.size())];
    RETURN(gene);
}

Gene * get_current_gene(Infection * infection) {
    BEGIN();
    Gene * gene = infection->expression_order[infection->expression_index];
    RETURN(gene);
}

void gain_immunity(Host * host, Gene * gene) {
    BEGIN();
    ImmuneHistory * immune_history = host->immune_history;
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
    update_next_immunity_loss_time(host);
    RETURN();
}

void update_next_immunity_loss_time(Host * host) {
    BEGIN();
    if(SELECTION_MODE != SPECIFIC_IMMUNITY) {
        RETURN();
    }
    uint64_t immune_allele_count = get_immune_allele_count(host);
    host->next_immunity_loss_time = draw_exponential_after_now(IMMUNITY_LOSS_RATE * immune_allele_count);
    immunity_loss_queue.update(host);
    RETURN();
}

void lose_random_immunity(Host * host) {
    BEGIN();
    
    ImmuneHistory * immune_history = host->immune_history;
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
    update_next_immunity_loss_time(host);
    RETURN();
}

uint64_t get_immune_allele_count(Host * host) {
    BEGIN();
    uint64_t immune_allele_count = 0;
    ImmuneHistory * immune_history = host->immune_history;
    for(uint64_t i = 0; i < N_LOCI; i++) {
        LocusImmunity * immunity = immune_history->immunity_by_locus[i];
        immune_allele_count += immunity->immunity_level_by_allele.size();
    }
    RETURN(immune_allele_count);
}

double get_specific_immunity_level(Host * host, Gene * gene) {
    BEGIN();
    
    uint64_t immunity_count = 0;
    for(uint64_t i = 0; i < N_LOCI; i++) {
        auto immunity_level_by_allele = host->immune_history->immunity_by_locus[i]->immunity_level_by_allele; 
        auto itr = immunity_level_by_allele.find(gene->alleles[i]);
        if(itr != immunity_level_by_allele.end() && itr->second > 0) {
            immunity_count +=  1;
        }
    }
    
    RETURN(immunity_count / (double)N_LOCI);
}

uint64_t get_active_infection_count(Host * host) {
    BEGIN();
    uint64_t count = 0;
    for(Infection * infection : host->infections) {
        if(infection->expression_index >= 0) {
            count += 1;
        }
    }
    
    RETURN(count);
}

void infect_host(Host * host, Strain * strain) {
    BEGIN();
    
    Infection * infection = infection_manager.create();
    infection->host = host;
    infection->strain = strain;
    
    for(uint64_t i = 0; i < strain->genes_sorted.size(); i++) {
        infection->expression_order[i] = strain->genes_sorted[i];
    }
    std::shuffle(
        infection->expression_order.begin(), infection->expression_order.end(), rng
    );
    
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
    
    n_infections_cumulative++;
    
    RETURN();
}

void perform_infection_transition(Infection * infection) {
    BEGIN();
    
    Host * host = infection->host;
    if(infection->expression_index == -1) {
        infection->expression_index = 0;
    }
    else {
        assert(
            infection->expression_index >= 0 &&
            infection->expression_index < N_GENES_PER_STRAIN
        );
        Gene * gene = get_current_gene(infection);
        
        if(SELECTION_MODE == SPECIFIC_IMMUNITY) {
            gain_immunity(host, gene);
        }
        
        infection->expression_index++;
    }
    
    if(infection->expression_index == N_GENES_PER_STRAIN) {
        clear_infection(infection);
    }
    else {
        update_host_infection_times(infection->host);
    }
    
    RETURN();
}

void clear_infection(Infection * infection) {
    BEGIN();
    
    Host * host = infection->host;
    host->infections.erase(infection);
    
    destroy_infection(infection);
    host->completed_infection_count++;
    
    if(SELECTION_MODE == GENERAL_IMMUNITY) {
        update_host_infection_times(host);
    }
    
    RETURN();
}

void update_host_infection_times(Host * host) {
    BEGIN();
    for(Infection * infection : host->infections) {
        update_infection_times(infection);
    }
    RETURN();
}

void update_infection_times(Infection * infection) {
    BEGIN();
    update_transition_time(infection, false);
    RETURN();
}

void update_transition_time(Infection * infection, bool initial) {
    BEGIN();
    if(infection->expression_index == -1) {
        return;
    }
    
    assert(infection->expression_index < N_GENES_PER_STRAIN);
    
    Host * host = infection->host;
    double rate;
    if(SELECTION_MODE == SPECIFIC_IMMUNITY) {
        Gene * active_gene = get_current_gene(infection);
        double immunity_level =  get_specific_immunity_level(host, active_gene);
        assert(immunity_level >= 0.0 && immunity_level <= 1.0);
        rate = TRANSITION_RATE_IMMUNE * TRANSITION_RATE_NOT_IMMUNE / (
            TRANSITION_RATE_IMMUNE * (1.0 - immunity_level) +
            TRANSITION_RATE_NOT_IMMUNE * immunity_level
        );
    }
    else {
        rate = TRANSITION_RATE_NOT_IMMUNE;
    }
    infection->transition_time = draw_exponential_after_now(rate);
    
    if(initial) {
        transition_queue.add(infection);
    }
    else {
        transition_queue.update(infection);
    }
    RETURN();
}

void update_mutation_time(Infection * infection, bool initial) {
    BEGIN();
    infection->mutation_time = draw_exponential_after_now(
        MUTATION_RATE * N_GENES_PER_STRAIN * N_LOCI
    );
    if(initial) {
        mutation_queue.add(infection);
    }
    else {
        mutation_queue.update(infection);
    }
    RETURN();
}

void update_recombination_time(Infection * infection, bool initial) {
    BEGIN();
    infection->recombination_time = draw_exponential_after_now(
        ECTOPIC_RECOMBINATION_RATE * N_GENES_PER_STRAIN * (N_GENES_PER_STRAIN - 1) / 2.0
    );
    
    if(initial) {
        recombination_queue.add(infection);
    }
    else {
        recombination_queue.update(infection);
    }
    RETURN();
}

void update_clearance_time(Infection * infection, bool initial) {
    BEGIN();
    
    assert(SELECTION_MODE == GENERAL_IMMUNITY);
    
    // TODO: replace with Qixin's function based on vector of parameters
    double rate;
    if(infection->expression_index == -1) {
        rate = 0.0;
    }
    else {
        Host * host = infection->host;
        double a = GENERAL_IMMUNITY_PARAMS[0];
        double b = GENERAL_IMMUNITY_PARAMS[1];
        double c = GENERAL_IMMUNITY_PARAMS[2];
        double d = GENERAL_IMMUNITY_PARAMS[3];
         
        if(host->completed_infection_count < N_INFECTIONS_FOR_GENERAL_IMMUNITY) {
            rate = 1.0 / (a +
                b * exp(
                    -c * host->completed_infection_count
                ) /
                pow(
                    d * host->completed_infection_count + 1.0,
                    d
                )
            );
        }
        else {
            rate = CLEARANCE_RATE_IMMUNE;
        }
    }
    
    infection->clearance_time = draw_exponential_after_now(rate);
    if(initial) {
        clearance_queue.add(infection);
    }
    else {
        clearance_queue.update(infection);
    }
    RETURN();
}


#pragma mark \
*** Main simulation loop ***

enum class EventType {
    NONE,
    VERIFICATION,
    HOST_SAMPLING,
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

void run(bool override_seed, uint64_t random_seed) {
    BEGIN();
    
    initialize(override_seed, random_seed);
    while(do_next_event()) {  }
    clean_up();
    
    RETURN();
}

bool do_next_event() {
    BEGIN();
    double next_event_time = INF;
    EventType next_event_type = EventType::NONE;
    
    // Find the event type with the lowest-timed event (simple linear search)
    // This is the least elegant part of the reimplementation, but it should
    // be worth it for simplicity and reduced memory usage.
    // If this section is a speed bottleneck, re-evaluate the choice.
    // It shouldn't come even close.
    if(next_sampling_time < next_event_time) {
        next_event_time = next_sampling_time;
        next_event_type = EventType::HOST_SAMPLING;
    }
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
    
    PRINT_DEBUG(1, "next_event_time: %f", next_event_time);
    PRINT_DEBUG(1, "next_event_type: %d", next_event_type);
    
    // Execute the next event unless it's past T_END, in which case just advance time
    assert(next_event_time >= now);
    if(next_event_time > T_END) {
        now = T_END;
        return false;
    }
        
    now = next_event_time;
    switch(next_event_type) {
        case EventType::NONE:
            break;
        case EventType::VERIFICATION:
            do_verification_event();
            break;
        case EventType::HOST_SAMPLING:
            do_sampling_event();
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
    RETURN(true);
}

void do_verification_event() {
    BEGIN();
    verify_simulation_state();
    next_verification_time += VERIFICATION_PERIOD;
    RETURN();
}

void do_sampling_event() {
    BEGIN();
    write_summary();
    sample_hosts();
    next_sampling_time += HOST_SAMPLING_PERIOD;
    RETURN();
}

void do_checkpoint_event() {
    BEGIN();
    save_checkpoint();
    next_checkpoint_time += CHECKPOINT_SAVE_PERIOD;
    RETURN();
}

void do_biting_event() {
    BEGIN();
    assert(biting_queue.size() > 0);
    Population * pop = biting_queue.head();
    PRINT_DEBUG(1, "biting event pop: %llu", pop->id);
    Host * src_host = draw_random_source_host(pop); 
    Host * dst_host = draw_random_destination_host(pop);
    
    transmit(src_host, dst_host);
    
    // Update biting event time
    update_biting_time(pop, false);
    
    RETURN();
}

Host * draw_random_source_host(Population * pop) {
    BEGIN();
    Host * host = pop->hosts.object_at_index(draw_uniform_index(pop->hosts.size()));
    RETURN(host);
}

Host * draw_random_destination_host(Population * src_pop) {
    BEGIN();
    
    Host * host = nullptr;
    
    // Uniform random across populations using linear search
    // TODO: add distance-based method back in
    uint64_t n_hosts_total = accumulate(N_HOSTS);
    uint64_t index = draw_uniform_index(n_hosts_total);
    uint64_t offset = 0;
    for(Population * pop : population_manager.objects()) {
        if(index < offset + pop->hosts.size()) {
            host = pop->hosts.object_at_index(index - offset);
            break;
        }
        offset += pop->hosts.size();
    }
    assert(host != nullptr);
    
    RETURN(host);
}

void transmit(Host * src_host, Host * dst_host) {
    BEGIN();
    if(src_host->infections.size() > 0) {
        std::vector<Strain *> src_strains;
        for(Infection * infection : src_host->infections) {
            if(infection->expression_index >= 0) {
                if(draw_bernoulli(get_transmission_probability(infection))) {
                    src_strains.push_back(infection->strain);
                }
            }
        }
        
        // Form set of strains to transmit: some recombinants; some unmodified 
        std::vector<Strain *> strains_to_transmit(src_strains.size());
        
        // Produce a set of strains of the same size as src_strains
        for(uint64_t i = 0; i < src_strains.size(); i++) {
            Strain * strain1 = src_strains[draw_uniform_index(src_strains.size())];
            Strain * strain2 = src_strains[draw_uniform_index(src_strains.size())];
            
            // If they're the same, use them unchanged. Otherwise, recombine.
            if(strain1 == strain2) {
                strains_to_transmit[i] = strain1;
            }
            else {
                strains_to_transmit[i] = recombine_strains(strain1, strain2);
            }
        }
        
        for(Strain * strain : strains_to_transmit) {
            infect_host(dst_host, strain);
        }
    }    
    RETURN();
}

double get_transmission_probability(Infection * infection) {
    BEGIN();
    if(COINFECTION_REDUCES_TRANSMISSION) {
        RETURN(
            GENE_TRANSMISSIBILITY /
            get_active_infection_count(infection->host)
        );
    }
    RETURN(GENE_TRANSMISSIBILITY);
}

void do_immigration_event() {
    BEGIN();
    assert(IMMIGRATION_ON);
    assert(immigration_queue.size() > 0);
    Population * pop = immigration_queue.head();
    
    PRINT_DEBUG(1, "immigration event pop: %llu", pop->id);
    
    uint64_t n_new_genes;
    if(draw_bernoulli(P_IMMIGRATION_INCLUDES_NEW_GENES)) {
        n_new_genes = N_IMMIGRATION_NEW_GENES;
    }
    else {
        n_new_genes = 0;
    }
    
    Strain * strain = generate_random_strain(n_new_genes, SOURCE_IMMIGRATION);
    uint64_t index = draw_uniform_index(pop->hosts.size());
    Host * host = pop->hosts.object_at_index(index);
     
    infect_host(host, strain);
    
    // Update immigration event time
    update_immigration_time(pop, false);
    RETURN();
}

void do_immunity_loss_event() {
    BEGIN();
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
    RETURN();
}

void do_death_event() {
    BEGIN();
    assert(death_queue.size() > 0);
    
    Host * host = death_queue.head();
    Population * pop = host->population;
    destroy_host(host);
    Host * new_host = create_host(pop);
    PRINT_DEBUG(1, "Created new host id %llu in population %llu", new_host->id, pop->id);
    RETURN();
}

void do_transition_event() {
    BEGIN();
    assert(transition_queue.size() > 0);
    
    Infection * infection = transition_queue.head();
    perform_infection_transition(infection);
    RETURN();
}

void do_mutation_event() {
    BEGIN();
    assert(mutation_queue.size() > 0);
    
    Infection * infection = mutation_queue.head();
    infection->strain = mutate_strain(infection->strain);
    
    update_mutation_time(infection, false);
    update_transition_time(infection, false);
    update_next_immunity_loss_time(infection->host);
    
    RETURN();
}

void do_recombination_event() {
    BEGIN();
    assert(recombination_queue.size() > 0);
    
    Infection * infection = recombination_queue.head();
    recombine_infection(infection);
    
    update_recombination_time(infection, false);
    update_transition_time(infection, false);
    update_next_immunity_loss_time(infection->host);
    
    RETURN();
}

void do_clearance_event() {
    BEGIN();
    assert(SELECTION_MODE == GENERAL_IMMUNITY);
    Infection * infection = clearance_queue.head();
    clear_infection(infection);
    RETURN();
}

#pragma mark \
*** Verification function implementations ***

void verify_simulation_state() {
    BEGIN();
    
    assert(biting_queue.verify_heap());
    assert(immigration_queue.verify_heap());
    assert(immunity_loss_queue.verify_heap());
    assert(death_queue.verify_heap());
    assert(transition_queue.verify_heap());
    assert(mutation_queue.verify_heap());
    assert(recombination_queue.verify_heap());
    assert(clearance_queue.verify_heap());
    
    RETURN();
}

#pragma mark \
*** Sampling implementations ***

void sample_hosts() {
    BEGIN();
    
    std::vector<Host *> infected_hosts;
    for(Host * host : host_manager.objects()) {
        if(host->infections.size() > 0) {
            infected_hosts.push_back(host);
        }
    }
    
    std::unordered_set<Host *> sampled_hosts;
    uint64_t sample_size = infected_hosts.size() < HOST_SAMPLE_SIZE ? infected_hosts.size() : HOST_SAMPLE_SIZE;
    
    for(uint64_t index : draw_uniform_indices_without_replacement(infected_hosts.size(), sample_size)) {
        Host * host = infected_hosts[index];
        write_sampled_host(host);
        
        for(Infection * infection : host->infections) {
            write_sampled_infection(host, infection);
        }
    }
    
    // Wrap up the transaction, including all the hosts, strains, genes created since the last sampling event
    sqlite3_exec(sample_db, "COMMIT", NULL, NULL, NULL);
    
    // Begin a new transaction to include all the hosts, strains, genes created before the next sampling event
    sqlite3_exec(sample_db, "BEGIN TRANSACTION", NULL, NULL, NULL);
    
    RETURN();
}

void write_summary() {
    BEGIN();
    
    uint64_t n_infected = 0;
    uint64_t n_infections = 0;
    
    std::unordered_set<Strain *> distinct_strains;
    std::unordered_set<Gene *> distinct_genes;
    std::array<std::unordered_set<uint64_t>, N_LOCI> distinct_alleles; 
    
    for(Population * pop : population_manager.objects()) {
        for(Host * host : pop->hosts.as_vector()) {
            if(host->infections.size() > 0) {
                n_infected++;
            }
            n_infections += host->infections.size();
            
            for(Infection * infection : host->infections) {
                Strain * strain = infection->strain;
                distinct_strains.insert(strain);
                for(Gene * gene : strain->genes_sorted) {
                    distinct_genes.insert(gene);
                    for(uint64_t i = 0; i < N_LOCI; i++) {
                        distinct_alleles[i].insert(gene->alleles[i]);
                    }
                }
            }
        }
    }
    
    printf("Summary at t = %f:\n", now);
    printf("               n_infections: %llu\n", n_infections);
    printf("                 n_infected: %llu\n", n_infected);
    printf("    n_infections_cumulative: %llu\n", n_infections_cumulative);
    printf("      n_circulating_strains: %lu\n", distinct_strains.size());
    printf("        n_circulating_genes: %lu\n", distinct_genes.size());
    
    sqlite3_bind_double(summary_stmt, 1, now); // time
    sqlite3_bind_int64(summary_stmt, 2, n_infections); // n_infections
    sqlite3_bind_int64(summary_stmt, 3, n_infected); // n_infected
    sqlite3_bind_int64(summary_stmt, 4, n_infections_cumulative); // n_infections_cumulative
    sqlite3_bind_int64(summary_stmt, 5, distinct_strains.size()); // n_circulating_strains
    sqlite3_bind_int64(summary_stmt, 6, distinct_genes.size()); // n_circulating_genes
    sqlite3_step(summary_stmt);
    sqlite3_reset(summary_stmt);
    
    for(uint64_t i = 0; i < N_LOCI; i++) {
        sqlite3_bind_double(summary_alleles_stmt, 1, now); // time
        sqlite3_bind_int64(summary_alleles_stmt, 2, i); // locus
        sqlite3_bind_int64(summary_alleles_stmt, 3, distinct_alleles[i].size()); // n_circulating_alleles
        sqlite3_step(summary_alleles_stmt);
        sqlite3_reset(summary_alleles_stmt);
    }
    
    RETURN();
}

void write_host(Host * host) {
    sqlite3_bind_int64(host_stmt, 1, host->id);
    sqlite3_bind_int64(host_stmt, 2, host->population->id);
    sqlite3_bind_double(host_stmt, 3, host->birth_time);
    sqlite3_bind_double(host_stmt, 4, host->death_time);
    sqlite3_step(host_stmt);
    sqlite3_reset(host_stmt);
}

void write_sampled_host(Host * host) {
    sqlite3_bind_double(sampled_host_stmt, 1, now);
    sqlite3_bind_int64(sampled_host_stmt, 2, host->id);
    sqlite3_bind_int64(sampled_host_stmt, 3, host->population->id);
    sqlite3_bind_double(sampled_host_stmt, 4, host->birth_time);
    sqlite3_bind_double(sampled_host_stmt, 5, host->death_time);
    sqlite3_step(sampled_host_stmt);
    sqlite3_reset(sampled_host_stmt);
}

void write_sampled_infection(Host * host, Infection * infection) {
    Strain * strain = infection->strain;
    
    sqlite3_bind_double(sampled_inf_stmt, 1, now);
    sqlite3_bind_int64(sampled_inf_stmt, 2, host->id);
    sqlite3_bind_int64(sampled_inf_stmt, 3, infection->id);
    sqlite3_bind_int64(sampled_inf_stmt, 4, strain->id);
    sqlite3_step(sampled_inf_stmt);
    sqlite3_reset(sampled_inf_stmt);
    
    write_strain(strain, sampled_strain_stmt, sampled_gene_stmt, sampled_allele_stmt);
}

void write_strain(Strain * strain, sqlite3_stmt * s_stmt, sqlite3_stmt * g_stmt, sqlite3_stmt * a_stmt) {
    for(uint64_t i = 0; i < N_GENES_PER_STRAIN; i++) {
        Gene * gene = strain->genes_sorted[i];
        sqlite3_bind_int64(s_stmt, 1, strain->id);
        sqlite3_bind_int64(s_stmt, 2, i);
        sqlite3_bind_int64(s_stmt, 3, gene->id);
        sqlite3_step(s_stmt);
        sqlite3_reset(s_stmt);
        
        if(g_stmt != NULL) {
            write_gene(gene, g_stmt, a_stmt);
        }
    }
}

void write_gene(Gene * gene, sqlite3_stmt * g_stmt, sqlite3_stmt * a_stmt) {
    sqlite3_bind_int64(g_stmt, 1, gene->id);
    sqlite3_bind_int64(g_stmt, 2, gene->source);
    sqlite3_bind_int64(g_stmt, 3, gene->is_functional);
    sqlite3_step(g_stmt);
    sqlite3_reset(g_stmt);
    
    for(uint64_t j = 0; j < N_LOCI; j++) {
        sqlite3_bind_int64(a_stmt, 1, gene->id);
        sqlite3_bind_int64(a_stmt, 2, j);
        sqlite3_bind_int64(a_stmt, 3, gene->alleles[j]);
        sqlite3_step(a_stmt);
        sqlite3_reset(a_stmt);
    }

}

#pragma mark \
*** Checkpoint function implementations ***

void save_checkpoint() {
    BEGIN();
    std::string old_checkpoint_filename = CHECKPOINT_SAVE_FILENAME + "-old";
    if(file_exists(CHECKPOINT_SAVE_FILENAME)) {
        assert(!rename(CHECKPOINT_SAVE_FILENAME.c_str(), old_checkpoint_filename.c_str()));
    }
    
    sqlite3 * db;
    assert(!file_exists(CHECKPOINT_SAVE_FILENAME));
    sqlite3_open(CHECKPOINT_SAVE_FILENAME.c_str(), &db);
    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, NULL);
    
    save_global_state_to_checkpoint(db);
    
    strain_manager.save_to_checkpoint(db);
    gene_manager.save_to_checkpoint(db);
    population_manager.save_to_checkpoint(db);
    host_manager.save_to_checkpoint(db);
    infection_manager.save_to_checkpoint(db);
    immune_history_manager.save_to_checkpoint(db);
    locus_immunity_manager.save_to_checkpoint(db);
    
    sqlite3_exec(db, "COMMIT", NULL, NULL, NULL);
    sqlite3_close(db);
    
    if(file_exists(old_checkpoint_filename)) {
        assert(!unlink(old_checkpoint_filename.c_str()));
    }
    RETURN();
}

void load_checkpoint(bool should_load_rng_state) {
    BEGIN();
    
    assert(file_exists(CHECKPOINT_LOAD_FILENAME));
    
    sqlite3 * db;
    sqlite3_open(CHECKPOINT_LOAD_FILENAME.c_str(), &db);
    
    load_global_state_from_checkpoint(db, should_load_rng_state);
    
    // Load all objects, minus references to other objects
    strain_manager.load_from_checkpoint(db);
    gene_manager.load_from_checkpoint(db);
    population_manager.load_from_checkpoint(db);
    host_manager.load_from_checkpoint(db);
    infection_manager.load_from_checkpoint(db);
    immune_history_manager.load_from_checkpoint(db);
    locus_immunity_manager.load_from_checkpoint(db);
    
    // Resolve references to other objects
    strain_manager.resolve_references(db, gene_manager);
    gene_manager.resolve_references(db); // Does nothing unless referenes are added to Gene
    population_manager.resolve_references(db, host_manager);
    host_manager.resolve_references(db, population_manager, immune_history_manager, infection_manager);
    infection_manager.resolve_references(db, strain_manager, host_manager, gene_manager);
    immune_history_manager.resolve_references(db, locus_immunity_manager);
    locus_immunity_manager.resolve_references(db); // Does nothing unless referenes are added to LocusImmunity
    
    sqlite3_close(db);
    
    initialize_event_queues_from_state();
    
    RETURN();
}

void save_global_state_to_checkpoint(sqlite3 * db) {
    sqlite3_exec(db,
        "CREATE TABLE global_state ("
            "rng TEXT, "
            "now REAL, "
            "next_verification_time REAL, "
            "next_checkpoint_time REAL, "
            "next_info_time REAL, "
            "n_infections_cumulative INTEGER"
        ");",
        NULL, NULL, NULL
    );
    
    sqlite3_stmt * stmt;
    sqlite3_prepare_v2(db, "INSERT INTO global_state VALUES (?,?,?,?,?,?);", -1, &stmt, NULL);
    std::string rng_str = get_rng_as_string();
    sqlite3_bind_text(stmt, 1, rng_str.c_str(), (int)(rng_str.size() + 1), SQLITE_STATIC);
    sqlite3_bind_double(stmt, 2, now);
    sqlite3_bind_double(stmt, 3, next_verification_time);
    sqlite3_bind_double(stmt, 4, next_checkpoint_time);
    sqlite3_bind_double(stmt, 5, next_info_time);
    sqlite3_bind_int64(stmt, 6, n_infections_cumulative);
    sqlite3_step(stmt);
    sqlite3_finalize(stmt);
}

void load_global_state_from_checkpoint(sqlite3 * db, bool should_load_rng_state) {
    sqlite3_stmt * stmt;
    sqlite3_prepare_v2(db, "SELECT * FROM global_state LIMIT 1;", -1, &stmt, NULL);
    sqlite3_step(stmt);
    
    if(should_load_rng_state) {
        // Load rng from string representation
        std::string rng_str = (const char *)sqlite3_column_text(stmt, 0);
        PRINT_DEBUG(1, "rng read: %s", rng_str.c_str());
        set_rng_from_string(rng_str);
        assert(rng_str == get_rng_as_string());
    }
    
    now = sqlite3_column_double(stmt, 1);
    next_verification_time = sqlite3_column_double(stmt, 2);
    next_checkpoint_time = sqlite3_column_double(stmt, 3);
    next_info_time = sqlite3_column_double(stmt, 4);
    n_infections_cumulative = sqlite3_column_int(stmt, 5);
    
    sqlite3_finalize(stmt);
}

void initialize_event_queues_from_state() {
    for(Population * pop : population_manager.objects()) {
        biting_queue.add(pop);
        immigration_queue.add(pop);
    }
    
    for(Host * host : host_manager.objects()) {
        immunity_loss_queue.add(host);
        death_queue.add(host);
    }
    
    for(Infection * infection : infection_manager.objects()) {
        transition_queue.add(infection);
        mutation_queue.add(infection);
        recombination_queue.add(infection);
        clearance_queue.add(infection);
    }
}

std::string get_rng_as_string() {
    std::stringstream ss;
    ss << rng;
    return ss.str();
}

void set_rng_from_string(std::string const & rng_str) {
    std::stringstream ss;
    ss << rng_str;
    ss >> rng;
}


#pragma mark \
*** Random draw helper function implementations *** 

double draw_exponential_after_now(double lambda) {
    BEGIN();
    double time;
    if(lambda == 0.0) {
        time = INF;
    }
    else {
        time = now + draw_exponential(lambda);
    }
    RETURN(time);
}

double draw_exponential(double lambda) {
    BEGIN();
    RETURN(std::exponential_distribution<>(lambda)(rng));
}

uint64_t draw_uniform_index(uint64_t size) {
    BEGIN();
    RETURN(std::uniform_int_distribution<uint64_t>(0, size - 1)(rng));
}

uint64_t draw_uniform_index_except(uint64_t size, uint64_t except_index) {
    BEGIN();
    uint64_t index = draw_uniform_index(size - 1);
    if(index >= except_index) {
        index++;
    }
    return(index);
}

double draw_uniform_real(double min, double max) {
    BEGIN();
    RETURN(std::uniform_real_distribution<double>(min, max)(rng));
}

std::vector<uint64_t> draw_uniform_indices_without_replacement(uint64_t n, uint64_t k) {
    BEGIN();
    
    assert(k <= n);
    
    std::uniform_int_distribution<uint64_t> index_dist(0, n - 1);
    
    // If we're drawing less than half, draw indices to *include*
    std::vector<uint64_t> indices;
    if(k < n / 2) {
        std::unordered_set<uint64_t> include_indices;
        for(uint64_t i = 0; i < k; i++) {
            while(true) {
                uint64_t index = index_dist(rng);
                if(include_indices.find(index) == include_indices.end()) {
                    include_indices.insert(index);
                    break;
                }
            }
        }
        assert(include_indices.size() == k);
        indices.insert(indices.end(), include_indices.begin(), include_indices.end());
        std::sort(indices.begin(), indices.end());
    }
    // Otherwise draw indices to *exclude*
    else {
        std::unordered_set<uint64_t> exclude_indices;
        for(uint64_t i = 0; i < n - k; i++) {
            while(true) {
                uint64_t index = index_dist(rng);
                if(exclude_indices.find(index) == exclude_indices.end()) {
                    exclude_indices.insert(index);
                    break;
                }
            }
        }
        assert(exclude_indices.size() == n - k);
        for(uint64_t i = 0; i < n; i++) {
            if(exclude_indices.find(i) == exclude_indices.end()) {
                indices.push_back(i);
            }
        }
    }
    assert(indices.size() == k);
    RETURN(indices);
}

bool draw_bernoulli(double p) {
    BEGIN();
    RETURN(std::bernoulli_distribution(p)(rng));
}

void update_biting_time(Population * pop, bool initial) {
    BEGIN();
    double biting_rate = BITING_RATE_MEAN[pop->ind] * (
        1.0 + BITING_RATE_RELATIVE_AMPLITUDE[pop->ind] * cos(
            2 * M_PI * ((now / T_YEAR) - BITING_RATE_PEAK_PHASE[pop->ind])
        )
    );
    pop->next_biting_time = draw_exponential_after_now(biting_rate);
    if(initial) {
        biting_queue.add(pop);
    }
    else {
        biting_queue.update(pop);
    }
    RETURN();
}

void update_immigration_time(Population * pop, bool initial) {
    BEGIN();
    pop->next_immigration_time = draw_exponential_after_now(IMMIGRATION_RATE[pop->ind]);
    if(initial) {
        immigration_queue.add(pop);
    }
    else {
        immigration_queue.update(pop);
    }
    RETURN();
}

} // namespace varmodel
