RANDOM_SEED = 1

N_GENES_INITIAL = 100
N_GENES_PER_STRAIN = 60
N_LOCI = 2
N_ALLELES_INITIAL = N_LOCI * [1200]

SELECTION_MODE = 'SPECIFIC_IMMUNITY'
N_INFECTIONS_FOR_GENERAL_IMMUNITY = 0
GENERAL_IMMUNITY_PARAMS = [-1, -1, -1, -1]

TRANSITION_RATE_NOT_IMMUNE = 1.0 / 12.0
CLEARANCE_RATE_IMMUNE = 0.0001

BITING_RATE_MEAN = [0.2]

PRINT_FUNCTION_TRACE = False
PRINT_DEBUG_LEVEL = 0

T_YEAR = 360.0
T_BURNIN = 100.0 * T_YEAR
T_END = 200.0 * T_YEAR

SAMPLE_DB_FILENAME = '"sample_db.sqlite"'

PRINT_INFO_PERIOD = 360.0

VERIFICATION_ON = True
VERIFICATION_PERIOD = 50.0 * T_YEAR

SAVE_TO_CHECKPOINT = True
CHECKPOINT_SAVE_FILENAME = '"checkpoint_db.sqlite"'
CHECKPOINT_SAVE_PERIOD = T_YEAR * 20

LOAD_FROM_CHECKPOINT = False
CHECKPOINT_LOAD_FILENAME = '""'

OUTPUT_HOSTS = True
OUTPUT_GENES = True
OUTPUT_STRAINS = True

HOST_SAMPLING_PERIOD = 30.0

GENE_TRANSMISSIBILITY = 0.5
COINFECTION_REDUCES_TRANSMISSION = True

ECTOPIC_RECOMBINATION_RATE = 1.8e-07
P_ECTOPIC_RECOMBINATION_IS_CONVERSION = 0.0

IMMUNITY_LOSS_RATE = 0.001

MUTATION_RATE = 1.42e-08

T_LIVER_STAGE = 14.0

TRANSITION_RATE_IMMUNE = 1000.0

CLEARANCE_RATE_IMMUNE = 0.041628833086078301

MEAN_HOST_LIFETIME = 30.0 * T_YEAR
MAX_HOST_LIFETIME = 80.0 * T_YEAR

N_POPULATIONS = 1

N_HOSTS = [10000]
N_INITIAL_INFECTIONS = [20]

BITING_RATE_RELATIVE_AMPLITUDE = [0.0]
BITING_RATE_PEAK_PHASE = [0.0]

IMMIGRATION_ON = True
IMMIGRATION_RATE = [1.0]
