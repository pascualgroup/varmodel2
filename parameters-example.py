PRINT_FUNCTION_TRACE = True
PRINT_DEBUG_LEVEL = 0

SELECTION_MODE = 'NEUTRALITY'

RANDOM_SEED = 100
T_YEAR = 360.0
T_BURNIN = 100.0 * T_YEAR
T_END = 200.0 * T_YEAR

SAMPLE_DB_FILENAME = '"sample_db.sqlite"'

VERIFICATION_PERIOD = 1.0

SAVE_TO_CHECKPOINT = False
CHECKPOINT_SAVE_FILENAME = '"checkpoint_db.sqlite"'
CHECKPOINT_SAVE_PERIOD = T_YEAR * 20
LOAD_FROM_CHECKPOINT = False
CHECKPOINT_LOAD_FILENAME = '""'

OUTPUT_HOSTS = True
OUTPUT_GENES = True
OUTPUT_STRAINS = True

HOST_SAMPLING_PERIOD = 30.0

SAMPLE_TRANSMISSION_EVENTS = False
TRANSMISSION_EVENT_SAMPLING_SKIP = 1000

N_GENES_IN_POOL = 50
N_GENES_PER_STRAIN = 5
N_LOCI = 10
N_ALLELES = N_LOCI * [5]

GENE_TRANSMISSIBILITY = 1.0
IMMUNITY_LOSS_RATE = 1.0

P_MUTATION = 0.01
P_STRAIN_RECOMBINATION = 0.05

T_LIVER_STAGE = 7.0

ACTIVATION_RATE = 100.0
DEACTIVATION_RATE_NOT_IMMUNE = 25.0/300.0
DEACTIVATION_RATE_IMMUNE = 1000.0
CLEARANCE_RATE = 0.0001

USE_ALLELE_IMMUNITY = False

USE_HOST_LIFETIME_PDF = False
HOST_LIFETIME_PDF = [0.3] * 20 + [0.4] * 10
MEAN_HOST_LIFETIME = 30.0 * T_YEAR
MAX_HOST_LIFETIME = 80.0 * T_YEAR

N_POPULATIONS = 1

N_HOSTS = [10000]
N_INITIAL_INFECTIONS = [100]

BITING_RATE_MEAN = [1.4]
BITING_RATE_RELATIVE_AMPLITUDE = [0.05]
BITING_RATE_PEAK_PHASE = [0.0]

IMMIGRATION_RATE = [0.1]
