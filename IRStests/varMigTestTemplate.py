RANDOM_SEED={RANDOM_SEED}
N_GENES_INITIAL={N_GENES_INITIAL}
N_GENES_PER_STRAIN=60
N_LOCI=2
N_ALLELES_INITIAL=N_LOCI*[{N_ALLELES_PER_LOCUS}]
SELECTION_MODE='{SELECTION_MODE}'
N_INFECTIONS_FOR_GENERAL_IMMUNITY=245
GENERAL_IMMUNITY_PARAMS=[51.4043638136979,284.518878105621,-0.0016944082357847,0.277409542527348]
CLEARANCE_RATE_IMMUNE=0.00511272898496179
TRANSITION_RATE_NOT_IMMUNE=1.0 / 6.0 # switching rate. every 6 days
TRANSITION_RATE_IMMUNE=1000
PRINT_FUNCTION_TRACE=False
PRINT_DEBUG_LEVEL=0
T_YEAR=360.0
T_BURNIN={T_BURNIN}
T_END={T_END}
EXPECTED_EQUILIBRIUM={EXPECTED_EQUILIBRIUM}
SAMPLE_DB_FILENAME='"{SAMPLE_DB_FILENAME}"'
PRINT_INFO_PERIOD=120.0
VERIFICATION_ON=False
VERIFICATION_PERIOD=39960
SAVE_TO_CHECKPOINT={SAVE_TO_CHECKPOINT}
CHECKPOINT_SAVE_FILENAME='"{CHECKPOINT_SAVE_FILENAME}"'
CHECKPOINT_SAVE_PERIOD=T_END
LOAD_FROM_CHECKPOINT={IRS}
CHECKPOINT_LOAD_FILENAME='"{CHECKPOINT_LOAD_FILENAME}"'
OUTPUT_HOSTS=True
OUTPUT_GENES=False # Set to false to save space
OUTPUT_STRAINS=False # Set to false to save space
HOST_SAMPLING_ON=True
HOST_SAMPLING_PERIOD=30.0
HOST_SAMPLE_SIZE=100
GENE_TRANSMISSIBILITY=0.5
COINFECTION_REDUCES_TRANSMISSION=True
ECTOPIC_RECOMBINATION_RATE=1.8e-07
P_ECTOPIC_RECOMBINATION_IS_CONVERSION=0
IMMUNITY_LOSS_RATE=0.00001
MUTATION_RATE=1.42e-08
T_LIVER_STAGE=14.0
MEAN_HOST_LIFETIME=30.0 * T_YEAR
MAX_HOST_LIFETIME=80.0 * T_YEAR
N_POPULATIONS=1
DISTANCE_MAT=[[1]]
DIST_POWER=1
N_HOSTS=[10000]
N_INITIAL_INFECTIONS=[20]
#BITING_RATE_MEAN=[0.0005]
BITING_RATE_MEAN=[{BITING_RATE_MEAN}]
BITING_RATE_RELATIVE_AMPLITUDE=[0]
BITING_RATE_PEAK_PHASE=[0]
DAILY_BITING_RATE_DISTRIBUTION=[{DAILY_BITING_RATE_DISTRIBUTION}]
IRS_ON={IRS}
IRS_START_TIMES=[{IRS_START}] # A vector with start times for the IRS interventions. Each element is the starting time of an IR,29160]
BITING_RATE_FACTORS=[[{BITING_RATE_FACTORS}]] # An array in which each element is a vector of values of mosquito numbers as the ones in DAILY_BITING_RATE_DISTRIBUTION. These are obtained from the model externally run in Mathematic,
IRS_IMMIGRATION_RATE_FACTORS=[1] # A vector with length as number of IRS events. Each element is a factor to multiply the IMMIGRATION_RATE. For example, 0.3 will reduce the usual immigration rate to 30% of its original size,0]
MDA_ON=False
MDA_START_TIMES=[]
HOST_FAIL_RATE=[] # % of hosts that did not take the drug
DRUG_EFF_DURATION=[] # How long the drug remains effective in the body
MDA_IMMIGRATION_RATE_FACTORS=[]
IMMIGRATION_ON=True
IMMIGRATION_RATE=[0.0026]
REGION_TO_LOCAL_POP_SIZE_RATIO = 1