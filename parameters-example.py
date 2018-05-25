RANDOM_SEED=7335537
N_GENES_INITIAL=12000
N_GENES_PER_STRAIN=60
N_LOCI=2
N_ALLELES_INITIAL=N_LOCI*[1200]
SELECTION_MODE='SPECIFIC_IMMUNITY'
N_INFECTIONS_FOR_GENERAL_IMMUNITY=245
GENERAL_IMMUNITY_PARAMS=[51.4043638136979,284.518878105621,-0.0016944082357847,0.277409542527348]
CLEARANCE_RATE_IMMUNE=0.00511272898496179
TRANSITION_RATE_NOT_IMMUNE=1.0 / 6.0 # switching rate. every 6 days
TRANSITION_RATE_IMMUNE=1000
PRINT_FUNCTION_TRACE=False
PRINT_DEBUG_LEVEL=0
T_YEAR=360.0
T_BURNIN=15
T_END=2000
SAMPLE_DB_FILENAME='"test.sqlite"'
PRINT_INFO_PERIOD=120.0
VERIFICATION_ON=False
VERIFICATION_PERIOD=10800
SAVE_TO_CHECKPOINT=False
CHECKPOINT_SAVE_FILENAME='"tesck.sqlite"'
CHECKPOINT_SAVE_PERIOD=1000
LOAD_FROM_CHECKPOINT=True
CHECKPOINT_LOAD_FILENAME='"tesck.sqlite"'
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
IMMUNITY_LOSS_RATE=0.001
MUTATION_RATE=1.42e-08
T_LIVER_STAGE=14.0
MEAN_HOST_LIFETIME=30.0 * T_YEAR
MAX_HOST_LIFETIME=80.0 * T_YEAR
N_POPULATIONS=1
N_HOSTS=[10000]
N_INITIAL_INFECTIONS=[20]
BITING_RATE_MEAN=[0.00005]
BITING_RATE_RELATIVE_AMPLITUDE=[0.0]
BITING_RATE_PEAK_PHASE=[0.0]
DAILY_BITING_RATE_DISTRIBUTION=[184.982060636238,200.255404185268,190.091657432075,176.936673923227,164.343031052803,152.793760078616,142.280309311918,132.721990757893,124.033625709693,116.13590946434,108.956491569867,102.429591067436,96.4954777083097,91.0998775063187,86.1935081934926,81.7316148582299,77.6735473752113,73.9823941296359,70.624631409777,67.569823745941,64.7903263402798,62.261040675501,59.9591757095815,57.8640360180417,55.9568281889751,54.2204871570074,52.6395139643939,51.1998352707491,49.8886665452749,48.6943983409682,47.6064821131221,46.6153343553555,45.7122530802748,44.889322061416,44.13935303385,43.4558081247075,42.8327483785961,42.2647680297853,41.7469538312647,41.2748368576276,40.8443501562017,40.4517970480741,40.0938070759306,39.7673175268291,39.4695383049637,39.1979296104896,38.9501783120656,38.7241775187336,38.5180090184007,38.3299237028482,38.1583293639103,38.0017737769016,37.8589341338113,37.7286047949206,37.609686843514,37.50117832764,37.4021653181591,37.3118147288408,37.2293672736381,37.1541298330001,37.0854713133807,37.0228146684226,36.9656352249567,36.9134526965649,36.8658304874747,36.8223689332749,36.7827047634528,36.7465053716365,36.7134680370484,36.683316586589,36.6557984022823,36.6306826034548,36.6077597837153,36.586838896635,36.5677447120854,36.5503175349244,36.5344117986761,36.5198946799105,36.5066448234627,36.4945517554252,36.4835144475775,36.4734406657602,36.4642461883492,36.4558541345434,36.4481944766219,36.4412033547379,36.4348224858472,36.4289985845382,36.423682949872,36.4188312037933,36.4144029157693,36.4103610895911,36.4066719764741,36.4033048155977,36.4002316266374,36.3974265546817,36.3948664137349,36.3925295614191,36.390396747356,36.3884499478837,36.3866731020119,36.3850512618356,36.383570977466,36.3822198706082,36.3809866531596,36.3798611099218,36.3788337913255,36.3778961228354,36.3770403186877,36.376259127581,36.3755461792137,36.3748953467457,36.3743013727248,36.3737591652055,36.3732643065013,36.3728126044509,36.3724003301053,36.3720240316874,36.3716805650306,36.3713671195684,36.3710809645286,36.3708198308674,36.3705814728698,36.3703638797406,36.370165316919,36.3699840324805,36.369818578311,36.3696675718257,36.3695297102287,36.3694039040808,36.3692890594064,36.3691842391427,36.3690885778915,36.3690012449681,36.368921567071,36.3688488589252,36.3687824611981,36.3687218882497,36.3686666214112,36.368616137142,36.368570066055,36.3685280375853,36.3684896516034,36.3684545980557,36.3684226195349,36.3683934248545,36.3683667563531,36.3683424240487,36.3683202157831,36.3682999305956,36.368281422106,36.3682645324877,36.368249105087,36.3682350279662,36.3682221849002,36.3682104538153,36.3681997446467,36.3681899755132,36.3681810545177,36.3681729107047,36.3851737913722,36.6273491430708,37.4853628842298,39.3391521070744,42.5139152093929,47.287707648272,53.9085580483502,62.6109722095744,73.6276001726413,87.1941215922183,103.547494101901,122.919379658258,145.527298129489,171.565847125294,201.199569715256,234.558182047779,271.734187710335,312.782500620207,357.721545950115,406.535310033193,459.175902687809,515.566288170545,575.603037004223,639.158875127664,706.085079258744,776.213619996515,849.359092975813,925.320508432896,1003.88286873039,1084.81866407249,1167.88921717579,1252.84603432857,1339.43197501728,1427.38246518551,1516.42665703923,1606.288540236,1696.68804794844,1787.34213872482,1877.96586367714,1968.27341928461,2057.97918561558,2146.79874988722,2234.44991407981,2320.65368564903,2405.13524893337,2487.62491634739,2567.859055995,2645.58099505399,2720.54189505184,2792.50159831076,2861.22944219628,2926.50503922685,2988.11902160115,3045.87374666899,3099.58396234427,3149.0774305219,3194.19550599262,3234.79366968277,3270.74201453417,3301.92568255082,3328.24525171735,3349.61707166055,3365.97354703317,3377.26336782105,3383.45168589349,3384.52023734929,3380.46741043488,3371.30825831918,3357.07445753911,3337.81421221584,3313.59210333082,3284.48888537285,3250.60122960876,3212.04141547099,3168.93697156399,3121.4302662414,3069.67805102919,3013.85095644985,2954.13294358876,2890.72071254956,2823.82306974722,2753.66025685211,2680.46324274245,2604.472983022,2525.93965668524,2445.12188176186,2362.28573918503,2277.70407924451,2191.65575876302,2104.42450326994,2016.29812452819,1927.56762428701,1838.52625981138,1749.46862102118,1660.68973822427,1572.48417420318,1485.14512218572,1398.96353218762,1314.22725642948,1231.22022187782,1150.22162457876,1071.50518305081,995.338420555275,921.981996241624,851.689079603363,784.704748808794,721.265363911465,661.597778657233,605.918111912624,554.429624489276,507.29104204215,464.393290873619,425.407641599444,389.985626692259,357.802524436665,328.562007641641,301.994626224167,277.855605119689,255.922545136195,235.993412158679,217.884657199491,201.429550269167,186.476615850599,172.888251531252,160.539461254969,149.316695756803,139.116797656341,129.846076373408,121.419412190101,113.759501549716,106.796114867289,100.465476414138,94.709648242056,89.4760212533734,84.7168040863485,80.3886048480315,76.4520068339136,72.8712249533315,69.6137536857966,66.6500844545312,63.9534135854316,61.4994098098765,59.2659720392144,57.2330380312933,55.3823836622487,53.69746513125,52.1632546823583,50.7661037923442,49.4936221913513,48.3345467621114,47.2786509575201,46.316642834647,45.4400794217592,44.6412870359932,43.9132945674952,43.2497596404901,42.6449197270064,42.0935343856084,41.5908356082067,41.1324876618248,40.714549150696,40.3334264715257,39.9858542704671,39.6688594745203,39.379735021602,39.1160158594602,38.8754565646507,38.6560144266936,38.455825228768,38.2731932483092,38.1065711198888,37.9545499808064,37.8158473626304,37.6892902381191,37.5738142292774,37.4684449625939,37.3722953501391,37.2845577495937,37.2044940467375,37.1314316812032,37.0647574694356,37.0039118006962,36.948384199192,36.8977096264088,36.8514623263637,36.8092571039936,36.770737875963,36.7355843878313,36.7035008724317,36.6742197285944,36.6474959600555,36.6231057548678,36.6008445488391,36.5805270876016,36.5619839390151,36.5450597868185,36.5296131069445,36.5155148897651,36.5026474377478,36.4909032774043]
IRS_ON=False
IRS_START_TIMES=[]
BITING_RATE_FACTORS=[]
IRS_IMMIGRATION_RATE_FACTORS=[]
MDA_ON=False
MDA_START_TIMES=[250,300]
HOST_FAIL_RATE=[0.2,0.2] # % of hosts that did not take the drug
DRUG_EFF_DURATION=[15,15] # How long the drug remains effective in the body
MDA_IMMIGRATION_RATE_FACTORS=[0.2,0.6]
IMMIGRATION_ON=True
IMMIGRATION_RATE=[1.0]
P_IMMIGRATION_INCLUDES_NEW_GENES=0.5
N_IMMIGRATION_NEW_GENES=0
