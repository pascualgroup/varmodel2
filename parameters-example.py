RANDOM_SEED=2800712
N_GENES_INITIAL=1200
N_GENES_PER_STRAIN=60
N_LOCI=2
N_ALLELES_INITIAL=N_LOCI*[120]
SELECTION_MODE='SPECIFIC_IMMUNITY'
N_INFECTIONS_FOR_GENERAL_IMMUNITY=0
GENERAL_IMMUNITY_PARAMS=[-1, -1, -1, -1]
CLEARANCE_RATE_IMMUNE=0.041628833086078301 # For generalized immunity
TRANSITION_RATE_NOT_IMMUNE=1.0 / 6.0 # switching rate. every 6 days
TRANSITION_RATE_IMMUNE=1000
PRINT_FUNCTION_TRACE=False
PRINT_DEBUG_LEVEL=0
T_YEAR=360.0
T_BURNIN=15
T_END=35000
SAMPLE_DB_FILENAME='""'
PRINT_INFO_PERIOD=120.0
VERIFICATION_ON=False
VERIFICATION_PERIOD=35000
SAVE_TO_CHECKPOINT=False
CHECKPOINT_SAVE_FILENAME='"test2cp.sqlite"'
CHECKPOINT_SAVE_PERIOD=0
LOAD_FROM_CHECKPOINT=True
CHECKPOINT_LOAD_FILENAME='"/Users/Qixin/Dropbox/varmodel2/build/PS36_S_E000_R1_CP.sqlite"'
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
BITING_RATE_MEAN=[0.0005]
BITING_RATE_RELATIVE_AMPLITUDE=[0]
BITING_RATE_PEAK_PHASE=[0]
DAILY_BITING_RATE_DISTRIBUTION=[184.982060636238,200.255404185268,190.091657432075,176.936673923227,164.343031052803,152.793760078616,142.280309311918,132.721990757893,124.033625709693,116.13590946434,108.956491569867,102.429591067436,96.4954777083097,91.0998775063187,86.1935081934926,81.7316148582299,77.6735473752113,73.9823941296359,70.624631409777,67.569823745941,64.7903263402798,62.261040675501,59.9591757095815,57.8640360180417,55.9568281889751,54.2204871570074,52.6395139643939,51.1998352707491,49.8886665452749,48.6943983409682,47.6064821131221,46.6153343553555,45.7122530802748,44.889322061416,44.13935303385,43.4558081247075,42.8327483785961,42.2647680297853,41.7469538312647,41.2748368576276,40.8443501562017,40.4517970480741,40.0938070759306,39.7673175268291,39.4695383049637,39.1979296104896,38.9501783120656,38.7241775187336,38.5180090184007,38.3299237028482,38.1583293639103,38.0017737769016,37.8589341338113,37.7286047949206,37.609686843514,37.50117832764,37.4021653181591,37.3118147288408,37.2293672736381,37.1541298330001,37.0854713133807,37.0228146684226,36.9656352249567,36.9134526965649,36.8658304874747,36.8223689332749,36.7827047634528,36.7465053716365,36.7134680370484,36.683316586589,36.6557984022823,36.6306826034548,36.6077597837153,36.586838896635,36.5677447120854,36.5503175349244,36.5344117986761,36.5198946799105,36.5066448234627,36.4945517554252,36.4835144475775,36.4734406657602,36.4642461883492,36.4558541345434,36.4481944766219,36.4412033547379,36.4348224858472,36.4289985845382,36.423682949872,36.4188312037933,36.4144029157693,36.4103610895911,36.4066719764741,36.4033048155977,36.4002316266374,36.3974265546817,36.3948664137349,36.3925295614191,36.390396747356,36.3884499478837,36.3866731020119,36.3850512618356,36.383570977466,36.3822198706082,36.3809866531596,36.3798611099218,36.3788337913255,36.3778961228354,36.3770403186877,36.376259127581,36.3755461792137,36.3748953467457,36.3743013727248,36.3737591652055,36.3732643065013,36.3728126044509,36.3724003301053,36.3720240316874,36.3716805650306,36.3713671195684,36.3710809645286,36.3708198308674,36.3705814728698,36.3703638797406,36.370165316919,36.3699840324805,36.369818578311,36.3696675718257,36.3695297102287,36.3694039040808,36.3692890594064,36.3691842391427,36.3690885778915,36.3690012449681,36.368921567071,36.3688488589252,36.3687824611981,36.3687218882497,36.3686666214112,36.368616137142,36.368570066055,36.3685280375853,36.3684896516034,36.3684545980557,36.3684226195349,36.3683934248545,36.3683667563531,36.3683424240487,36.3683202157831,36.3682999305956,36.368281422106,36.3682645324877,36.368249105087,36.3682350279662,36.3682221849002,36.3682104538153,36.3681997446467,36.3681899755132,36.3681810545177,36.3681729107047,36.3851737913722,36.6273491430708,37.4853628842298,39.3391521070744,42.5139152093929,47.287707648272,53.9085580483502,62.6109722095744,73.6276001726413,87.1941215922183,103.547494101901,122.919379658258,145.527298129489,171.565847125294,201.199569715256,234.558182047779,271.734187710335,312.782500620207,357.721545950115,406.535310033193,459.175902687809,515.566288170545,575.603037004223,639.158875127664,706.085079258744,776.213619996515,849.359092975813,925.320508432896,1003.88286873039,1084.81866407249,1167.88921717579,1252.84603432857,1339.43197501728,1427.38246518551,1516.42665703923,1606.288540236,1696.68804794844,1787.34213872482,1877.96586367714,1968.27341928461,2057.97918561558,2146.79874988722,2234.44991407981,2320.65368564903,2405.13524893337,2487.62491634739,2567.859055995,2645.58099505399,2720.54189505184,2792.50159831076,2861.22944219628,2926.50503922685,2988.11902160115,3045.87374666899,3099.58396234427,3149.0774305219,3194.19550599262,3234.79366968277,3270.74201453417,3301.92568255082,3328.24525171735,3349.61707166055,3365.97354703317,3377.26336782105,3383.45168589349,3384.52023734929,3380.46741043488,3371.30825831918,3357.07445753911,3337.81421221584,3313.59210333082,3284.48888537285,3250.60122960876,3212.04141547099,3168.93697156399,3121.4302662414,3069.67805102919,3013.85095644985,2954.13294358876,2890.72071254956,2823.82306974722,2753.66025685211,2680.46324274245,2604.472983022,2525.93965668524,2445.12188176186,2362.28573918503,2277.70407924451,2191.65575876302,2104.42450326994,2016.29812452819,1927.56762428701,1838.52625981138,1749.46862102118,1660.68973822427,1572.48417420318,1485.14512218572,1398.96353218762,1314.22725642948,1231.22022187782,1150.22162457876,1071.50518305081,995.338420555275,921.981996241624,851.689079603363,784.704748808794,721.265363911465,661.597778657233,605.918111912624,554.429624489276,507.29104204215,464.393290873619,425.407641599444,389.985626692259,357.802524436665,328.562007641641,301.994626224167,277.855605119689,255.922545136195,235.993412158679,217.884657199491,201.429550269167,186.476615850599,172.888251531252,160.539461254969,149.316695756803,139.116797656341,129.846076373408,121.419412190101,113.759501549716,106.796114867289,100.465476414138,94.709648242056,89.4760212533734,84.7168040863485,80.3886048480315,76.4520068339136,72.8712249533315,69.6137536857966,66.6500844545312,63.9534135854316,61.4994098098765,59.2659720392144,57.2330380312933,55.3823836622487,53.69746513125,52.1632546823583,50.7661037923442,49.4936221913513,48.3345467621114,47.2786509575201,46.316642834647,45.4400794217592,44.6412870359932,43.9132945674952,43.2497596404901,42.6449197270064,42.0935343856084,41.5908356082067,41.1324876618248,40.714549150696,40.3334264715257,39.9858542704671,39.6688594745203,39.379735021602,39.1160158594602,38.8754565646507,38.6560144266936,38.455825228768,38.2731932483092,38.1065711198888,37.9545499808064,37.8158473626304,37.6892902381191,37.5738142292774,37.4684449625939,37.3722953501391,37.2845577495937,37.2044940467375,37.1314316812032,37.0647574694356,37.0039118006962,36.948384199192,36.8977096264088,36.8514623263637,36.8092571039936,36.770737875963,36.7355843878313,36.7035008724317,36.6742197285944,36.6474959600555,36.6231057548678,36.6008445488391,36.5805270876016,36.5619839390151,36.5450597868185,36.5296131069445,36.5155148897651,36.5026474377478,36.4909032774043]
IRS_ON=True
IRS_START_TIMES=[28950]
BITING_RATE_FACTORS=[[184.982060636238,200.255404185268,190.091657432075,176.936673923228,164.343031052801,152.793760078616,142.280309311919,132.721990757895,124.033625709624,116.135909464249,108.956491569859,102.429591067481,96.495477708316,91.09987750652,86.1935081935306,81.7316148582278,77.6735473751643,73.982394129563,70.6246314093291,67.5698237458452,64.7903263401945,62.2610406754092,59.959175709491,57.8640360180816,55.9568281889067,54.2204871570316,52.6395139631895,51.1998352706814,49.8886665452167,48.6943983409033,47.6064821130883,46.6153343546024,45.7122530802607,44.8893220613562,44.1393530338366,43.4558081239425,42.8327483785986,42.2647680298509,41.7469538312049,41.2748368576288,40.8443501571467,40.4517970485615,40.0938070754175,39.7673175269324,39.4695383050475,39.1979296104861,38.9501783120655,38.7241775188069,38.5180090180686,38.3299237023424,38.1583293638247,38.0017737777167,37.8589341342821,37.728604790006,37.6096867825927,37.5011783298152,37.4021652598463,37.3118146846574,37.2293669768166,37.1541296114915,37.0854711505477,37.0228149555499,36.9656353525861,36.913452979037,36.8658306163296,36.8223691638958,36.7827049309392,36.7465054185517,36.7134679470576,36.6833160910902,36.6557977411604,36.6306826484332,36.6077606529641,36.5868393540273,36.5677448788829,36.550317768031,36.5344121322457,36.5198948627396,36.5066448312893,36.4945514532639,36.4835137623846,36.4734395880891,36.4642447904162,36.4558525662252,36.4481928471122,36.4412017067385,36.4348207200189,36.4289967538628,36.423681251215,36.418829782654,36.4144018931781,36.410360410137,36.4066714668058,36.4033044595958,36.400231142602,36.3974260524796,36.3948658142862,36.3925290236655,36.3903962884307,36.3884495769372,36.3866728789688,36.3850510765042,36.383570881991,36.3822197469098,36.3809865721771,36.379860955045,36.3788336216379,36.3778958948846,36.3770401233274,36.3762589913056,36.3755460643695,36.3748954189662,36.3743014212894,36.3737593690929,36.3732645236708,36.3728128552858,36.3724006250497,36.3720242862518,36.371680841067,36.3713673356634,17.2697952096329,9.39196131291392,6.10972886177342,4.70727310593465,4.07482548544797,3.76076165964537,3.58241727162781,3.46625281852498,3.38242651641713,3.31829469707124,3.26789437552109,3.22788854286972,3.19605375224608,3.1707272375712,3.15059791346112,3.1346144983667,3.12193240596048,3.1118749561071,3.10390168394546,3.09758212127039,3.09257400145883,3.0886055474542,3.08546114242121,3.0829697669556,3.08099581236387,3.07943193040232,3.07819290603465,3.0772112725815,3.0764335240385,3.07581735672273,3.07532920059806,3.07494244726987,3.07463601547071,3.07439328660884,3.07420099426183,3.07404864122134,3.0739279372193,3.07383230477958,3.07375654224199,3.07369652416994,3.07769395218549,3.1414533008597,3.35806373687594,3.77223470313668,4.38467625913802,5.18647966875562,6.17887945685692,7.37742780578422,8.80900000190005,10.5074837774453,12.5103891822764,14.8565391658069,17.5843505669733,20.7302612326687,24.3271037296742,28.40243156225,32.9769202251286,38.0630085019713,43.6639323245414,49.7732798680477,56.3750943913381,63.4445305935157,70.94896939707,78.8494349777068,87.1022092929547,95.6604105121745,104.475496583282,113.498533640787,122.681214668054,131.976609726492,141.339671975874,150.727514781171,160.099535667819,169.417412833078,178.64502111292,187.748293425694,196.695063764887,205.454917082164,213.999049814447,222.300131175919,230.332211063558,238.070634881311,245.491994325881,252.574108947918,259.295979071357,265.637814605436,271.581036070854,277.108306768877,282.203557384218,286.852027441008,291.040300681182,294.756341980308,297.989546088764,300.730767169596,302.972354062101,304.708189351719,305.933707692064,306.645921692841,306.843442814603,306.526487851299,305.696889815546,304.358095224473,302.515166941085,300.17476058871,297.345121420922,294.036058324874,290.258912280739,286.026524769175,281.35320127989,276.254682970789,270.748065428004,264.851775411959,258.585496802312,251.970124567997,245.027681979472,237.781262827156,230.254950331981,222.473751125494,214.463508965274,206.250829267724,197.862982918402,189.327833268519,180.673745973972,171.929489527161,163.124163919476,154.287095655281,145.447751653483,136.635653887765,127.880284833197,119.211003508537,110.656959540389,102.247011256504,94.0096471504235,85.9729126413586,78.1643427128047,70.6109026702607,63.3389385010651,56.3741404726545,49.7415233510851,43.4654280711278,37.5695526261245,32.0770204031651,27.010498011497,22.3923724013311,18.2449956389214,14.5909628176001,11.4532868748982,8.85501385253545,6.81706251472276,5.3519080073952,4.4426775957753,3.95269806888622,3.68907590491369,3.53308195166586,3.42921251952119,3.35347098612022,3.29530175008022,3.24955606026402,3.21326329196119,3.18440586266922,3.16146433903613,3.14324048899571,3.1287756181447,3.11730138075282,3.10820333192519,3.10099144008778,3.09527573612939,3.09074638947321,3.08715743569868,3.08431370010734,3.08206066829187,3.08027558752257,3.07886131770865,3.07774080109642,3.07685305463296,3.07614973819868,3.07559250861233,3.07515105004917,3.07480138883057,3.07452425636413,5.39939745079435,7.5427980350675,9.57653997253084,11.5384597040887,13.4274765978439,15.2280133019628,16.9264216139503,18.5153432637351,19.9929218319564,21.3610551621119,22.6239342168921,23.7870303790438,24.8564377104292,25.8384601464683,26.739358396113,27.565199511139,28.3217709787534,29.0145355210663,29.6486119583209,30.2287714083495,30.7594432216491,31.244727321492,31.6884109195993,32.0939845024455,32.4646627221379,32.8033981615401,33.1129072573759,33.3956785577724,33.6539991105958,33.8899609202755,34.1054843661704,34.3023244464606,34.4820913764015,34.6462566043052,34.7961667484731,34.9330530769389,35.0580433353263,35.1721661006144,35.2763637048275,35.3714962006854,35.4583497787315,35.5376433063482,35.6100332607399,35.6761186578941,35.7364484986768,35.7915228243772,35.8417986436087,35.8876934410222,35.9295885206493,35.9678324445515,36.0027414918681,36.0346080792957,36.0636959088595,36.0902477747719,36.1144840571803,36.1366067945173,36.1568004997859,36.1752336731943,36.1920589484922,36.2074164093063,36.2214343166083,36.2342291996656,36.2459078875182,36.2565676368486,36.2662977928293,36.2751789428576,36.283285351754,36.2906844377947,36.2974379341557,36.3036027549186,36.3092298380934,36.3143657934768,36.3190534867752,36.323332137481,36.3272374822276,36.3308021605487,36.3340557618111,36.3370254601239,36.3397359281382,36.3422097555697,36.3444677701747,36.346528735797,36.3484099163099,36.3501270181968,36.3516942755006,36.3531248893912,36.3544305038186,36.3556222092393,36.3567098443219,36.3577024820841,36.3586086196342,36.3594354697369,36.3601903420851,36.3608793272087,36.3615081068508,36.3620822265273,36.3626061645142,36.3630844126565,36.363520987855,36.3639193990332,36.3642831092652,36.3646150321236,36.3649179900056,36.3651945489402,36.365446913549,36.3656773213366,36.3658876070168,36.3660795282115,36.3662547569169,36.3664146640794,36.3665606353301,36.3666938757998,36.3668154712902,36.3669264800648,36.3670277824717,36.3671202516154,36.367204649343,36.3672816433555,36.3673519222192,36.3674160935183,36.3674746217306,36.3675280376187,36.3675768394854,36.3676213571142,36.3676619757138,36.3676990871535,36.3677329640273,36.3677638711049,36.3677920956757,36.3678178732229,36.3678413913377,36.3678628588352,36.3678824668462,36.3679003550672,36.3679166793435,36.367931589919,36.3679451845769,36.3679575769521,36.3679688874272,36.3679792355777,36.3679886753144,36.3679972655248,36.368005099468,36.3680122704028,36.3680188402327,36.3680248171383,36.3680302625223,36.3680352400722,36.3680398078233,36.3680439753741,36.3680477726617,36.3680512381671,36.3680544102865,36.3680573101312,36.368059953281,36.3680623640281,36.3680645666651,36.3680665804797,36.3680684138318,36.3680700844734,36.3680716107491,36.3680730088859,36.368074275187,36.3680754198126,36.3680764581623,36.3680774056357,36.3680782776321,36.3680790890484,36.368079824877,36.3680804811109,36.3680810710449,36.3680816079739,36.3680821051925,36.3680825759956,36.3680830147514,36.3680834048103,36.3680837545401,36.3680840725102,36.3680843672898,36.3680846474484,17.2684579157822,9.39140229675402,6.10949016607538,4.70716638650129,4.07477352318296,3.7607329632542,3.58239907197854,3.46623991757722,3.382416734516,3.31828703448764,3.26788829636959,3.22788370304217,3.19604989914545,3.17072417327415,3.15059547918904,3.1346125661922,3.12193087333724,3.11187374090934,3.1039007206971,3.09758135788537,3.09257339652845,3.08860506813412,3.08546076261593,3.08296946604558,3.08099557395487,3.07943174151294,3.07819275637893,3.07721115401072,3.07643343011304,3.07581728231007,3.07532914163495,3.07494240057377,3.07463597841273,3.07439325721441,3.07420097102455,3.07404862281231,3.07392792264219,3.07383229323919,3.07375653308905,3.07369651702775,3.56740404670625,5.27368512632685,7.63395909554388,10.4033133707667,13.7132496456331,17.8269420094311,23.0193036445392,29.5501573190855,37.66370382497,47.5793934121649,59.4653004938241,73.3980358191748,89.3193075672435,107.002848066263,126.046344154312,145.898512213887,165.920896251031,185.470550440752,203.980958464054,221.019660155333,236.31151189727,249.729937685316,261.267444412196,270.998730172022,279.04635015526,285.553967207728,290.668189395351,294.527793643036,297.258301113687,298.97008770424,299.75856977228,299.705480282341,298.880632483498,297.343780109252,295.146411229424,292.333338533457,288.944121169385,285.014265032566,280.576179247222,275.660032697379,270.294387040926,264.506737647013,258.323931546874,251.772492416424,244.878863916704,237.669593834846,230.171461162824,222.411565382622,214.417372977657,206.216737615636,197.837897887596,189.309458040103,180.66034594566,171.919767277854,163.117145014201,154.282055106699,145.444153399417,136.633100814219,127.878485181087,119.209743747121,110.656084335153,102.246408124083,94.0092352208327,85.972633913363,78.1641561599394,70.6107793069711,63.3388578881835,56.3740885209366,49.7414903993698,43.4654075619396,37.5695401308553,32.0770129460267,27.0104936776252,22.3923699501541,18.2449942946487,14.5909621065666,11.4532865087461,8.85501366197802,6.81706241262013,5.35190651512511,4.45121293965082,4.07231562996567,4.13735201299551,4.55235861120669,5.23476769881445,6.13222855111179,7.22632041776383,8.52561369343738,10.0555567401185,11.8500893810539,13.9463524269934,16.3818041645356,19.1925751302986,22.4122106538111,26.0704140171443,30.1917267459348,34.7942488343629,39.8885453567539,45.4768812133202,51.5528754367029,58.1016131454852,65.100201575583,72.5186970312611,80.3212983460249,88.467688924454,96.9143912768449,105.616050421983,114.526550298415,123.599943208845,132.791122076246,91.7679074094182,45.8319819458187,22.1517793931536,11.5947901077148,7.07660836641783,5.1421472925011,4.28509566125057,3.8751635842895,3.65382664073258,3.51629250588903,3.42014162639874,3.34776970089742,3.29126689093473,3.24650067493214,3.21087739577162,3.18252081689958,3.15996959143529,3.142054428746,3.12783471990606,3.11655522433494,3.10761178359949,3.10052257337219,3.09490417507722,3.0904519714558,3.08692407714441,3.08412889839465,3.08191423392779,3.08015958548557,3.07876936972556,3.07766795931623,5.40399718053844,7.54832862406546,9.5825419156512,11.5445710602633,13.4334584979913,15.2337308596365,16.931808351523,18.5203728027436,19.9975904025072,21.3653715393157,22.627913950535,23.7906924361,24.859802500713,25.8415483580957,26.7421903136956,27.5677946260601,28.3241477777106,29.0167113834799,29.6506031196573,30.2305929706354,30.7611091821734,31.2462506206399,31.689803497517,32.0952573569916,32.4658259696974,32.8044611026656,33.1138784264057,33.39656578922,33.6548095831673,33.8907012170213,34.1061605131531,34.3029419647307,34.4826553138664,34.6467715877999,34.7966370047588,34.9334824851288,35.0584349359445,35.1725240992778,35.2766898687721,35.371793862809,35.4586215500229,35.5378914033311,35.6102599154886,35.6763256852797,35.7366373415926,35.791695097616,35.8419559261741,35.8878368459925,35.9297195101969,35.9679513746967,36.0028506685056,36.0347073804,36.0637868791858,36.0903308570288,36.1145601652005,36.1366768110492,36.1568643352708,36.1752913734978,36.1921113895996,36.2074639907384]]
IRS_IMMIGRATION_RATE_FACTORS=[0]
MDA_ON=False
MDA_START_TIMES=[]
HOST_FAIL_RATE=[] # % of hosts that did not take the drug
DRUG_EFF_DURATION=[] # How long the drug remains effective in the body
MDA_IMMIGRATION_RATE_FACTORS=[]
POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION = False
IMMIGRATION_ON=True
IMMIGRATION_RATE=[1]
P_IMMIGRATION_INCLUDES_NEW_GENES=0.5
N_IMMIGRATION_NEW_GENES=0
