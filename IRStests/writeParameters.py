#require python 3.6 or above
# usage
# python writeParameters.py -i <paramExampleFile> -p <paraList> -n <run number>
import csv,random
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--parameters", dest="paraList",\
					help="csv file for all the parameter combinations", metavar="FILE")
parser.add_option("-i", "--input", dest="paramExampleFile",\
					help="input template file", metavar="FILE")
parser.add_option("-n", "--number", dest="outNumber",type = "int",\
					help="the run number")
parser.add_option("-r", "--rep", dest="repNumber",type = "int",\
					help="define how many replications of control is run")
parser.add_option("-x", "--prefix", dest="prefix",type = "string",\
					help="prefix for output filenames")
                  
(options, args) = parser.parse_args()


if __name__ == "__main__":
	#step 1 read in parameter combination file
	allParamFile = open(options.paraList,newline = '')
	allParam = csv.DictReader(allParamFile)
	
	#step 2 read in input file template
	prototype = open(options.paramExampleFile,"r").read()
	
	#0 is control, 1,2,3 corresponds to 2,5,and 10yr control
	scenario = [0,1,2,3]
	
	#iterate through the list and find the correct run number
	for row in allParam:
		#print(row)
		if row['NO'] == str(options.outNumber):
			f = open(row['DAILY_BITING_RATE_DISTRIBUTION'], newline = '')
			db = csv.reader(f)
			dbVec = []
			for y in db:
				dbVec.append(y[1])
			
			row['DAILY_BITING_RATE_DISTRIBUTION'] = ','.join(dbVec)
			
			t_end = row['T_END']
			
			#write scenario 0 control first
			for sc in scenario:
				if sc == 0:
					row['IRS'] = 'False'
					row['SAMPLE_DB_FILENAME'] = options.prefix + "_" + row['NO'] + "_s" + str(sc) + "_sd.sqlite"
					row['SAVE_TO_CHECKPOINT'] = 'True'
					row['CHECKPOINT_SAVE_FILENAME'] = options.prefix + "_" + row['NO'] + "_s" + str(sc) + "_cp.sqlite"
					row['CHECKPOINT_LOAD_FILENAME'] = ''
					row['T_END']=row['IRS_START']
					row['BITING_RATE_FACTORS']=''
					out = open(options.prefix + "_" + row['NO'] + "_s" + str(sc) +  "_input.py", "w")
					out.write(prototype.format(**row))
					out.close()
				else:
					row['IRS'] = 'True'
					row['T_BURNIN'] = '0'
					row['T_END'] = t_end
					f2 = open(row['IRS_BITING']+str(sc)+".csv", newline = '')
					db2 = csv.reader(f2)
					db2Vec = []
					for y in db2:
						db2Vec.append(y[1])
					
					row['BITING_RATE_FACTORS']=','.join(db2Vec)
					
					for r in range(options.repNumber):
						row['RANDOM_SEED'] = random.randint(1,10000000000)
						row['SAMPLE_DB_FILENAME'] = options.prefix + "_" + row['NO'] + "_s" + str(sc) + "_r" + str(r) +  "_sd.sqlite" 
						row['SAVE_TO_CHECKPOINT'] = 'False'
						row['CHECKPOINT_SAVE_FILENAME'] = ''
						row['CHECKPOINT_LOAD_FILENAME'] = options.prefix + "_" + row['NO'] + "_s0_cp.sqlite"
						out = open(options.prefix + "_" + row['NO'] + "_s" + str(sc) + "_r" + str(r) + "_input.py", "w")
						out.write(prototype.format(**row))
						out.close()
			
			
			break
	

