#require python 3.6 or above
# usage
# python writeParameters.py -i <paramExampleFile> -p <paraList> -n <run number>
import csv,random
from optparse import OptionParser
import sqlite3 as lite
import numpy as np
import sys
from scipy import stats
from scipy.optimize import curve_fit

parser = OptionParser()
parser.add_option("-p", "--parameters", dest="paraList",\
					help="csv file for all the parameter combinations", metavar="FILE")
parser.add_option("-i", "--input", dest="paramExampleFile",\
					help="input template file", metavar="FILE")
parser.add_option("-n", "--number", dest="outNumber",type = "int",\
					help="the run number")
parser.add_option("-g", "--GI", dest="generalImm",action = "store_true",\
					help="whether output generalized immunity parameters", default = False)
parser.add_option("-r", "--rep", dest="repNumber",type = "int",\
					help="define how many replications of control is run")
parser.add_option("-x", "--prefix", dest="prefix",type = "string",\
					help="prefix for output filenames")
                  
(options, args) = parser.parse_args()


##fitting function
def func(x,a,b,c,d):
	return b*np.exp(-c*x)/((d*x+1)**d)+a


def GIParam(fname):
	con = lite.connect(fname)
	with con:
		cur = con.cursor()
		cur.execute("SELECT duration, infection_id FROM sampled_duration WHERE (time-duration) > 5000")
		rows = cur.fetchall()
			
	xd = np.array([x[1] for x in rows])-1  #number of infections
	yd = np.array([x[0] for x in rows])  #infection duration
	
	#get range of infection times
	newxd=range(xd.min(),xd.max()+2)
	# mean of infection duration as a function of infection times
	xys = stats.binned_statistic(xd,yd,'mean',bins=newxd)
	bc = np.bincount(xd)
	flt = np.where(bc>1)  #take only bins that have more than 1 incidents
	x = flt[0]
	y = xys[0][flt[0]]-14 
	if len(x)>500:
		x = x[1:500]
		y = y[1:500]
	
	infectionTimesToImmune = max(x)
	if len(yd[xd>max(x)])>0:
		clearanceRateConstantImmune = 1/(np.mean(yd[xd>max(x)]) - 14)
	else:
		clearanceRateConstantImmune = 1/y[-1]
	
	while True:
		try:
			popt,pcov = curve_fit(func,x,y,p0 = (0.01,50,0.0017,0.8))
			break
		except:
			pass
			
	generalImmunityParams = list(popt)
	return (infectionTimesToImmune,clearanceRateConstantImmune,generalImmunityParams)


if __name__ == "__main__":
	#step 1 read in parameter combination file
	allParamFile = open(options.paraList,newline = '')
	allParam = csv.DictReader(allParamFile)
	
	#step 2 read in input file template
	prototype = open(options.paramExampleFile,"r").read()
	
	#0 is control, 1,2,3 corresponds to 2,5,and 10yr control
	if options.generalImm:
		scenario = [4,5,6,7]
	else:
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
			
			if min(scenario) >3:
				row['SELECTION_MODE'] = 'GENERAL_IMMUNITY'
				#generalized immunity
				SIfname = "../" + options.prefix + "_" + row['NO'] + "_s0_sd.sqlite"
				row['N_INFECTIONS_FOR_GENERAL_IMMUNITY'],row['CLEARANCE_RATE_IMMUNE'],general_imm_param = GIParam(SIfname)
				row['GENERAL_IMMUNITY_PARAMS'] = ','.join([str(x) for x in general_imm_param])
				row['CHECKPOINT_LOAD_FILENAME'] = options.prefix + "_" + row['NO'] + "_s4_cp.sqlite"
			else:
				row['CHECKPOINT_LOAD_FILENAME'] = options.prefix + "_" + row['NO'] + "_s0_cp.sqlite"
				row['N_INFECTIONS_FOR_GENERAL_IMMUNITY'] = 0
				row['CLEARANCE_RATE_IMMUNE'] = 1.0
				row['GENERAL_IMMUNITY_PARAMS'] = ''
							
			#write scenario 0 control first
			for sc in scenario:
				if sc == 0 or sc == 4:
					row['IRS'] = 'False'
					row['SAMPLE_DB_FILENAME'] = options.prefix + "_" + row['NO'] + "_s" + str(sc) + "_sd.sqlite"
					row['SAVE_TO_CHECKPOINT'] = 'True'
					row['CHECKPOINT_SAVE_FILENAME'] = options.prefix + "_" + row['NO'] + "_s" + str(sc) + "_cp.sqlite"
					#row['CHECKPOINT_LOAD_FILENAME'] = ''
					row['T_END']=int(row['EXPECTED_EQUILIBRIUM'])+720
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
						out = open(options.prefix + "_" + row['NO'] + "_s" + str(sc) + "_r" + str(r) + "_input.py", "w")
						out.write(prototype.format(**row))
						out.close()
			
			break
	

