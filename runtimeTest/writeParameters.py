#require python 3.6 or above
# usage
# python writeParameters.py -i <paramExampleFile> -p <paraList> -n <run number> -x <PREFIX>
import csv,random
from optparse import OptionParser
import sqlite3 as lite
import numpy as np
import sys

parser = OptionParser()
parser.add_option("-p", "--parameters", dest="paraList",\
					help="csv file for all the parameter combinations", metavar="FILE")
parser.add_option("-i", "--input", dest="paramExampleFile",\
					help="input template file", metavar="FILE")
parser.add_option("-n", "--number", dest="outNumber",type = "int",\
					help="the run number")
parser.add_option("-x", "--prefix", dest="prefix",type = "string",\
					help="prefix for output filenames")
                  
(options, args) = parser.parse_args()


if __name__ == "__main__":
	#step 1 read in parameter combination file
	allParamFile = open(options.paraList,newline = '')
	allParam = csv.DictReader(allParamFile)
	
	#step 2 read in input file template
	prototype = open(options.paramExampleFile,"r").read()
	
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
			
			row['SAMPLE_DB_FILENAME'] = options.prefix + "_" + row['NO'] + "_sd.sqlite"
			
			out = open(options.prefix + "_" + row['NO'] + "_input.py", "w")
			out.write(prototype.format(**row))
			out.close()
			
			break
	

