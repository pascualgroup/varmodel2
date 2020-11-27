#!/bin/bash
#SBATCH --job-name=runtimeTest
#SBATCH --time=36:00:00
#SBATCH --output=/scratch/midway2/heqixin/runInfo/runtimeTest_%A_%a.out
#SBATCH --error=/scratch/midway2/heqixin/runInfo/runtimeTest_%A_%a.err
#SBATCH --array=1-10
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30000
#SBATCH --partition=broadwl
#SBATCH --account=pi-pascualmm
#SBATCH --mail-type=ALL
#SBATCH --mail-user=heqixin@uchicago.edu

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
module load gcc/6.1
module load python/3.7.0
module load R

filePrefix=runtimeTest

#here define where the actual run will be stored
remoteDir=/scratch/midway2/heqixin/${filePrefix}_${SLURM_ARRAY_TASK_ID}
#copy original Rdiv code to the running folder
cp -r /home/heqixin/varmodel2 $remoteDir

cd $remoteDir
cd $filePrefix
python writeParameters.py -p ${filePrefix}_param.csv -i ${filePrefix}_Template.py -n $SLURM_ARRAY_TASK_ID -x $filePrefix
cd ..
#build the model, run preIRS first
./build.py -p $filePrefix/${filePrefix}_${SLURM_ARRAY_TASK_ID}_input.py -d n_{$SLURM_ARRAY_TASK_ID}
#execute the run
./n_{$SLURM_ARRAY_TASK_ID}/bin/varMig
