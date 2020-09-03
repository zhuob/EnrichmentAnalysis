#!/bin/sh

# Give the job a name
#$ -N simulation

#$ -S /bin/sh

# set working directory on all host to
# directory where the job was started
#$ -cwd

# send output to job.log (STDOUT + STDERR)
# #$ -o Enrichment_testStatistic.txt
# #$ -o SimulationLot.txt
# #$ -o emprical_power.txt
#$ -o meacaResult_test.txt
#$ -j y



# to request 100G free memory
#$ -l mem_free=100G


# email information
#$ -m e
# Just change the email address.  You will be emailed when the job has finished.
#$ -M rickyb.zh@gmail.com

#don't change the source
source /etc/profile.d/modules.sh

# USE STATISTICS DEPARTMENT CLUSTER
#$ -q di.q@di003 
# #$ -q all.q@cosine003
# #$ -q shared.q@di003

#Change which version of R you want to load on the Compute Nodes
module load R/3.6.3
module unload gcc/5.1.0
module load gcc/7.3.0

# command to run.  ONLY CHANGE THE NAME OF YOUR APPLICATION  
# Rscript real_data_correlation_simulation.R
Rscript run-simulation.R
# Rscript padog-real-data-type1error.R
