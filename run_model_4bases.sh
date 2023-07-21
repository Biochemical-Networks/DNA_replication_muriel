#PBS -N ATiGiC_rep_equi
#PBS -q highcpu 
#PBS -l walltime=156:00:00
#PBS -lselect=1:ncpus=1:mem=4gb
#PBS -lplace=pack
#PBS -j oe
#PBS -o /home/ipausers/louman/job_outs_eo


echo Variables defined for this job

#run Jennys model having these parameters, NB: you should always first compile the code in the right directory
python3 Documents/programming/DNA_replication_muriel/230413_read_parameters_4bases_forloop.py



