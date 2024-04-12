
spack load gcc@9.3.0%gcc@9.3.0
spack load root
#qsub -q i12h -F "$1" root.sh 
sbatch -p INTEL_HASWELL --time=120:00:00 --cpus-per-task=1 --nodes=1 --job-name="fit2Dfull" --mail-type=END --mail-user=adam.szabelski@ncbj.gov.pl --output="selection$1" root.sh $1 
