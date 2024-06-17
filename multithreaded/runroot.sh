
spack load gcc@8.3.0%gcc@8.3.0
spack load root
#qsub -q i12h -F "$1" root.sh 
sbatch -p INTEL_HASWELL --time=10:00:00 --cpus-per-task=40 --nodes=1 --job-name="fit2Dfull" --mail-type=END --mail-user=adam.szabelski@ncbj.gov.pl --output="selection" root.sh 
