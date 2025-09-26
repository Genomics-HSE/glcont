sbatch -c 16 --wrap "sh run_sim_parallel.sh pairs_0_20.samples.txt"
sbatch -c 16 --wrap "sh run_sim_parallel.sh pairs_21_30.samples.txt"
sbatch -c 16 --wrap "sh run_sim_parallel.sh pairs_31_50.samples.txt"
