#!/usr/bin/env bash
#
#
# ==== SLURM part (resource manager part) ===== #
#
# > Modify the following options based on your job's needs.
#   Remember that better job specifications mean better usage of resources,
#   which then means less time waiting for your job to start.
#   So, please specify as many details as possible.
#   A description of each option is available next to it.
#   SLURM cheatsheet:
# 
#     https://slurm.schedmd.com/pdfs/summary.pdf
# 
#
# ---- Metadata configuration ----
#
#SBATCH --job-name=KMH.dmft                 # The name of your job, you'll se it in squeue.
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL). Sends you an email when the job begins, ends, or fails; you can combine options.
#SBATCH --mail-user=gbellomi@sissa.it       # Where to send the mail
#
# ---- CPU resources configuration  ----    |  Clarifications at https://slurm.schedmd.com/mc_support.html
#
#SBATCH --ntasks=1                          # Number of MPI ranks (1 for MPI serial job)
#SBATCH --cpus-per-task=40                  # Number of threads per MPI rank (MAX: 2x32 cores on _partition_2, 2x20 cores on _partition_1) 
#[optional] #SBATCH --nodes=1               # Number of nodes
#[optional] #SBATCH --ntasks-per-node=1     # How many tasks on each node
#[optional] #SBATCH --ntasks-per-socket=1   # How many tasks on each socket
#[optional] #SBATCH --ntasks-per-core=1     # How many tasks on each core (set to 1 to be sure that different tasks run on different cores on multi-threaded systems)
#[optional] #SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically on nodes and sockets. For other options, read the docs.
#
# ---- Other resources configuration (e.g. GPU) ----
#
#[optional] #SBATCH --gpus=2                # Total number of GPUs for the job (MAX: 2 x number of nodes, only available on gpu1 and gpu2)
#[optional] #SBATCH --gpus-per-node=2       # Number of GPUs per node (MAX: 2, only available on gpu1 and gpu2)
#[optional] #SBATCH --gpus-per-task=1       # Number of GPUs per MPI rank (MAX: 2, only available on gpu1 and gpu2); to be used with --ntasks
#
# ---- Memory configuration ----
#
#SBATCH --mem=0                             # Memory per node (MAX: 63500 on the new ones, 40000 on the old ones); incompatible with --mem-per-cpu.
#[optional] #SBATCH --mem-per-cpu=4000mb    # Memory per thread; incompatible with --mem
#
# ---- Partition, Walltime and Output ----
#
#[unconfig] #SBATCH --array=01-10           # Create a job array. Useful for multiple, similar jobs. To use, read this: https://slurm.schedmd.com/job_array.html
#SBATCH --partition=regular1                # Partition (queue). Avail: regular1, regular2, long1, long2, wide1, wide2, gpu1, gpu2. Multiple partitions are possible.
#SBATCH --time=12:00:00                     # Time limit hrs:min:sec
#SBATCH --output=sLOG_%x_out%j.txt          # Standard output log -- WARNING: %x requires a new enough SLURM. Use %j for regular jobs and %A-%a for array jobs
#SBATCH --error=sLOG_%x_err%j.txt           # Standard error  log -- WARNING: %x requires a new enough SLURM. Use %j for regular jobs and %A-%a for array jobs
#
# ==== End of SLURM part (resource manager part) ===== #
#
#
# ==== Modules part (load all the modules) ===== #
#
# ---- ITCS-mantained modules ----
#
module load gnu8/8.3.0
module load mkl/19.1.3.304
module load openmpi3/3.1.4
module load matlab
#
# ---- QcmPlab stuff ----
#
module load scifor/gnu
module load dmft_tools/gnu
module load dmft_ed/gnu
#
# ==== End of Modules part (load all the modules) ===== #
#
#
# ==== Info part (say things) ===== #
#
# > DO NOT MODIFY. This part prints useful info on your output file.
#
NOW=`date +%H:%M-%a-%d/%b/%Y`
echo '------------------------------------------------------'
echo 'This job is allocated on '$SLURM_JOB_CPUS_PER_NODE' cpu(s)'
echo 'Job is running on node(s): '
echo  $SLURM_JOB_NODELIST
echo '------------------------------------------------------'
echo 'WORKINFO:'
echo 'SLURM: job starting at           '$NOW
echo 'SLURM: sbatch is running on      '$SLURM_SUBMIT_HOST
echo 'SLURM: executing on cluster      '$SLURM_CLUSTER_NAME
echo 'SLURM: executing on partition    '$SLURM_JOB_PARTITION
echo 'SLURM: working directory is      '$SLURM_SUBMIT_DIR
echo 'SLURM: current home directory is '$(getent passwd $SLURM_JOB_ACCOUNT | cut -d: -f6)
echo ""
echo 'JOBINFO:'
echo 'SLURM: job identifier is         '$SLURM_JOBID
echo 'SLURM: job name is               '$SLURM_JOB_NAME
echo ""
echo 'NODEINFO:'
echo 'SLURM: number of nodes is        '$SLURM_JOB_NUM_NODES
echo 'SLURM: number of cpus/node is    '$SLURM_JOB_CPUS_PER_NODE
echo 'SLURM: number of gpus/node is    '$SLURM_GPUS_PER_NODE
echo '------------------------------------------------------'
#
# ==== End of Info part (say things) ===== #
#
cd $SLURM_SUBMIT_DIR # Brings the shell into the directory from which you’ve submitted the script.
#
# ==== JOB COMMANDS ===== #
#
# > The part that actually executes all the operations you want to do.
#   Just fill this part as if it was a regular Bash script that you want to
#   run on your computer.
#
# >> DMFT-Workflow
#matlab -batch KMH-DMFT_dry		#-----------------
#matlab -batch KMH-DMFT_autostop	# Uncomment just
#matlab -batch KMH-DMFT_autoupdate	# one of these...
#matlab -batch KMH-DMFT_livemixing	#-----------------
# >> Post-Analysis
matlab -batch PostDMFT
mkdir postData_$SLURM_JOB_NAME
find . -name '*.mat' | cpio -pdm ./postData_$SLURM_JOB_NAME
#
#
# ==== END OF JOB COMMANDS ===== #
#
#
# Wait for processes, if any.
echo "Waiting for all the processes to finish..."
wait
