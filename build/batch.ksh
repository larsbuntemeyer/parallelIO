#!/bin/ksh
#
#SBATCH --job-name="parallelIO"      # Specify job name
#SBATCH --partition=compute,compute2
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --mail-user=lars.buntemeyer@hzg.de # Set your e−mail address
#SBATCH --account=ch0636       # Charge resources on this project account

# in bash or ksh script
source /sw/rhel6-x64/etc/profile.mistral

# Environment settings to run a MPI parallel program compiled with BullxMPI and Mellanox libraries
# Load environment
module load allinea-forge
module load intel/15.0.3
module load mxm/3.3.3002
module load fca/2.5.2393
module load bullxmpi_mlx/bullxmpi_mlx-1.2.8.3
# Settings for Open MPI and MXM (MellanoX Messaging) library
export OMPI_MCA_pml=cm
export OMPI_MCA_mtl=mxm
export OMPI_MCA_mtl_mxm_np=0
export MXM_RDMA_PORTS=mlx5_0:1
export MXM_LOG_LEVEL=ERROR
# Disable GHC algorithm for collective communication
export OMPI_MCA_coll=^ghc
# limit stacksize ... adjust to your programs need
ulimit -s 102400

# Environment settings to run a MPI parallel program compiled with Intel MPI
#module load allinea-forge
#module load intelmpi
#export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

export NETCDFF=/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-parallel-bullxmpi-intel14
export NETCDFC=/sw/rhel6-x64/netcdf/netcdf_c-4.4.0-parallel-bullxmpi-intel14
export LD_LIBRARY_PATH=${NETCDFC}/lib:${NETCDFF}/lib:${LD_LIBRARY_PATH}

# Use srun (not mpirun or mpiexec) command to launch programs compiled with any MPI library
#map --profile srun -l --propagate=STACK --cpu_bind=cores --distribution=block:cyclic main
#map --profile srun --cpu_bind=cores --distribution=block:cyclic main
#time srun -l --propagate=STACK --cpu_bind=cores --distribution=block:cyclic main
time srun ./main 
