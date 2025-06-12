# Welcome to the weird and wonderful world of AREPO on COSMA!

This guide will help you avoid issues commonly encountered in setting up you first AREPO simulations.
It assumes that you have access to the AREPO repository (or that you have access to a local copy) 
and that you are compiling and running the code on the Distributed Research utilising Advanced 
Computing (DiRAC) facility's COSmological MAchine 8 (COSMA8).
This guide was written for galactic simulations, though the general setup can be followed for other 
use cases.

Good general advice is to familiarise yourself with the documentation folder of the AREPO repository:
the "getting_started.md" file will be particularly useful at first.

The main steps required to run and work with AREPO simulations are as follows:  
I. Get access to AREPO  
II. Compile AREPO  
III. Run AREPO  
IV. Read and analyse AREPO outputs  

These steps are detailed further below.

Getting access to AREPO requires either permissions to the bitbucket repository, a local copy of the 
repository or access to the public version of the code. If you have a local copy of the code, you can
skip section I.

Compiling AREPO requires a list of configuration options and produces an executable.

Running AREPO requires the executable obtained from compiling the code, an initial conditions file,
a parameter file and may require additional dependencies based on the chosen configuration options.
AREPO produces outputs including snapshots, restart files and text file outputs useful for 
diagnosing potential issues.

Reading and analysing AREPO outputs is facilitated by tools which can be found here: 
https://github.com/rowanjsmith/ArepoISM_tools  
These tools are currently designed to work with HDF5 outputs, but can be modified to work for
binary outputs if necessary.

----------------------------------------------------------------------------------------------------
## I. GET ACCESS TO AREPO

To get access to the AREPO repository from bitbucket using SSH, first go to: 
https://bitbucket.org/account/settings/ssh-keys/  
and add your *public* SSH key (usually ~/.ssh/id_rsa.pub). If you can't find the key or the .ssh 
directory does not exist in your COSMA home directory, you might have to generate a key pair using:  

    ssh-keygen

Once you have added your public SSH key to your bitbucket account, go to your chosen host directory
(e.g. your home directory: ~) then clone the AREPO repository:  

    git clone git@bitbucket.org:volkerspringel/arepo.git

And checkout your chosen branch:  

    git checkout <branch_name>

You can find the status of any branch relative to the main (or 'master') branch on bitbucket:
https://bitbucket.org/volkerspringel/arepo/branches/  

Recommendation as of June 2025:
If using AREPO 1, I recommend using nbrucy/EcoGal_MW if using star particles with additional physics.
Otherwise, the main branch is free of issues for most other configurations.
If using AREPO 2, the Arepo2 branch should work for most configuration (including star particles 
with additional physics).

----------------------------------------------------------------------------------------------------
## II. COMPILE AREPO

1.  From your host directory (e.g. ~/arepo), find and edit the Template-Config.sh file to enable the 
    configuration options you require and disable the configuration options you do not require. 
    Below is an example of config flags for AREPO 1 with MHD, SGChem network 5 (Nelson & Langer 1997), 
    TreeCol radiative cooling, SF_ECOGAL star particles and feedback:

        NTYPES=6
        USE_ENTROPY_FOR_COLD_FLOWS
        ENTROPY_MACH_THRESHOLD=1.05
        IMPOSE_PINNING
        VORONOI
        MHD
        MHD_POWELL
        MHD_SEEDFIELD
        RIEMANN_HLLD
        REGULARIZE_MESH_CM_DRIFT
        REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED
        REGULARIZE_MESH_FACE_ANGLE
        TREE_BASED_TIMESTEPS
        ENLARGE_DYNAMIC_RANGE_IN_TIME
        REFINEMENT_SPLIT_CELLS
        REFINEMENT_MERGE_CELLS
        REFINEMENT_VOLUME_LIMIT
        SELFGRAVITY
        GRAVITY_NOT_PERIODIC
        EVALPOTENTIAL
        ADAPTIVE_HYDRO_SOFTENING
        CHUNKING
        DOUBLEPRECISION=1
        DOUBLEPRECISION_FFTW
        VORONOI_DYNAMIC_UPDATE
        NO_MPI_IN_PLACE
        NO_ISEND_IRECV_IN_DOMAIN
        FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
        FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
        OVERRIDE_PEANOGRID_WARNING
        OUTPUTPOTENTIAL
        HAVE_HDF5
        HOST_MEMORY_REPORTING
        MHD_DONT_PRINT_BMAGSUM
        SNE_FEEDBACK
        TREECOLV2
        TREECOLV2_CO
        TREECOLV2_H2
        NSIDE=2
        SGCHEM
        CHEMISTRYNETWORK=5
        SGCHEM_TEMPERATURE_FLOOR
        SF_ECOGAL
        SF_ECOGAL_CHECKS
        SF_ECOGAL_ADDITIONAL_OUTPUT
        SF_ECOGAL_FEEDBACK
        SF_ECOGAL_FEEDBACK_REINJECT_MASS
        SF_ECOGAL_FEEDBACK_PHOTOION

2.  Uncomment the following line in Makefile.systype:  
    
        SYSTYPE="cosma"  
    
    Make sure all other lines are commented out (ie. starting with #)

3.  In makefiles/systypes.make, make sure the compile options for SYSTYPE "cosma" matches the 
    following.

    FOR AREPO 1:  
    
        ifeq ($(SYSTYPE), "cosma")
            CC       =  mpicc
            FC       =  mpif90 -nofor-main
            OPTIMIZE =  -std=c99 -O2 -g -DH5_USE_16_API
            GSL_INCL =
            GSL_LIBS =  -lgsl
            FFTW_INCL=
            FFTW_LIBS=  -lfftw
            HDF5INCL =  -DH5_USE_16_API
            HDF5LIB  =  -lhdf5
            MPICHLIB = -lmpi
            HWLOC_INCL= -I/usr/include
            HWLOC_LIB = $(LDFLAGS) -lhwloc
            LINKER   = $(FC)
        endif

    FOR AREPO 2:  

        ifeq ($(SYSTYPE), "cosma")
            CC       = mpicc
            CPPC     = mpicxx -std=c++11
            FC       = mpif90 -nofor-main
            OPTIMIZE = -g -O2 -DH5_USE_16_API
            GSL_INCL =
            GSL_LIBS =  -lgsl
            FFTW_INCL=
            FFTW_LIBS=  -lfftw
            HDF5INCL =  -DH5_USE_16_API
            HDF5LIB  =  -lhdf5
            MPICHLIB = -lmpi
            HWLOC_INCL= -I/usr/include
            HWLOC_LIB = $(LDFLAGS) -lhwloc
            LINKER   = $(FC)
        endif
    
4.  Purge existing modules then load the required modules on COSMA:  
    
        module purge # purge all existing modules

        module load intel_comp/2024.2.0
        module load compiler-rt tbb compiler mpi
        module load gsl/2.8
        module load fftw/3.3.10cosma8 
        module load hdf5/1.14.4
        module load cosma/2024 # cosma/2018
        module load python/3.12.4
        module load armforge/23.1.0
        module load hdfview/3.3.2 # hdfview/3.1.4 
        module load gadgetviewer/1.1.4
        module load utils/202402
        module load hwloc/2.11.1
        module load allinea/ddt/23.1.0

    (You can copy this into a bash file (e.g. load_modules.sh) and run it using ./load_modules.sh)

5.  Clean existing build files using:  

        make clean

6.  Compile AREPO using:  

        make CONFIG=<your_config_file.sh> EXEC=<your_executable_name>

    This will generate many warnings, but hopefully no errors. If compilation is successful, you 
    will find your executable in the current directory.

----------------------------------------------------------------------------------------------------
## III. RUN AREPO

In your data directory (typically "/cosma8/data/<project_ID>/<user_ID>/"), create a new directory to 
run AREPO from, or copy the example directory found here: 
https://github.com/rowanjsmith/ArepoISM_tools/arepo_run_x

This directory should contain the following:
-   an AREPO executable obtained from compiling the code
-   a parameter file (typically "in.param")
-   a "snapshots" directory containing your initial conditions or initial snapshot
-   a batch job submission file (typically "batchsub")
-   If your executable was compiled with TreeCol enabled (option "TREECOLV2"), you will also need a 
    file called "TreeCol_lookup.dat". You can copy this file from the example directory.

NOTE: the example directory does not contain an executable. Copy or move the executable obtained 
      from compiling the code to your run directory.

Go through the parameter file to specify your chosen parameters. Read the example parameter file for 
information about each parameter.

Go through the batchsub file to specify your batch job submission options. Read the example batchsub
file for information about each option.

Once the parameters and batchsub options have been set, you can submit your batch job using:  

    sbatch batchsub

To find your job in the queue and check its status, use:  

    squeue --me

This will show you the following information:
-   JOBID is your job identifier. This can be used to cancel your job if submitted by mistake 
        using: "scancel <JOBID>". Make sure you have the correct job ID before using this command.
-   PARTITION is where the job is/will be running. If you are on COSMA8, this should be "cosma8".
-   NAME is the job name you specified in the batchsub file using "-J".
-   USER is the user who submitted the job. This should be your user name.
-   ST is the status of the job. "PD" stands for pending, indicating that you job is in the 
        queue, while "R" stands for running. 
-   TIME is the amount of time your job has been running for. It will start when you job status
        changes from "PD" to "R". If present, the prefix indicates the number of days (e.g. 
        2-01:00:00 means your job has been running for 2 days and 1 hour).
-   NODES is the number of nodes allocated to your job. One node has 128 processors (ntasks).
-   NODELIST (REASON) lists either the node ID(s) allocated to your job if it is running, or
        it indicates the reason your job is pending (e.g. "(Resources)" or "(Priority)")
        
If the "squeue --me" command does not show anything, it either means your job has terminated 
(successfully or unsuccessfully) or it wasn't submitted correctly. 

If your job has run, you will find the output file (e.g. "arepo_out.txt") in your run directory. 
The start of this file lists all of the configuration options your executable was compiled with as 
well as all of the input parameters used when the code was run. If your job terminates early or is 
unsuccessful, the end of this file will give you information about the cause of the termination, 
usually with an error message containing a traceback. If your job was unsuccessful and the end of 
your output file does not contain an error message, the last line will give you a starting point to 
find where the issue might have occured. Keep in mind that your job will terminate when it reaches 
the "TimeLimitCPU" parameter value - this will cause the output to cease abruptly which can be 
mistaken for a failure, as opposed to the program working as intended. The easiest way to find
whether that is the case is to check the email notification for the job termination and compare 
the run time to the "TimeLimitCPU" value.

Your first job will most likely fail within the first few minutes due to a missing a parameter or 
and invalid parameter in your parameter file. If you have a missing parameter, the output text file
will specify its name so you can add it to the parameter file, along with a reasonable value. If you 
have an invalid parameter, either remove that line from the parameter file or comment it out.

Note: In some cases, determining a reasonable value can be difficult if you are unsure what the 
parameter does or how it is used. The first steps would be to check the example parameter file for 
any helpful comments, and the documentation for the module to which the parameter belongs (see the 
documentation directory in the arepo codebase). If that is unhelpful, further sleuthing will be 
necessary. If the parameter is for an optional module, the paper introducing that module might 
contain useful information about required parameters. If the parameter is for the base code, the 
original AREPO paper (Springel 2010 "E pur si muove") or the public release paper (Weinberger, 
Springel and Pakmor 2020) might contain useful information. If a reasonable value cannot be 
determined from the literature, search for the parameter name in the codebase to understand how it 
is used.

----------------------------------------------------------------------------------------------------
## IV. READ AND ANALYSE AREPO OUTPUTS

Coming soon to a ~~movie theatre~~ repository near you...

In the meantime, have a look at snapshot_tools.py, plotter.py and plotter_example.py here:
https://github.com/rowanjsmith/ArepoISM_tools

----------------------------------------------------------------------------------------------------

Written by Zoe Faes with wisdom from Rowan Smith, Sansith Hewapathirana, David Whitworth, 
Kammy Bogue, Philipp Girichidis, Robin Tress, Junia Goeller and Ruediger Pakmor.
Last edited: 10 Jun 2025