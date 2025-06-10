#!/bin/bash            # this line only there to enable syntax highlighting in this file

#--------------------------------------- Basic operation mode of the code
NTYPES=6                       # number of particle types
ISOTHERM_EQS
USE_ENTROPY_FOR_COLD_FLOWS
ENTROPY_MACH_THRESHOLD=1.05

#--------------------------------------- MPI/Threading Hybrid
IMPOSE_PINNING

#--------------------------------------- Mesh Type
VORONOI

#--------------------------------------- Mesh motion and regularization
REGULARIZE_MESH_CM_DRIFT
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED
REGULARIZE_MESH_FACE_ANGLE

#--------------------------------------- Time integration options
TREE_BASED_TIMESTEPS      # non-local time step criterion (take “signal speed” into account)
ENLARGE_DYNAMIC_RANGE_IN_TIME  # This extends the dynamic range of the integer timeline from 32 to 64 bit

#--------------------------------------- Refinement and derefinement
REFINEMENT_SPLIT_CELLS
REFINEMENT_MERGE_CELLS
REFINEMENT_VOLUME_LIMIT

#--------------------------------------- Mesh-relaxing or mesh-adding (this will not carry out a simulation)
ADDBACKGROUNDGRID

#--------------------------------------- Gravity treatment
SELFGRAVITY                    # switch on for self-gravity
GRAVITY_NOT_PERIODIC          # if gravity is not to be treated periodically

#--------------------------------------- Gravity softening
ADAPTIVE_HYDRO_SOFTENING

#--------------------------------------- Things that are always recommended
CHUNKING                 # will calculate the gravity force in interleaved blocks.
                         # This can reduce imbalances in case multiple iterations due to
                         # insufficient buffer size need to be done.

#--------------------------------------- Single/Double Precision
DOUBLEPRECISION=1
DOUBLEPRECISION_FFTW

#-------------------------------------------- Things for special behavior
VORONOI_DYNAMIC_UPDATE          # keeps track of mesh connectivity, which speeds up mesh construction
SHIFT_BY_HALF_BOX
NO_MPI_IN_PLACE
NO_ISEND_IRECV_IN_DOMAIN
FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES  # this can be used to load SPH ICs that contain identical particle coordinates
OVERRIDE_PEANOGRID_WARNING

#--------------------------------------- Output/Input options
HAVE_HDF5                     # needed when HDF5 I/O support is desired

#--------------------------------------- Testing and Debugging options
HOST_MEMORY_REPORTING          # reports the available system memory after start-up by analyzing /proc/meminfo

#--------------------------------------- Miscellaneous
UPDATE_GRADIENTS_FOR_OUTPUT # New in Arepo2
