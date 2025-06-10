----------------------------------------------------------------------------------------------------
Last updated: 10 Jun 2025
----------------------------------------------------------------------------------------------------

This is a guide to create initial conditions (ICs) for isolated galaxy simulations in AREPO with 
pNbody (Revaz 2013) on the Distributed Research utilising Advanced Computing (DiRAC) facility's 
COSmological MAchine 8 (COSMA8).

pNbody documentation can be found here: https://www.astro.unige.ch/~revaz/pNbody/

This guide includes installation instructions for pNbody on COSMA8 and a tutorial to create initial 
conditions. Example files mentioned in this file can be found here:
https://github.com/rowanjsmith/ArepoISM_tools/pNbody_ics/

Note: any string contained in arrows (<some_string>) indicates a placeholder. For example, the
specific path "cosma/home/dp000/dc-aaaa1/" can be represented as "cosma/home/<project_ID>/<user_ID>/"
where <project_ID> is "dp000" and <user_ID> is "dc-aaaa1".

----------------------------------------------------------------------------------------------------

I. Installation instructions for pNbody on COSMA8

Install pNbody in a dedicated virtual environment using "venv". 
(If using VS Code, press F1 then select: "Python: Create Environment ...", then select: "venv")

Create a new virtual environment using: 
    $ python3 -m venv <pNbody_env>

Activate it using: 
    $ source <pNbody_env>/bin/activate

Your prompt should now be preceded by (<pNbody_env>)
 
Clone the pNbody repository using: 
    $ git clone https://gitlab.com/revaz/pNbody.git 
 
Check that you have python/3.12.4 using: 
    $ python --version

Otherwise, upgrade to python/3.12.4 using: 
    $ pip install python=3.12.4

You must use the correct cblas library: 
    $ export LD_PRELOAD=/usr/lib64/libcblas.so
 
Now, install pNbody with: 
    $ pip install ./pNbody
 
This should also install all dependencies required by pNbody.

Check installation using: 
    $ pNbody_test

Output should finish with:
    ########################################################################
    Good News ! pNbody with format swift, gadget and gh5 is working !
    ########################################################################

----------------------------------------------------------------------------------------------------

II. Tutorial to generate initial conditions with pNbody on COSMA8

See: https://www.astro.unige.ch/~revaz/pNbody/rst/GeneratingAMulitComponentsGalaxy.html

To generate ICs, pNbody requires a parameter file (see: example_params.yml).

Follow the steps below to create a valid initial conditions file.

1.  Activate your pNbody virtual environment: 
    $ source <pNbody_env>/bin/activate

2.  Configure the cblas library: 
    $ export LD_PRELOAD=/usr/lib64/libcblas.so

3.  Generate ICs with pNbody: 
    $ ic_makeGalaxy -p <params>.yml -o <ic_name>.hdf5

pNbody was not designed for AREPO, so the output needs to be modified to work with AREPO. 

4.  From ic_tools.py, use fix_ics to correct formatting errors:
    fix_ics(filepath=/<path_to>/<ic_name>.hdf5, 
            newfilename=None, 
            save_path=None, 
            boxsize=<your_boxsize>, 
            unit_length=<your_unit_length>, 
            unit_mass=<your_unit_mass>, 
            unit_velocity=<your_unit_velocity>)

    This will create a new IC file (by default: new_<ic_name>.hdf5 in the same directory as your 
    original IC file). 
    
This new IC file needs to be pre-processed with AREPO to add a background grid
to aid the mesh construction method and to translate from centered coordinates in a box with domain
[-boxsize/2, boxsize/2] to a box with domain [0, boxsize]

5.  To pre-process your ICs, compile AREPO with the following options enabled in the config file:
        ADDBACKGROUNDGRID
        SHIFT_BY_HALF_BOX

6.  Then, modify your parameter file to add the parameter "GridSize" set to 16. 

7.  Run AREPO with this parameter file and your new_<ic_name>.hdf5 initial conditions file. 
    If this is successful, you will find a new snapshot in your "snapshots" directory with the suffix 
    "-with-grid" (e.g. new_<ic_name>-with-grid.hdf5).

This is the initial conditions file you can use to start a new run. 

At this stage, it does not contain data for chemistry or magnetic fields. If required, these can be 
added using add_chem_abundances, add_dust_temperature and add_magnetic_field from ic_tools.py

--------------------
If using the example files, the pre-processing steps 5, 6 and 7 can be carried out as follows:

Compile the executable from your AREPO directory using the pre-processing configuration options:
    $ make CONFIG=config_preproc.sh EXEC=Arepo2_preproc

Use the batch job submission file "batchsub_preproc" linked to the parameter file "preproc.param" to 
add the background grid and translate the coordinates:
    $ sbatch batchsub_preproc
--------------------

Once the final "new_<ic_name>-with-grid.hdf5" IC file is obtained, you can check it using check_ics
from ic_tools.py:
    check_ics(filepath=/<path_to>/new_<ic_name>-with-grid.hdf5, 
              plots=True, 
              boxsize=<your_boxsize>, 
              unit_length=<your_unit_length>, 
              unit_mass=<your_unit_mass>, 
              unit_velocity=<your_unit_velocity>)

This function will check whether the boxsize and units correspond to the given values, make basic 2D
scatter plots of the gas particle positions and give you the option to remove gas particles located 
outside of the boxsize (this can sometimes happen when adding the background grid).

Compile AREPO with your chosen options and run with the required parameters and your checked
initial conditions file.
----------------------------------------------------------------------------------------------------

Written by Zoe Faes with wisdom from Robin Tress, David Whitworth, and the pNbody documentation
