# -*- coding: utf-8 -*-
"""
Created on Fri May 23 19:28:01 2025

@author: zoefaes
"""

import h5py
import os
import warnings
import numpy as np
import astropy.units as u
from astropy import constants
import matplotlib.pyplot as plt


def fix_ics(filepath, newfilename=None, save_path=None, boxsize=1000, unit_length=3.0856e20, unit_mass=1.991e33, unit_velocity=1.0e5):
    """
    Fixes pNbody initial conditions by removing the first gas particle with id = 0 and fixing the mass table.
    """

    if newfilename is None:
        filename = os.path.basename(filepath)
        newfilename = f'new_{filename}'

    if save_path is None:
        path = os.path.dirname(filepath)
    else:
        path = save_path

    # create new file
    newfilepath = path + '/' + newfilename
    newsnap = h5py.File(newfilepath, 'w')
    Header = newsnap.create_group('Header')
    part0 = newsnap.create_group('PartType0')
    part1 = newsnap.create_group('PartType1')
    part2 = newsnap.create_group('PartType2')

    # read ics
    f = h5py.File(filepath, 'r')
    print(newsnap)

    # Fix NumPart_ThisFile
    NumPart_ThisFile = f['Header'].attrs['NumPart_ThisFile']
    NumPart_ThisFile[2] = NumPart_ThisFile[4]
    NumPart_ThisFile[4] = 0
    NumPart_ThisFile[0] = NumPart_ThisFile[0] - 1 # skip first gas particle with id = 0
        
    for item in f['Header'].attrs:
        # Fix MassTable
        if item == 'MassTable':
            MassTable = f['Header'].attrs['MassTable']
            MassTable[2] = 0.0
            Header.attrs.create('MassTable', MassTable)
        # Fix NumPart_ThisFile
        elif item == 'NumPart_ThisFile':
            Header.attrs.create('NumPart_ThisFile', NumPart_ThisFile)
        # Fix NumPart_Total
        elif item == 'NumPart_Total':
            Header.attrs.create('NumPart_Total', NumPart_ThisFile)
        # Fix units
        elif item == 'UnitLength_in_cm':
            Header.attrs.create(item, unit_length)
        elif item == 'UnitMass_in_g':
            Header.attrs.create(item, unit_mass)
        elif item == 'UnitVelocity_in_cm_per_s':
            Header.attrs.create(item, unit_velocity)
        # Fix BoxSize
        elif item == 'BoxSize':
            Header.attrs.create(item, boxsize)
        else:
            Header.attrs.create(item, f['Header'].attrs[item])
                
        
    for item in f['PartType0'].keys():
        part0.create_dataset(item, data=f['PartType0'][item][1:]) # skip the first gas particle with id = 0
            
    for item in f['PartType1'].keys():
        part1.create_dataset(item, data=f['PartType1'][item][:])
                    
    # Fix PartType mixup
    for item in f['PartType4'].keys():
        part2.create_dataset(item, data=f['PartType4'][item][:])
            
    f.close()
    newsnap.close()

    # Check
    snap = h5py.File(newfilepath, 'r')

    print('\nCheck new file header:')
    for item in snap['Header'].attrs:
        print(item, snap['Header'].attrs[item])

    snap.close()

    print('\nFixed initial conditions and saved to ', newfilepath)


def check_ics(filepath, plots=True, boxsize=1000, unit_length=3.0856e20, unit_mass=1.991e33, unit_velocity=1.0e5):
    """
    Check if the initial conditions are correct.
    """
    # load file
    f = h5py.File(filepath, 'r+')

    # Check header
    print("\nHEADER:")
    for item in f['Header'].attrs:
        print(item, f['Header'].attrs[item])

    try:
        ic_boxsize = f['Header'].attrs['BoxSize']
        if ic_boxsize != boxsize:
            print(f"Boxsize is {ic_boxsize}, not {boxsize}. Check the file!")
    except KeyError:
        raise Exception("BoxSize not found in Header")

    try:
        ic_unit_length = f['Header'].attrs['UnitLength_in_cm']
        ic_unit_mass = f['Header'].attrs['UnitMass_in_g']
        ic_unit_velocity = f['Header'].attrs['UnitVelocity_in_cm_per_s']
    except KeyError:
        ic_unit_length = f['Parameters'].attrs['UnitLength_in_cm']
        ic_unit_mass = f['Parameters'].attrs['UnitMass_in_g']
        ic_unit_velocity = f['Parameters'].attrs['UnitVelocity_in_cm_per_s']

    if ic_unit_length != unit_length:
        print(f"Unit length is {ic_unit_length}, not {unit_length}. Check the file!")
    if ic_unit_mass != unit_mass:
        print(f"Unit mass is {ic_unit_mass}, not {unit_mass}. Check the file!")
    if ic_unit_velocity != unit_velocity:
        print(f"Unit velocity is {ic_unit_velocity}, not {unit_velocity}. Check the file!")

    # Check gas particles positions
    print("\nGAS PARTICLES:")
    print('Number of gas particles: ', len(f['PartType0']['ParticleIDs']))

    x, y, z = np.array(f['PartType0']['Coordinates']).T
    print("Gas particle coordinate ranges:")
    print("x: ", np.min(x), np.max(x))
    print("y: ", np.min(y), np.max(y))
    print("z: ", np.min(z), np.max(z))

    if plots == True:
        # quick plot gas particles
        plt.scatter(x, y)
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('Gas particle positions')
        plt.show()

        plt.scatter(x, z)
        plt.xlabel('X')
        plt.ylabel('Z')
        plt.title('Gas particle positions')
        plt.show()

    # Check out of bounds gas particles
    out_of_bounds = np.argwhere((np.array(f['PartType0']['Coordinates']) < 0).any(axis=1) | (np.array(f['PartType0']['Coordinates']) > boxsize).any(axis=1))

    print(f'\nNumber of gas particles outside boxsize: {len(out_of_bounds)}')

    if len(out_of_bounds) > 0:
        while True:
            answer = input(f"Do you want to remove {len(out_of_bounds)} gas particles out of bounds? (y/n): ").strip().lower()
            if answer in ['y', 'yes']:

                print("Removing out of bounds gas particles...")
                # Remove out of bounds gas particles
                group = f['PartType0']
                for dataset_name in group:
                    data = group[dataset_name][:]
                    
                    # Delete the specified index (axis=0 → row)
                    new_data = np.delete(data, out_of_bounds, axis=0)
                    
                    # Overwrite the dataset
                    del group[dataset_name]  # delete old dataset
                    group.create_dataset(dataset_name, data=new_data)
                    
                # Update NumParts
                f['Header'].attrs['NumPart_ThisFile'][0] = len(f['PartType0']['ParticleIDs'])
                f['Header'].attrs['NumPart_Total'][0] = len(f['PartType0']['ParticleIDs'])

                print(f"Removed {len(out_of_bounds)} gas particles out of bounds.")
                # Recalculate the number of gas particles
                out_of_bounds = np.argwhere((np.array(f['PartType0']['Coordinates']) < 0).any(axis=1) | (np.array(f['PartType0']['Coordinates']) > boxsize).any(axis=1))
                print(f'Number of gas particles outside boxsize: {len(out_of_bounds)}')

                if plots == True:
                    x, y, z = np.array(f['PartType0']['Coordinates']).T

                    # quick plot gas particles
                    plt.scatter(x, y)
                    plt.xlabel('X')
                    plt.ylabel('Y')
                    plt.title('Gas particles')
                    plt.show()

                    plt.scatter(x, z)
                    plt.xlabel('X')
                    plt.ylabel('Z')
                    plt.title('Gas particles')
                    plt.show()

                break
            elif answer in ['n', 'no']:
                print("Exiting...")
                break
            else:
                print("Please enter 'y' or 'n'.")

    print(f"\nMean dark matter particle mass: {np.mean(f['PartType1']['Masses'])}")
    print(f"BH particle mass: {np.max(f['PartType1']['Masses'])}")

    f.close()


def add_snapshot_property(filepath: str,
                          snap_group: str,
                          property: str,
                          data):
    # write to hdf5 file
    with h5py.File(filepath, 'a') as f:
        if snap_group in f.keys():

            if snap_group[0:8] == 'PartType':
                if len(data) != len(f.get(snap_group)['ParticleIDs']):
                    warnings.warn(f"Length of {property} data does not match number of particles for {snap_group}.", UserWarning)
                    
            f[snap_group].create_dataset(property, data=data)
            print(f"{property} successfully added to {snap_group}.")
        else:
            raise Exception(f"{snap_group} not found in the HDF5 file.")
        

def update_snapshot_property(filepath: str,
                             snap_group: str,
                             property: str,
                             data):
    # write to hdf5 file
    with h5py.File(filepath, 'r+') as f:
        if snap_group in f.keys():
            if property in f[snap_group]:
                if not len(f[snap_group][property]) == len(data):
                    if snap_group == ('Config' or 'Header' or 'Parameters'):
                        warnings.warn(f"Existing {property} data and new {property} data have different lengths.", UserWarning)
                    else:
                        raise Exception(f"Existing {property} data and new {property} data must have the same length (corresponding to the number of cells in snapshot).")
                # check for same shape
                elif not np.shape(f[snap_group][property]) == np.shape(data):
                    warnings.warn(f"Existing {property} data and new {property} data have different shapes. {property} cannot be modified in place and will be overwritten.", UserWarning)
                    del f[snap_group][property]
                    f[snap_group].create_dataset(property, data=data)

                    print(f"{property} successfully updated.")
                else:
                    # Overwrite existing dataset with new data
                    f[snap_group][property][:] = data  # Modify in place

                    # f.flush()

                    print(f"{property} successfully updated.")
            else:
                raise Exception(f"{property} not found in {snap_group}.")
        else:
            raise Exception(f"{snap_group} not found in the HDF5 file.")
        

def add_bh_particle(filepath, bh_mass_in_Msun=None):
    """Add mass to a dark matter particle at the 
    center of the box to account for SMBH potential.
    If unspecified, the SMBH mass is set in relation to 
    the stellar mass according to Reines and Volonteri 2015 
    (https://arxiv.org/pdf/1508.06274)

    Parameters
    ----------
    filepath : _type_
        _description_
    bh_mass : _type_, optional
        _description_, by default None
    """
    f = h5py.File(filepath, 'r+')

    if bh_mass_in_Msun==None:
        # set black hole mass according to Reines and Volonteri 2015 (https://arxiv.org/pdf/1508.06274)
        stellar_mass = sum(f['PartType2']['Masses'])*(f['Header'].attrs['UnitMass_in_g']/1.991E33)
        bh_mass_in_Msun = 10**(7.45 + 1.05 * np.log10(stellar_mass / 1e11))
        print(f"Black hole mass set to {bh_mass_in_Msun:.2e} M_sun")

    # Find dm particle closest to the center
    hx, hy, hz = np.array(f['PartType1']['Coordinates']).T
    boxsize = f['Header'].attrs['BoxSize']
    center = boxsize / 2
    bh = np.argmin(np.abs(hx - center) + np.abs(hy - center) + np.abs(hz - center))

    # Re-center the dm particle
    f['PartType1']['Coordinates'][bh] = [center, center, center]
    # set the velocity of the black hole to zero
    f['PartType1']['Velocities'][bh] = [0.0, 0.0, 0.0]
    # Add mass to the dm particle to account for the mass of the black hole
    f['PartType1']['Masses'][bh] += bh_mass_in_Msun/(f['Header'].attrs['UnitMass_in_g']/1.991E33)
    
    f.close()

    print('Hacky black hole fix completed!')



def add_dust_temperature(filepath, dust_temperature_in_K=2.0):
    """
    Add dust temperature to the initial conditions file.
    """
    # load file
    f = h5py.File(filepath, 'r+')

    # set the same dust temperature for all particles
    dust_temperature = np.array([dust_temperature_in_K] * len(f['PartType0']['ParticleIDs']))
    add_snapshot_property(filepath, 'PartType0', 'DustTemperature', dust_temperature)

    f.close()


def add_chem_abundances(filepath, chem_abundances='NL97'):
    """
    Add chemical abundances to the initial conditions file.
    """
    # load file
    f = h5py.File(filepath, 'r+')

    # chemical abundances are specified directly
    if isinstance(chem_abundances, list):
        chem_abundances = np.array([chem_abundances] * len(f['PartType0']['ParticleIDs']))
        add_snapshot_property(filepath, 'PartType0', 'ChemicalAbundances', chem_abundances)

    # use standard chemical abundances for a given network
    else:
        if isinstance(chem_abundances, str):
            chem_abundances = chem_abundances.lower().strip()
        elif isinstance(chem_abundances, int):
            if chem_abundances == 5:
                chem_abundances = 'nl97'
            elif chem_abundances == 16:
                chem_abundances = 'gong'
            else:
                raise ValueError("Unknown chemical network")

        # NL97 (network 5)
        if chem_abundances in ['nl97', 'nl1997', 'nelson langer 1997', 'nelson langer 97', 'nelson and langer 1997', 'nelson and langer 97', 'network 5', 'network5', 'sgchem network 5']:
            chem_abundances = np.array([[0, 1.0e-02, 0]] * len(f['PartType0']['ParticleIDs']))
            add_snapshot_property(filepath, 'PartType0', 'ChemicalAbundances', chem_abundances)

        # Gong (network 16)
        elif chem_abundances in ['gong', 'gong18', 'gong2018', 'gong et al 2018', 'gong et al 18', 'gong et al. 2018', 'gong et al. 18', 'network 16', 'network16', 'sgchem network 16']:
            chem_abundances = np.array([[0, 1.0e-02, 1.5e-05, 0, 0, 0, 0, 0, 0]] * len(f['PartType0']['ParticleIDs']))
            add_snapshot_property(filepath, 'PartType0', 'ChemicalAbundances', chem_abundances)

        else:
            raise ValueError("Unknown chemical network: Available networks are NL97 (network 5) and Gong (network 16). Or, specify the abundances directly as a list.")

    f.close()


def add_magnetic_field(filepath, B0_in_gauss=1e-6, direction='toroidal'):
    """
    Add magnetic field to the initial conditions file.
    """
    # load file
    f = h5py.File(filepath, 'r+')

    boxsize = f['Header'].attrs['BoxSize']

    # generate magnetic field
    mag_field = generate_magnetic_field(f['PartType0']['Coordinates'], boxsize=boxsize, B0=B0_in_gauss, direction=direction)
    add_snapshot_property(filepath, 'PartType0', 'MagneticField', mag_field)

    f.close()


def generate_magnetic_field(points: np.ndarray,
                            boxsize: float,
                            B0: float, 
                            direction: str = 'toroidal') -> np.ndarray:
    """
    Generates standard initial seed magnetic field configurations.

    Parameters
    ----------
    points : np.ndarray
        Position coordinates for which the magnetic field is to be defined on.
    boxsize : float
        Length of box
    B0 : float
        Magnitude of magnetic field
    direction : str, optional
        Direction of the magnetic field. Choose from 'toroidal', 'poloidal', 'x', 'y', or 'z'. 
        Direction can be reversed by setting magnitude 'B0' to a negative value. Default is 'toroidal'

    Returns
    -------
    np.ndarray
        Magnetic field defined on 'points'.
    """    

    # input validation
    if not isinstance(B0, (int, float)):
        raise Exception(f"'B0' should be a float or an integer, not {type(B0)}")
    if not isinstance(direction, str):
        raise Exception(f"'direction' should be a string, not {type(direction)}")
    else:
        direction = direction.strip().lower()

    # generate magnetic seedfield
    x, y, z = np.array(points).T
    c = boxsize/2
    if direction == 'toroidal':
        r = np.sqrt((x-c)**2 + (y-c)**2)
        r[r == 0] = 1e-10           # Prevent singularity at r = 0
        bx = -B0 * (y-c)/r       # sin(arctan(y/x)) = y/sqrt((x**2 + y**2)
        by = B0 * (x-c)/r        # cos(arctan(y/x)) = x/sqrt(x**2 + y**2)
        Bx = np.full_like(x, bx)
        By = np.full_like(y, by)
        Bz = np.zeros_like(z)
    # TODO: make 'poloidal' generate an actual poloidal field
    elif (direction == 'poloidal') or (direction == 'z'): # initial seed field for poloidal is same as z-direction?
        Bx = np.zeros_like(x)
        By = np.zeros_like(y)
        Bz = np.full_like(z, B0)
    elif direction == 'x':
        Bx = np.full_like(x, B0)
        By = np.zeros_like(y)
        Bz = np.zeros_like(z)
    elif direction == 'y':
        Bx = np.zeros_like(x)
        By = np.full_like(y, B0)
        Bz = np.zeros_like(z)
    else:
        raise Exception("Invalid string given for 'direction'. Please choose from 'toroidal', 'poloidal', 'x', 'y', or 'z'.")
    
    B_field = np.array([Bx, By, Bz]).T
    
    return B_field

def disk_mask(coordinates, boxsize, unit_length_in_kpc, radii_in_kpc = [0, 15], half_heights_in_kpc = [0, 0.5]):

    xc, yc, zc = (np.array(coordinates).T - boxsize/2) * unit_length_in_kpc
    r = np.sqrt(xc*xc + yc*yc)
    theta = np.arctan2(yc, xc)

    disk = ((r > radii_in_kpc[0]) & (r < radii_in_kpc[1])) & ((np.abs(zc) < half_heights_in_kpc[1]) & (np.abs(zc) > half_heights_in_kpc[0]))
    print(f"Disk mask: {np.sum(disk)} particles in disk with radii {radii_in_kpc} kpc and half heights {half_heights_in_kpc} kpc")

    return disk

def cgm_mask(coordinates, boxsize, unit_length_in_kpc, radii_in_kpc = [20, 100], half_heights_in_kpc = [2, 100]):

    xc, yc, zc = (np.array(coordinates).T - boxsize/2) * unit_length_in_kpc
    r = np.sqrt(xc*xc + yc*yc)
    theta = np.arctan2(yc, xc)

    cgm = ((r > radii_in_kpc[0]) & (r < radii_in_kpc[1])) | ((np.abs(zc) < half_heights_in_kpc[1]) & (np.abs(zc) > half_heights_in_kpc[0]))
    print(f"CGM mask: {np.sum(cgm)} particles in region with radii {radii_in_kpc} kpc and half heights {half_heights_in_kpc} kpc")

    return cgm

def dens_mask(densities, unit_density, density_threshold_cgs = 1e-37):

    mask = ((densities * unit_density) < density_threshold_cgs)
    print(f"Density mask: {np.sum(mask)} particles in region with density below {density_threshold_cgs} g/cm^3")

    return mask

def temp_mask(temperatures, temperature_threshold = 1e10):

    mask = (temperatures > temperature_threshold)
    print(f"Temperature mask: {np.sum(mask)} particles in region with temperature above {temperature_threshold} K")

    return mask

def get_number_density():
    """
    Calculates number density using 
    :math:'n = \\frac{\\rho}{(1+4x_\\mathrm{He})m_\\mathrm{p}}'.

    Parameters
    ----------

    density: '~astropy.units.Quantity'
        Density :math:'\\rho'.

    He_abundance: float, optional
        Helium abundance (by number). Assumed value is '0.1'.

    Returns
    -------

    ndensity : '~astropy.units.Quantity'
        Number density :math:'n'.

    """
    density = self.density.to(u.g/u.cm**3)
    xHe = 0.1 if (self.xHe is None) else self.xHe
    ndensity = density/((1.0 + 4.0 * xHe) * constants.m_p.cgs)

    return ndensity

def get_disk_mass_ndensity_temp(filepath, radii_in_kpc = [0, 15], half_heights_in_kpc = [0, 0.5], xHe=0.1):
    """
    Get the number density per cubic cm and temperature in Kelvin of the disk.
    """
    # load file
    f = h5py.File(filepath, 'r+')

    # get units
    unit_mass = f['Parameters'].attrs['UnitMass_in_g']
    unit_mass_in_solmass = unit_mass / 1.9884e33
    unit_length = f['Parameters'].attrs['UnitLength_in_cm']
    unit_length_in_kpc = unit_length / 3.0856e21
    unit_velocity = f['Parameters'].attrs['UnitVelocity_in_cm_per_s']
    unit_density = unit_mass / unit_length**3

    # get gas properties
    coords = np.array(f['PartType0']['Coordinates'])
    mass = np.array(f['PartType0']['Masses'])
    density = np.array(f['PartType0']['Density'])
    internal_energy = np.array(f['PartType0']['InternalEnergy'])
    xH2, xHp, xCO = np.array(f['PartType0']['ChemicalAbundances']).T

    boxsize = f['Header'].attrs['BoxSize']

    f.close()

    radii_in_kpc = np.array(radii_in_kpc)
    half_heights_in_kpc = np.array(half_heights_in_kpc)

    # get disk
    disk = disk_mask(coords, boxsize, unit_length_in_kpc, radii_in_kpc, half_heights_in_kpc)

    # get mass of disk
    disk_mass = np.sum(mass[disk]) * unit_mass_in_solmass
    print('Total disk mass (M_sun): ', disk_mass)

    # get density of disk
    disk_ndensity = density[disk] * unit_density / ((1.0 + 4.0 * xHe) * constants.m_p.cgs.value)
    print('mean disk_density (g cm-3): ', np.mean(density[disk] * unit_density))
    print('mean disk_ndensity (cm-3): ', np.mean(disk_ndensity))

    # get temperature of disk
    xTOT = 1.0 + xHp[disk] - xH2[disk] + xHe
    nTOT = xTOT * disk_ndensity
    mean_molecular_weight = density[disk] * unit_density / nTOT
    # print(f"Mean Molecular Weight: \nMean: {np.mean((mean_molecular_weight/constants.m_p).decompose())}, \nMax: {np.max((mean_molecular_weight/constants.m_p).decompose())}, \nMin: {np.min((mean_molecular_weight/constants.m_p).decompose())}")
    disk_temp = (2.0/3.0) * internal_energy[disk] * unit_velocity**2 * mean_molecular_weight / constants.k_B.cgs.value
    print('mean disk_temp (K): ', np.mean(disk_temp))

    return disk_mass, disk_ndensity, disk_temp

def get_cgm_mass_ndensity_temp(filepath, radii_in_kpc = [20, 100], half_heights_in_kpc = [1, 50], xHe=0.1):
    """
    Get the number density per cubic cm and temperature in Kelvin of the disk.
    """
    # load file
    f = h5py.File(filepath, 'r+')

    # get units
    unit_mass = f['Parameters'].attrs['UnitMass_in_g']
    unit_mass_in_solmass = unit_mass / 1.9884e33
    unit_length = f['Parameters'].attrs['UnitLength_in_cm']
    unit_length_in_kpc = unit_length / 3.0856e21
    unit_velocity = f['Parameters'].attrs['UnitVelocity_in_cm_per_s']
    unit_density = unit_mass / unit_length**3

    # get gas properties
    coords = np.array(f['PartType0']['Coordinates'])
    mass = np.array(f['PartType0']['Masses'])
    density = np.array(f['PartType0']['Density'])
    internal_energy = np.array(f['PartType0']['InternalEnergy'])
    xH2, xHp, xCO = np.array(f['PartType0']['ChemicalAbundances']).T

    boxsize = f['Header'].attrs['BoxSize']

    f.close()

    radii_in_kpc = np.array(radii_in_kpc)
    half_heights_in_kpc = np.array(half_heights_in_kpc)

    # get cgm
    cgm = cgm_mask(coords, boxsize, unit_length_in_kpc, radii_in_kpc, half_heights_in_kpc)

    # get mass of cgm
    cgm_mass = np.sum(mass[cgm]) * unit_mass_in_solmass
    print('Total cgm mass (M_sun): ', cgm_mass)

    # get density of cgm
    cgm_ndensity = density[cgm] * unit_density / ((1.0 + 4.0 * xHe) * constants.m_p.cgs.value)
    print('mean cgm_density (g cm-3): ', np.mean(density[cgm] * unit_density))
    print('mean cgm_ndensity (cm-3): ', np.mean(cgm_ndensity))

    # get temperature of cgm
    xTOT = 1.0 + xHp[cgm] - xH2[cgm] + xHe
    nTOT = xTOT * cgm_ndensity
    mean_molecular_weight = density[cgm] * unit_density / nTOT
    # print(f"Mean Molecular Weight: \nMean: {np.mean((mean_molecular_weight/constants.m_p).decompose())}, \nMax: {np.max((mean_molecular_weight/constants.m_p).decompose())}, \nMin: {np.min((mean_molecular_weight/constants.m_p).decompose())}")
    cgm_temp = (2.0/3.0) * internal_energy[cgm] * unit_velocity**2 * mean_molecular_weight / constants.k_B.cgs.value
    print('mean cgm_temp (K): ', np.mean(cgm_temp))

    return cgm_mass, cgm_ndensity, cgm_temp

def set_CGM_ndensity(filepath, CGM_ndensity=1e-4, density_threshold_cgs= 1e-37, radii_in_kpc = [20, 100], half_heights_in_kpc = [2, 100], xHe=0.1):

    # load file
    f = h5py.File(filepath, 'r+')

    # get units
    unit_mass = f['Parameters'].attrs['UnitMass_in_g']
    unit_length = f['Parameters'].attrs['UnitLength_in_cm']
    unit_length_in_kpc = unit_length / 3.0856e21
    unit_density = unit_mass / unit_length**3

    # get gas properties
    coords = np.array(f['PartType0']['Coordinates'])
    density = np.array(f['PartType0']['Density'])

    boxsize = f['Header'].attrs['BoxSize']

    f.close()

    radii_in_kpc = np.array(radii_in_kpc)
    half_heights_in_kpc = np.array(half_heights_in_kpc)

    # get region outside disk
    CGM = dens_mask(density, unit_density, density_threshold_cgs) & cgm_mask(coords, boxsize, unit_length_in_kpc, radii_in_kpc, half_heights_in_kpc) 
    
    print(f"CGM mask: {np.sum(CGM)} particles outside galactic disk under density threshold {density_threshold_cgs} g/cm^3")

    # update densities
    density[CGM] = CGM_ndensity * ((1.0 + 4.0 * xHe) * constants.m_p.cgs.value) / unit_density

    # set densities
    update_snapshot_property(filepath, 'PartType0', 'Density', density)


def set_CGM_temperature(filepath, CGM_temp_K = 1e6, temp_threshold_K = 1e10, radii_in_kpc = [20, 100], half_heights_in_kpc = [2, 100], xHe=0.1):

    # load file
    f = h5py.File(filepath, 'r+')

    # get units
    unit_mass = f['Parameters'].attrs['UnitMass_in_g']
    unit_length = f['Parameters'].attrs['UnitLength_in_cm']
    unit_length_in_kpc = unit_length / 3.0856e21
    unit_velocity = f['Parameters'].attrs['UnitVelocity_in_cm_per_s']
    unit_density = unit_mass / unit_length**3

    # get gas properties
    coords = np.array(f['PartType0']['Coordinates'])
    density = np.array(f['PartType0']['Density'])
    internal_energy = np.array(f['PartType0']['InternalEnergy'])
    xH2, xHp, xCO = np.array(f['PartType0']['ChemicalAbundances']).T

    boxsize = f['Header'].attrs['BoxSize']

    f.close()

    # get temperature
    ndensity = density * unit_density / ((1.0 + 4.0 * xHe) * constants.m_p.cgs.value)
    xTOT = 1.0 + xHp - xH2 + xHe
    nTOT = xTOT * ndensity
    mean_molecular_weight = density * unit_density / nTOT
    temperature = (2.0/3.0) * internal_energy * unit_velocity**2 * mean_molecular_weight / constants.k_B.cgs.value

    radii_in_kpc = np.array(radii_in_kpc)
    half_heights_in_kpc = np.array(half_heights_in_kpc)

    # get region outside disk
    CGM = temp_mask(temperature, temp_threshold_K) & cgm_mask(coords, boxsize, unit_length_in_kpc, radii_in_kpc, half_heights_in_kpc)
    
    print(f"CGM mask: {np.sum(CGM)} particles outside galactic disk over temperature threshold {temp_threshold_K} K")

    # update temperature
    internal_energy[CGM] = (3.0/2.0) * CGM_temp_K * constants.k_B.cgs.value / (unit_velocity**2 * mean_molecular_weight[CGM])

    # set temperature
    update_snapshot_property(filepath, 'PartType0', 'InternalEnergy', internal_energy)

