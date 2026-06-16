# -*- coding: utf-8 -*-
"""
This file contains methods to read and modify AREPO snapshots in HDF5 format.
The Snapshot class contains tools to facilitate the computation of useful derived quantities.

Created on Wed Oct 30 12:42:00 2024

@author: zoefaes
"""

# Imports
from matplotlib import path
import numpy as np
import astropy.units as u
from astropy import constants
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.animation as animation
import h5py
import warnings
import time
from pathlib import Path
import os
import re

class Snapshot:
    def __init__(self, filepath):
        # Read snapshot
        t0 = time.time()
        file = h5py.File(filepath, 'r')

        self.filepath = filepath
        path = Path(self.filepath)
        self.name = path.stem
        self.nsnap = self.name.split('_')[-1]
        self.model = self.name.split('_')[0]

        # set units
        try:
            self.arepo_length = u.def_unit('arepo_length', file.get('Header')['UnitLength_in_cm'] << u.cm)
        except KeyError:
            self.arepo_length = u.def_unit('arepo_length', file['Parameters'].attrs['UnitLength_in_cm'] << u.cm)
        try:
            self.arepo_mass = u.def_unit('arepo_mass', file.get('Header')['UnitMass_in_g'] << u.g)
        except KeyError:
            self.arepo_mass = u.def_unit('arepo_mass', file['Parameters'].attrs['UnitMass_in_g'] << u.g)
        try:
            self.arepo_velocity = u.def_unit('arepo_velocity', file.get('Header')['UnitVelocity_in_cm_per_s'] << (u.cm/u.s))
        except KeyError:
            self.arepo_velocity = u.def_unit('arepo_velocity', file['Parameters'].attrs['UnitVelocity_in_cm_per_s'] << (u.cm/u.s))
        self.arepo_time = u.def_unit('arepo_time', self.arepo_length / self.arepo_velocity)
        self.arepo_density = u.def_unit('arepo_density', self.arepo_mass / self.arepo_length / self.arepo_length / self.arepo_length)
        self.arepo_energy = u.def_unit('arepo_energy', self.arepo_mass * self.arepo_velocity * self.arepo_velocity)
        self.arepo_mag = u.def_unit('arepo_mag', (self.arepo_mass ** (1/2) / (self.arepo_time * (self.arepo_length ** (1/2)))))

        # get snapshot info
        self.time = (file['Header'].attrs['Time'] << self.arepo_time).to(u.Myr)
        self.boxsize = file['Header'].attrs['BoxSize'][()]
        self.parameters = file['Parameters']
        self.config_flags = file['Config']

        # gas data
        if 'PartType0' in file:
            self.has_type_0 = True
            self.n_type_0 = file['PartType0']['ParticleIDs'].shape[0]

            self.coordinates = file['PartType0']['Coordinates'][()] # << self.arepo_length
            self.density = file['PartType0']['Density'][()] # << self.arepo_density
            self.mass = file['PartType0']['Masses'][()] # << self.arepo_mass
            self.velocity = file['PartType0']['Velocities'][()] # << self.arepo_velocity

            if 'MagneticField' in file['PartType0']:
                self.has_mag_field = True
                self.mag_field = file['PartType0']['MagneticField'][()] # << (self.arepo_mass ** (1/2) / (self.arepo_time * (self.arepo_length ** (1/2))))
                try:
                    self.mag_field_divergence = file['PartType0']['MagneticFieldDivergence']
                    self.mag_field_divergence_alternative = file['PartType0']['MagneticFieldDivergenceAlternative']
                except KeyError:
                    warnings.warn("No MagneticFieldDivergence or MagneticFieldDivergenceAlternative found in snapshot.")
            else:
                self.has_mag_field = False
                self.mag_field = None
                self.mag_field_divergence = None
                self.mag_field_divergence_alternative = None

            if 'InternalEnergy' in file['PartType0']:
                self.has_internal_energy = True
                self.internal_energy = file['PartType0']['InternalEnergy'][()] # << self.arepo_energy/self.arepo_mass
            else:
                self.has_internal_energy = False
                self.internal_energy = None

            if 'ChemicalAbundances' in file['PartType0']:
                self.has_chemistry = True
                self.xHe = 0.1 # Assumed in chemical network
                chem = file['PartType0']['ChemicalAbundances'][:]
                if chem.shape[1] == 3:
                    self.xH2, self.xHp, self.xCO = chem.T
                    self.chem_model = 'NL97'
                elif chem.shape[1] == 9:
                    self.chem_model = 'Gong'
                    self.xH2, self.xHp, self.xCp, self.xCH, self.xOH, self.xCO, self.xHCOp, self.xHep, self.xSip = chem.T
                else:
                    self.chem_model = 'Unknown'


            else:
                self.has_chemistry = False
                self.xHe = None
                self.xH2 = None
                self.xHp = None
                self.xCO = None

            if 'Potential' in file['PartType0']:
                self.has_potential = True
                self.potential = file['PartType0']['Potential'][()] # << self.arepo_energy/self.arepo_mass
            else:
                self.has_potential = False
                self.potential = None

            if ('TracerField' in file['PartType0']) or ('ZoomTracer' in file['PartType0']):
                self.has_tracer_field = True
                try:
                    self.tracer_field = file['PartType0']['ZoomTracer'][()]  # Eagerly read (i.e. load into memory)
                    try:
                        self.zoom_level = file['PartType0']['ZoomLevel'][()]
                    except KeyError:
                        self.zoom_level = None
                except KeyError:
                    self.tracer_field = file['PartType0']['TracerField'][()]  # Eagerly read
            else:
                self.has_tracer_field = False
                self.tracer_field = None

            if 'ZoomLevel' in file['PartType0']:
                self.has_zoom_level = True
                self.zoom_level = file['PartType0']['ZoomLevel'][()]
            else:
                self.has_zoom_level = False
                self.zoom_level = None

            # get derived quantities
            self.ndensity = self.get_number_density()

            self.temperature = self.get_temperature()

            # use reasonable units for division to prevent overflow runtime warning
            self.cell_volume = self.mass / self.density # << (self.arepo_length**3)
            self.effective_cell_radius = ((3 * self.cell_volume / (4 * np.pi)) ** (1/3)) # << self.arepo_length

        else:
            self.has_type_0 = False
            print('Could not find particles of type 0 (typically, gas particles) in snapshot.')

        # dark matter data
        if 'PartType1' in file:
            self.has_type_1 = True
            self.n_type_1 = file['PartType1']['ParticleIDs'].shape[0]
            self.coordinates_type_1 = file['PartType1']['Coordinates'] # << self.arepo_length
            self.velocity_type_1 = file['PartType1']['Velocities'] # << self.arepo_velocity
            self.mass_type_1 = file['PartType1']['Masses'] # << self.arepo_mass
        else:
            self.has_type_1 = False

        # disk data
        if 'PartType2' in file:
            self.has_type_2 = True
            self.n_type_2 = file['PartType2']['ParticleIDs'].shape[0]
            self.coordinates_type_2 = file['PartType2']['Coordinates'] # << self.arepo_length
            self.velocity_type_2 = file['PartType2']['Velocities'] # << self.arepo_velocity
            self.mass_type_2 = file['PartType2']['Masses'][()] # << self.arepo_mass
        else:
            self.has_type_2 = False

        # bulge data
        if 'PartType3' in file:
            self.has_type_3 = True
            self.n_type_3 = file['PartType3']['ParticleIDs'].shape[0]
            self.coordinates_type_3 = file['PartType3']['Coordinates'] # << self.arepo_length
            self.velocity_type_3 = file['PartType3']['Velocities'] # << self.arepo_velocity
            self.mass_type_3 = file['PartType3']['Masses'][()] # << self.arepo_mass
        else:
            self.has_type_3 = False

        # star data
        if 'PartType4' in file:
            self.has_type_4 = True
            self.n_type_4 = file['PartType4']['ParticleIDs'].shape[0]
            self.coordinates_type_4 = file['PartType4']['Coordinates'][()] # << self.arepo_length
            self.velocity_type_4 = file['PartType4']['Velocities'] # << self.arepo_velocity
            self.mass_type_4 = file['PartType4']['Masses'][()] # << self.arepo_mass

            if 'TracerField' in file['PartType4']:
                self.has_tracer_field_type_4 = True
                self.tracer_field_type_4 = file['PartType4']['TracerField'][()]  # Eagerly read (i.e. load into memory)
            else:
                self.has_tracer_field_type_4 = False
                self.tracer_field_type_4 = None
        else:
            self.has_type_4 = False

        # sink data
        if 'PartType5' in file:
            self.has_type_5 = True
            self.n_type_5 = file['PartType5']['ParticleIDs'].shape[0]
            self.coordinates_type_5 = file['PartType5']['Coordinates'][()] # << self.arepo_length
            self.velocity_type_5 = file['PartType5']['Velocities'] # << self.arepo_velocity
            self.mass_type_5 = file['PartType5']['Masses'] # << self.arepo_mass
        else:
            self.has_type_5 = False

        t1 = time.time()
        print(f"Time taken to read snapshot: {t1-t0:.2f} seconds")

        file.close()
    
    # The __repr__ method for a formal representation
    def __repr__(self):
        return f"Format class attribute nicely here."
    
    # The __str__ method for a more user-friendly representation
    def __str__(self):
        return f"Summarize snapshot properties here."
    
    def get_number_density(self):
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
        xHe = 0.1 if (self.xHe is None) else self.xHe
        ndensity = self.density/((1.0 + 4.0 * xHe) * constants.m_p.to_value(self.arepo_mass))
        if self.has_chemistry:
            xTOT = 1.0 + self.xHp - self.xH2 + xHe
        else:
            warnings.warn('No chemistry data available. Temperature calculated assuming neutral atomic gas with 0.1 He/H.', UserWarning)
            xTOT = 1.0 + 0.1
        self.mean_molecular_weight = (1.0 + 4.0 * xHe) * constants.m_p.to_value(self.arepo_mass) / xTOT

        return ndensity # << (1/(self.arepo_length**3))
    
    def get_temperature(self):
        """        
        Calculates temperature using
        :math:`T = \\frac{2}{3}\\frac{u_\\mathrm{therm} m_\\mathrm{avg} m_\\mathrm{p}}{k_\\mathrm{b}}`.

        Parameters
        ----------

        rho :  `~astropy.units.Quantity`
            Density :math:`\\rho`.

        utherm :  `~astropy.units.Quantity`
            Thermal energy per unit mass :math:`u_\\mathrm{therm}`.

        nTOT : `~astropy.units.Quantity`
            Total number density :math:`n_\\mathrm{tot}`.

        Returns
        -------

        temperature :  `~astropy.units.Quantity` or None
            Temperature :math:`T`.

        """
        if self.has_chemistry:
            xTOT = 1.0 + self.xHp - self.xH2 + self.xHe
        else:
            warnings.warn('No chemistry data available. Temperature calculated assuming neutral atomic gas with 0.1 He/H.', UserWarning)
            xTOT = 1.0 + 0.1  # Default value if no chemistry data is available
        nTOT = xTOT * self.ndensity # .to_value(1/(u.cm**3))
        mean_molecular_weight = self.density / nTOT
        temperature = (2.0/3.0) * (self.internal_energy * mean_molecular_weight) / constants.k_B.to_value(self.arepo_energy/u.K)

        return temperature # << u.K
    

    def check_particle_coordinate_limits(self) -> bool:
        """
        Check that all coordinates lie within the simulation box for each particle type.
        """
        success = True
        if self.has_type_0:
            if not np.all(np.logical_and(self.coordinates >= 0, 
                                         self.coordinates <= self.boxsize.to_value(self.arepo_length))):
                warnings.warn('Type 0 particle coordinates are not within the simulation box.')
                success = False
        if self.has_type_1:
            if not np.all(np.logical_and(self.coordinates_type_1 >= 0, 
                                         self.coordinates_type_1 <= self.boxsize.to_value(self.arepo_length))):
                warnings.warn('Type 1 particle coordinates are not within the simulation box.')
                success = False
        if self.has_type_2:
            if not np.all(np.logical_and(self.coordinates_type_2 >= 0, 
                                         self.coordinates_type_2 <= self.boxsize.to_value(self.arepo_length))):
                warnings.warn('Type 2 particle coordinates are not within the simulation box.')
                success = False
        if self.has_type_3:
            if not np.all(np.logical_and(self.coordinates_type_3 >= 0, 
                                         self.coordinates_type_3 <= self.boxsize.to_value(self.arepo_length))):
                warnings.warn('Type 3 particle coordinates are not within the simulation box.')
                success = False
        if self.has_type_4:
            if not np.all(np.logical_and(self.coordinates_type_4 >= 0, 
                                         self.coordinates_type_4 <= self.boxsize.to_value(self.arepo_length)    )):
                warnings.warn('Type 4 particle coordinates are not within the simulation box.')
                success = False
        if self.has_type_5:
            if not np.all(np.logical_and(self.coordinates_type_5 >= 0, 
                                         self.coordinates_type_5 <= self.boxsize.to_value(self.arepo_length))):
                warnings.warn('Type 5 particle coordinates are not within the simulation box.')
                success = False
        if success:
            print('All particle coordinates are within the simulation box.')

        return success
        
    def get_centered_coordinates(self, center=None, part_type=0) -> np.ndarray:
        """
        Re-centers simulation box such that the center of the simulation box lies at specified center or box center if center is None.

        Returns
        -------

        coordinates: np.ndarray
            re-centered coordinates [x, y, z].
        """
        if center is None:
            center = (self.boxsize / 2)
        else:
            if not all(isinstance(c, u.Quantity) and c.unit.is_length_unit() for c in center):
                raise ValueError('Center must be an array of astropy Quantities with length units.')
            center = np.array([c.to_value(self.arepo_length) for c in center])

        if part_type == 0:
            coordinates = (self.coordinates - center) # << self.arepo_length
        elif part_type == 1:
            coordinates = (self.coordinates_type_1 - center) # << self.arepo_length
        elif part_type == 2:
            coordinates = (self.coordinates_type_2 - center) # << self.arepo_length
        elif part_type == 3:
            coordinates = (self.coordinates_type_3 - center) # << self.arepo_length
        elif part_type == 4:
            coordinates = (self.coordinates_type_4 - center) # << self.arepo_length
        elif part_type == 5:
            coordinates = (self.coordinates_type_5 - center) # << self.arepo_length
        else:
            raise ValueError('Invalid particle type. Choose from 0, 1, 2, 3, 4, or 5.')

        return coordinates
    
    def set_disk(self, radius = 10 * u.kpc, half_height = 1 * u.kpc, part_type=0):

        if radius is None:
            radius = self.boxsize/2
        else:
            radius = radius.to_value(self.arepo_length)
        
        if half_height is None:
            half_height = self.boxsize/2
        else:
            half_height = half_height.to_value(self.arepo_length)
            
        x, y, z = self.get_centered_coordinates(part_type=part_type).T
        r = np.sqrt(x*x + y*y)
        theta = np.arctan2(y,x)

        disk = (r < radius) & (np.abs(z) < half_height)

        self.disk = disk

    def set_box(self, xmax=10 * u.kpc, ymax=10 * u.kpc, zmax=10 * u.kpc, part_type=0):

        xmax = xmax.to_value(self.arepo_length)
        ymax = ymax.to_value(self.arepo_length)
        zmax = zmax.to_value(self.arepo_length)
        x, y, z = self.get_centered_coordinates(part_type=part_type).T

        box = (np.abs(x) < xmax) & (np.abs(y) < ymax) & (np.abs(z) < zmax)

        self.box = box


def update_snapshot_property(filepath: str,
                             snap_group: str,
                             property: str,
                             data):
    # write to hdf5 file
    with h5py.File(filepath, 'r+') as f:
        if snap_group in f.keys():
            if property in f[snap_group]:
                if not f[snap_group][property].shape[0] == data.shape[0]:
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

                    print(f"{property} successfully updated.")
            else:
                raise Exception(f"{property} not found in {snap_group}.")
        else:
            raise Exception(f"{snap_group} not found in the HDF5 file.")
        

def add_snapshot_property(filepath: str,
                          snap_group: str,
                          property: str,
                          data):
    # write to hdf5 file
    with h5py.File(filepath, 'a') as f:
        if snap_group in f.keys():

            if snap_group[0:8] == 'PartType':
                if data.shape[0] != f.get(snap_group)['ParticleIDs'].shape[0]:
                    warnings.warn(f"Length of {property} data does not match number of particles for {snap_group}.", UserWarning)
                    
            f[snap_group].create_dataset(property, data=data)
            print(f"{property} successfully added to {snap_group}.")
        else:
            raise Exception(f"{snap_group} not found in the HDF5 file.")


def remove_snapshot_property(filepath: str,
                             snap_group: str,
                             property: str):
    # write to hdf5 file
    with h5py.File(filepath, 'r+') as f:
        if snap_group in f.keys():
            if property in f[snap_group]:
                del f[snap_group][property]
            else:
                raise Exception(f"{property} not found in {snap_group}.")
        else:
            raise Exception(f"{snap_group} not found in the HDF5 file.")
        if property in f[snap_group]:
            print(f"Failed to remove {property} from {snap_group}.")
        else:
            print(f"{property} successfully removed from {snap_group}.")
        

def stats_description(array):
    """
    Returns a string with the mean, median, min, and max of an array.
    
    Parameters
    ----------
    array : np.ndarray
        The array to calculate statistics for.
    
    Returns
    -------
    str
        A string with the statistics.
    """
    mean = np.mean(array)
    median = np.median(array)
    min_val = np.min(array)
    max_val = np.max(array)
    first_quartile = np.percentile(array, 25)
    third_quartile = np.percentile(array, 75)
    std_dev = np.std(array)
    
    return f'Mean: {mean}, Std Dev: {std_dev}, Median: {median}, Min: {min_val}, Max: {max_val}, First Quartile: {first_quartile}, Third Quartile: {third_quartile}'


def print_snap_info(filepath: str):
    """
    Prints the snapshot information from the HDF5 file.
    
    Parameters
    ----------
    filepath : str
        The path to the HDF5 file.
    """
    with h5py.File(filepath, 'r') as file:
        if 'Config' not in file.keys():
            print("\nNo 'Config' group found in the file.")
        else:
            print("\nCONFIG")
            for item in file['Config'].attrs:
                print(item, file['Config'].attrs[item])

        print("\nHEADER")
        for item in file['Header'].attrs:
            print(item, file['Header'].attrs[item])

        if 'Parameters' not in file.keys():
            print("\nNo 'Parameters' group found in the file.")
        else:
            print("\nPARAMETERS")
            for item in file['Parameters'].attrs:
                print(item, file['Parameters'].attrs[item])

        for key in file.keys():
            if key[0:8] == 'PartType':  # Skip 'Config', 'Header', and 'Parameters'
                try:
                    print(f"\n{key} keys: \n{file[key].keys()}")
                except AttributeError:
                    print(f"\n{key} doesn't have keys")

        if (np.array(file['PartType0']['ParticleIDs']) == 0).any():
            print("PartType0 ParticleID is 0 at index: ", np.where(np.array(file['PartType0']['ParticleIDs']) == 0))

        if (np.array(file['PartType0']['Masses']) == 0).any():
            print(f"PartType0 Mass is 0 for {np.sum(np.where(np.array(file['PartType0']['Masses']) == 0))} cells.")

            zero_mask = np.where(np.array(file['PartType0']['Masses']) == 0)
            coords = np.array(file['PartType0']['Coordinates'])

            # Plot XY and XZ planes
            plt.figure(figsize=(12, 6))
            plt.subplot(1, 2, 1)
            plt.scatter(coords[zero_mask, 0], coords[zero_mask, 1], s=2, color="red")
            plt.xlabel("X")
            plt.ylabel("Y")
            plt.title("Cells with mass zero (XY plane)")
            plt.axis("equal")

            plt.subplot(1, 2, 2)
            plt.scatter(coords[zero_mask, 0], coords[zero_mask, 2], s=2, color="red")
            plt.xlabel("X")
            plt.ylabel("Z")
            plt.title("Cells with mass zero (XZ plane)")
            plt.axis("equal")
            plt.show()

        if (np.array(file['PartType0']['Masses']) < 0).any():
            print(f"PartType0 Mass is negative for {np.sum(np.where(np.array(file['PartType0']['Masses']) < 0))} cells.")

        if (np.isnan(np.array(file['PartType0']['Masses']))).any():
            print(f"PartType0 Mass is NaN for {np.sum(np.where(np.isnan(np.array(file['PartType0']['Masses']))))} cells.")

        if (np.isinf(np.array(file['PartType0']['Masses']))).any():
            print(f"PartType0 Mass is Inf for {np.sum(np.where(np.isinf(np.array(file['PartType0']['Masses']))))} cells.")

        file.close()

def get_max_snapshot_number(rundir):

    # Regex to capture numbers in filenames
    pattern = re.compile(r"V\d+B?_(\d{3,4})\.hdf5")

    if not os.path.isdir(rundir):
        print(f"Directory {rundir} does not exist.")
    
    nums = []
    for fname in os.listdir(rundir):
        match = pattern.match(fname)
        if match:
            nums.append(int(match.group(1)))

    return max(nums) if nums else None