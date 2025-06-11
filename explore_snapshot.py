# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 12:42:00 2024

@author: zoefaes
"""

# Imports
import numpy as np
import astropy.units as u
from astropy import constants
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.animation as animation
import h5py
import warnings

class Snapshot:
    def __init__(self, filepath):
        # Read snapshot
        file = h5py.File(filepath, 'r')

        self.filepath = filepath
        self.name = filepath.split('/')[-1][0:-5]
        self.nsnap = self.name.split('_')[-1]

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

        # get snapshot info
        self.time = (file['Header'].attrs['Time'] << self.arepo_time).to(u.Myr)
        self.boxsize = file['Header'].attrs['BoxSize']
        self.parameters = file['Parameters']
        self.config_flags = file['Config']

        # gas data
        if 'PartType0' in file:
            self.has_type_0 = True
            self.ntype0 = len(file['PartType0']['ParticleIDs'])
            self.coordinates = file['PartType0']['Coordinates'] << self.arepo_length
            self.density = file['PartType0']['Density'] << self.arepo_density
            self.mass = file['PartType0']['Masses'] << self.arepo_mass
            self.velocity = file['PartType0']['Velocities'] << self.arepo_velocity
            
            if 'MagneticField' in file['PartType0']:
                self.has_mag_field = True
                self.mag_field = file['PartType0']['MagneticField'] << (self.arepo_mass ** (1/2) / (self.arepo_time * (self.arepo_length ** (1/2))))
                self.mag_field_divergence = file['PartType0']['MagneticFieldDivergence']
                self.mag_field_divergence_alternative = file['PartType0']['MagneticFieldDivergenceAlternative']
            else:
                self.has_mag_field = False
                self.mag_field = None
                self.mag_field_divergence = None
                self.mag_field_divergence_alternative = None

            if 'InternalEnergy' in file['PartType0']:
                self.has_internal_energy = True
                self.internal_energy = file['PartType0']['InternalEnergy'] << self.arepo_energy/self.arepo_mass
            else:
                self.has_internal_energy = False
                self.internal_energy = None

            if 'ChemicalAbundances' in file['PartType0']:
                self.has_chemistry = True
                self.xHe = 0.1 # Assumed in chemical network?
                self.xH2 = file['PartType0']['ChemicalAbundances'].T[0]
                self.xHp = file['PartType0']['ChemicalAbundances'].T[1]
                self.xCO = file['PartType0']['ChemicalAbundances'].T[2]
            else:
                self.has_chemistry = False
                self.xHe = None
                self.xH2 = None
                self.xHp = None
                self.xCO = None

            if 'Potential' in file['PartType0']:
                self.has_potential = True
                self.potential = file['PartType0']['Potential'] << self.arepo_energy/self.arepo_mass
            else:
                self.has_potential = False
                self.potential = None

            # get derived quantities
            self.ndensity = self.get_number_density()
            self.temperature = self.get_temperature()
            # use reasonable units for division to prevent overflow runtime warning
            self.cell_volume = self.mass.to_value(u.solMass) / self.density.to_value(u.solMass/u.pc**3)
            self.cell_volume = self.cell_volume * u.pc**3
            self.effective_cell_radius = (3 * self.cell_volume / (4 * np.pi)) ** (1/3)


        else:
            self.has_type_0 = False
            print('Could not find particles of type 0 (typically, gas particles) in snapshot.')

        # dark matter data
        if 'PartType1' in file:
            self.has_type_1 = True
            self.ntype1 = len(file['PartType1']['ParticleIDs'])
            self.coordinates_type1 = file['PartType1']['Coordinates'] << self.arepo_length
            self.velocity_type1 = file['PartType1']['Velocities'] << self.arepo_velocity
        else:
            self.has_type_1 = False

        # disk data
        if 'PartType2' in file:
            self.has_type_2 = True
            self.ntype2 = len(file['PartType2']['ParticleIDs'])
            self.coordinates_type2 = file['PartType2']['Coordinates'] << self.arepo_length
            self.velocity_type2 = file['PartType2']['Velocities'] << self.arepo_velocity
        else:
            self.has_type_2 = False

        # bulge data
        if 'PartType3' in file:
            self.has_type_3 = True
            self.ntype3 = len(file['PartType3']['ParticleIDs'])
            self.coordinates_type3 = file['PartType3']['Coordinates'] << self.arepo_length
            self.velocity_type3 = file['PartType3']['Velocities'] << self.arepo_velocity
        else:
            self.has_type_3 = False

        # star data
        if 'PartType4' in file:
            self.has_type_4 = True
            self.ntype4 = len(file['PartType4']['ParticleIDs'])
            self.coordinates_type4 = file['PartType4']['Coordinates'] << self.arepo_length
            self.velocity_type4 = file['PartType4']['Velocities'] << self.arepo_velocity
            self.mass_type4 = file['PartType4']['Masses'] << self.arepo_mass
        else:
            self.has_type_4 = False

        # sink data
        if 'PartType5' in file:
            self.has_type_5 = True
            self.ntype5 = len(file['PartType5']['ParticleIDs'])
            self.coordinates_type5 = file['PartType5']['Coordinates'] << self.arepo_length
            self.velocity_type5 = file['PartType5']['Velocities'] << self.arepo_velocity
        else:
            self.has_type_5 = False

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
        density = self.density.to(u.g/u.cm**3)
        xHe = 0.1 if (self.xHe is None) else self.xHe
        ndensity = density/((1.0 + 4.0 * xHe) * constants.m_p.cgs)

        return ndensity
    
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
        nTOT = xTOT * self.ndensity
        mean_molecular_weight = self.density/nTOT
        # print(f"Mean Molecular Weight: \nMean: {np.mean((mean_molecular_weight/constants.m_p).decompose())}, \nMax: {np.max((mean_molecular_weight/constants.m_p).decompose())}, \nMin: {np.min((mean_molecular_weight/constants.m_p).decompose())}")
        temperature = ((2.0/3.0) * self.internal_energy * mean_molecular_weight / constants.k_B).to(u.K)
        
        return temperature

    def check_particle_coordinate_limits(self) -> bool:
        """
        Check that all coordinates lie within the simulation box for each particle type.
        """
        success = True
        if self.has_type_0:
            if not np.all(np.logical_and(self.coordinates.to_value(self.arepo_length) >= 0, 
                                         self.coordinates.to_value(self.arepo_length) <= self.boxsize)):
                warnings.warn('Type 0 particle coordinates are not within the simulation box.')
                success = False
        if self.has_type_1:
            if not np.all(np.logical_and(self.coordinates_type1.to_value(self.arepo_length) >= 0, 
                                         self.coordinates_type1.to_value(self.arepo_length) <= self.boxsize)):
                warnings.warn('Type 1 particle coordinates are not within the simulation box.')
                success = False
        if self.has_type_2:
            if not np.all(np.logical_and(self.coordinates_type2.to_value(self.arepo_length) >= 0, 
                                         self.coordinates_type2.to_value(self.arepo_length) <= self.boxsize)):
                warnings.warn('Type 2 particle coordinates are not within the simulation box.')
                success = False
        if self.has_type_3:
            if not np.all(np.logical_and(self.coordinates_type3.to_value(self.arepo_length) >= 0, 
                                         self.coordinates_type3.to_value(self.arepo_length) <= self.boxsize)):
                warnings.warn('Type 3 particle coordinates are not within the simulation box.')
                success = False
        if self.has_type_4:
            if not np.all(np.logical_and(self.coordinates_type4.to_value(self.arepo_length) >= 0, 
                                         self.coordinates_type4.to_value(self.arepo_length) <= self.boxsize)):
                warnings.warn('Type 4 particle coordinates are not within the simulation box.')
                success = False
        if self.has_type_5:
            if not np.all(np.logical_and(self.coordinates_type5.to_value(self.arepo_length) >= 0, 
                                         self.coordinates_type5.to_value(self.arepo_length) <= self.boxsize)):
                warnings.warn('Type 5 particle coordinates are not within the simulation box.')
                success = False
        if success:
            print('All particle coordinates are within the simulation box.')

        return success
        
    def get_centered_coordinates(self, parttype=0) -> np.ndarray:
        """
        Re-centers simulation box such that the center of the simulation box lies at (0,0,0).

        Returns
        -------

        coordinates: np.ndarray
            re-centered coordinates [x, y, z].
        """
        center = (self.boxsize / 2) << self.arepo_length

        if parttype == 0:
            coordinates = self.coordinates - center
        elif parttype == 1:
            coordinates = self.coordinates_type1 - center
        elif parttype == 2:
            coordinates = self.coordinates_type2 - center
        elif parttype == 3:
            coordinates = self.coordinates_type3 - center
        elif parttype == 4:
            coordinates = self.coordinates_type4 - center
        elif parttype == 5:
            coordinates = self.coordinates_type5 - center
        else:
            raise ValueError('Invalid particle type. Choose from 0, 1, 2, 3, 4, or 5.')

        return coordinates
    
    def set_disk(self, radius = 10 * u.kpc, half_height = 1 * u.kpc):

        radius = radius.to_value(u.pc)
        half_height = half_height.to_value(u.pc)
        x, y, z = self.get_centered_coordinates().to_value(u.pc).T
        r = np.sqrt(x*x + y*y)
        theta = np.arctan2(y,x)

        disk = (r < radius) & (np.abs(z) < half_height)

        self.disk = disk


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
    
    return f'Mean: {mean:.2f}, Std Dev: {std_dev:.2f}, Median: {median:.2f}, Min: {min_val:.2f}, Max: {max_val:.2f}, First Quartile: {first_quartile:.2f}, Third Quartile: {third_quartile:.2f}'