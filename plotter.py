# -*- coding: utf-8 -*-
"""
This file contains plotting functions compatible with Snapshot objects.

Created on Fri May 30 11:43:00 2025

@author: zoefaes
"""

# Imports
import numpy as np
import astropy.units as u
from astropy import constants
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import matplotlib.animation as animation
import h5py
import warnings
from snapshot_tools import *


####################################################################################################
#                                         Helper Functions                                         #
####################################################################################################

def scale_array(arr, min_val, max_val):
    arr_min = np.min(arr)
    arr_max = np.max(arr)
    return min_val + (arr - arr_min) * (max_val - min_val) / (arr_max - arr_min)


def compute_galactic_longitudes(positions, observer):
    """
    Compute Galactic longitudes (in degrees) for 3D coordinates
    relative to an observer.

    Parameters:
    -----------
    positions : np.ndarray
        Array of shape (N, 3) with positions (x, y, z).
        Assumes Galactic center is at (0,0,0).
    observer : tuple or array-like
        Observer position (x_obs, y_obs, z_obs).

    Returns:
    --------
    longitudes : np.ndarray
        Galactic longitudes in degrees, in range [0, 360).
    """
    pos = np.asarray(positions)
    obs = np.asarray(observer)

    # Vectors from observer to each point
    rel_pos = pos - obs  # shape (N, 3)

    x_rel = rel_pos[:, 0]
    y_rel = rel_pos[:, 1]

    # Galactic longitude: angle in XY plane from observer to point
    longitudes_rad = np.arctan2(y_rel, x_rel)  # radians, range [-π, π]
    longitudes_deg = np.degrees(longitudes_rad)

    # Normalize to [0, 360)
    longitudes_deg = np.mod(longitudes_deg, 360)

    return longitudes_deg


def compute_los_velocity(positions, velocities, observer):
    """
    Compute the line-of-sight (LOS) velocity for a set of particles relative to an observer.

    Parameters
    ----------
    positions : (N, 3) array
        Cartesian coordinates of particles (x, y, z).
    velocities : (N, 3) array
        Velocity vectors (vx, vy, vz) of the particles.
    observer : (3,) array-like
        Position of the observer (x_obs, y_obs, z_obs).

    Returns
    -------
    v_los : (N,) array
        Line-of-sight velocities (positive = moving away from observer).
    """
    pos = np.asarray(positions)
    vel = np.asarray(velocities)
    obs = np.asarray(observer)

    # Relative position vectors from observer to each particle
    rel_pos = pos - obs  # shape (N, 3)

    # Normalize to get unit line-of-sight direction vectors
    distance = np.linalg.norm(rel_pos, axis=1)
    unit_los = rel_pos / distance[:, np.newaxis]  # shape (N, 3)

    # Project velocities onto line-of-sight unit vector
    v_los = np.sum(vel * unit_los, axis=1)  # scalar projection

    return v_los


####################################################################################################
#                                         Plotting Functions                                       #
####################################################################################################


def mass_hist(snap,
                   axis = 'z',
                   disk_radius = None,
                   disk_half_height = None,
                   bins = 500,
                   norm = 'log',
                   vmin = None,
                   vmax = None,
                   cmap = 'viridis',
                   show_stars = False,
                   save_path = None):
    
    # get disk mask if applicable
    if (disk_radius is not None) & (disk_half_height is not None) & (axis != 'z'):
        raise ValueError('Disk mask can only be applied to z-axis projection.')
    if (disk_radius is not None) & (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)
    
    # set color limits
    if vmin is None:
        vmin = np.min(snap.mass.to_value(u.solMass))
    if vmax is None:
        vmax = np.max(snap.mass.to_value(u.solMass))

    if norm == 'log':
        if vmin <= 0:
            warnings.warn("LogNorm requires vmin to be positive. Setting vmin to 1e-10.", UserWarning)
            vmin = 1e-10  # avoid log(0) issues        
        if vmax <= 0:
            raise ValueError("LogNorm requires vmin and vmax to be positive.")
        norm = LogNorm(vmin=vmin, vmax=vmax)

    # Get coordinates
    x, y, z = snap.get_centered_coordinates().T.to_value(u.kpc)

    if show_stars:
        if snap.has_type_4:
            xs, ys, zs = snap.get_centered_coordinates(parttype=4).T.to_value(u.kpc)
        else:
            warnings.warn(f"Snapshot {snap.name} does not have star particles.")

    # Start figure
    fig = plt.figure(figsize=(10,8))
    
    if axis == 'z':
        if (disk_radius is not None) & (disk_half_height is not None): # apply disk mask if applicable
            hist2d = plt.hist2d(x[snap.disk], y[snap.disk], bins=bins, norm=norm, weights=snap.mass[snap.disk].to_value(u.solMass), cmap=cmap)#, cmin=10)
        else: # no disk mask
            hist2d = plt.hist2d(x, y, bins=bins, norm=norm, weights=snap.mass.to_value(u.solMass), cmap=cmap)#, cmin=10)
        if show_stars:
            plt.scatter(xs, ys, s=3, c='r', marker='o')
        plt.xlabel(r'X [kpc]')
        plt.ylabel(r'Y [kpc]')
    elif axis == 'x':
        hist2d = plt.hist2d(y, z, bins=bins, norm=norm, weights=snap.mass.to_value(u.solMass), cmap=cmap)#, cmin=10)
        if show_stars:
            plt.scatter(ys, zs, s=3, c='r', marker='o')
        plt.xlabel(r'Y [kpc]')
        plt.ylabel(r'Z [kpc]')
    elif axis == 'y':
        hist2d = plt.hist2d(x, z, bins=bins, norm=norm, weights=snap.mass.to_value(u.solMass), cmap=cmap)
        if show_stars:
            plt.scatter(xs, zs, s=3, c='r', marker='o')
        plt.xlabel(r'X [kpc]')
        plt.ylabel(r'Z [kpc]')
    else:
        raise ValueError('Invalid axis. Choose from x, y, or z.')
       
    # colorbar
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'Mass [M$_{\odot}$]')
    # title
    plt.title(f'Mass-Weighted {str(axis).upper()}-Projection of {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)', fontsize=16)
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_{axis}_mass_hist.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_{axis}_mass_hist.jpg")
    plt.show()


def column_density_hist(snap,
                        axis = 'z',
                        disk_radius = None,
                        disk_half_height = None,
                        bins = 500,
                        norm = 'log',
                        vmin = None,
                        vmax = None,
                        cmap = 'viridis',
                        show_stars = False,
                        save_path = None):
    
    # get disk mask if applicable
    if (disk_radius is not None) & (disk_half_height is not None) & (axis != 'z'):
        raise ValueError('Disk mask can only be applied to z-axis projection.')
    if (disk_radius is not None) & (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    # Get coordinates
    x, y, z = snap.get_centered_coordinates().T.to_value(u.kpc)
    if show_stars:
        if snap.has_type_4:
            xs, ys, zs = snap.get_centered_coordinates(parttype=4).T.to_value(u.kpc)
        else:
            warnings.warn(f"Snapshot {snap.name} does not have star particles.")

    # Get area per bin
    if (axis == 'z') & (disk_radius is not None) & (disk_half_height is not None):
        length = 2 * disk_half_height.to_value(u.cm)
    else:
        length = (snap.boxsize * snap.arepo_length).to_value(u.cm)

    # Calculate mass per area
    weights = snap.ndensity.to_value(1/u.cm**3) * length

    # set color limits
    if vmin is None:
        vmin = np.min(weights)
    if vmax is None:
        vmax = np.max(weights)

    if norm == 'log':
        if vmin <= 0:
            warnings.warn("LogNorm requires vmin to be positive. Setting vmin to 1e-10.", UserWarning)
            vmin = 1e-10  # avoid log(0) issues        
        if vmax <= 0:
            raise ValueError("LogNorm requires vmin and vmax to be positive.")
        norm = LogNorm(vmin=vmin, vmax=vmax)

    # Start figure
    fig = plt.figure(figsize=(10,8))
    
    if axis == 'z':
        if (disk_radius is not None) & (disk_half_height is not None): # apply disk mask if applicable
            hist2d = plt.hist2d(x[snap.disk], y[snap.disk], bins=bins, norm=norm, weights=weights[snap.disk], cmap=cmap)#, cmin=10)
        else: # no disk mask
            hist2d = plt.hist2d(x, y, bins=bins, norm=norm, weights=weights, cmap=cmap)#, cmin=10)
        if show_stars:
            plt.scatter(xs, ys, s=3, c='r', marker='o')
        plt.xlabel(r'X [kpc]')
        plt.ylabel(r'Y [kpc]')

    elif axis == 'x':
        hist2d = plt.hist2d(y, z, bins=bins, norm=norm, weights=weights, cmap=cmap)#, cmin=10)
        if show_stars:
            plt.scatter(ys, zs, s=3, c='r', marker='o')
        plt.xlabel(r'Y [kpc]')
        plt.ylabel(r'Z [kpc]')

    elif axis == 'y':
        hist2d = plt.hist2d(x, z, bins=bins, norm=norm, weights=weights, cmap=cmap)
        if show_stars:
            plt.scatter(xs, zs, s=3, c='r', marker='o')
        plt.xlabel(r'X [kpc]')
        plt.ylabel(r'Z [kpc]')
    else:
        raise ValueError('Invalid axis. Choose from x, y, or z.')
       
    # colorbar
    cbar = plt.colorbar(hist2d[3], ax=plt.gca())
    cbar.ax.set_ylabel(r'Column Density [cm$^{-2}$]')
    # title
    plt.title(f'Column Density {str(axis).upper()}-Projection {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)', fontsize=16)
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_{axis}_column_density.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_{axis}_column_density.jpg")
    plt.show()


def temperature_hist(snap,
                   axis = 'z',
                   disk_radius = None,
                   disk_half_height = None,
                   bins = 500,
                   norm = 'log',
                   vmin = None,
                   vmax = None,
                   cmap = 'viridis',
                   save_path = None):
    
    # get disk mask if applicable
    if (disk_radius is not None) & (disk_half_height is not None) & (axis != 'z'):
        raise ValueError('Disk mask can only be applied to z-axis projection.')
    if (disk_radius is not None) & (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)
    
    temperature = snap.temperature.to_value(u.K)

    inf_count = 0
    for temp in temperature:
        if np.isinf(temp):
            inf_count += 1
    print(f"Percentage of infinite temperature values: {inf_count / len(snap.temperature) * 100:.2f}% \nReplacing with 1e8 K")

    temperature[np.isinf(temperature)] = 1e8

    # set color limits
    if vmin is None:
        vmin = np.min(temperature)
    if vmax is None:
        vmax = np.max(temperature)

    if norm == 'log':
        if vmin <= 0:
            warnings.warn("LogNorm requires vmin to be positive. Setting vmin to 1e-10.", UserWarning)
            vmin = 1e-10  # avoid log(0) issues        
        if vmax <= 0:
            raise ValueError("LogNorm requires vmin and vmax to be positive.")
        norm = LogNorm(vmin=vmin, vmax=vmax)

    # Get coordinates
    x, y, z = snap.get_centered_coordinates().T.to_value(u.kpc)

    # Start figure
    fig = plt.figure(figsize=(10,8))
    
    if axis == 'z':
        if (disk_radius is not None) & (disk_half_height is not None): # apply disk mask if applicable
            hist2d = plt.hist2d(x[snap.disk], y[snap.disk], bins=bins, norm=norm, weights=temperature[snap.disk], cmap=cmap)#, cmin=10)
        else: # no disk mask
            hist2d = plt.hist2d(x, y, bins=bins, norm=norm, weights=temperature, cmap=cmap)#, cmin=10)
        plt.xlabel(r'X [kpc]')
        plt.ylabel(r'Y [kpc]')
    elif axis == 'x':
        hist2d = plt.hist2d(y, z, bins=bins, norm=norm, weights=temperature, cmap=cmap)#, cmin=10)
        plt.xlabel(r'Y [kpc]')
        plt.ylabel(r'Z [kpc]')
    elif axis == 'y':
        hist2d = plt.hist2d(x, z, bins=bins, norm=norm, weights=temperature, cmap=cmap)
        plt.xlabel(r'X [kpc]')
        plt.ylabel(r'Z [kpc]')
    else:
        raise ValueError('Invalid axis. Choose from x, y, or z.')
       
    # colorbar
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'Temperature [K]')
    # title
    plt.title(f'Temperature {str(axis).upper()}-Projection of {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)', fontsize=16)
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_{axis}_temperature_hist.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_{axis}_temperature_hist.jpg")
    plt.show()


def resolution(snap,
                    disk_radius = None,
                    disk_half_height = None,
                    bins = 500,
                    norm = 'log',
                    cmap = 'viridis',
                    save_path = None):
    
    if (disk_radius is not None) & (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    log_r_eff = np.log10(np.abs(snap.effective_cell_radius.to_value(u.pc)))
    log_ndens = np.log10(np.abs(snap.ndensity.to_value(u.cm**(-3))) + 1e-10)  # add small value to avoid log(0)

    if norm == 'log':
        norm = LogNorm()

    # Plot resolutiom
    fig, ax = plt.subplots(figsize=(10, 8))
    if (disk_radius is not None) & (disk_half_height is not None):
        hist2d = plt.hist2d(log_ndens[snap.disk], log_r_eff[snap.disk], cmap=cmap, bins=bins, norm=norm)
    else:
        hist2d = plt.hist2d(log_ndens, log_r_eff, cmap=cmap, bins=bins, norm=norm)

    # Add colorbar
    cbar = plt.colorbar(hist2d[3], ax=ax)
    cbar.ax.set_ylabel('Counts')
    plt.title(f'Spatial Resolution of {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)')
    plt.xlabel(r'Log(Number density) [cm$^{-3}$]')
    plt.ylabel(r'Log(Effective radius) [pc]')
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_resolution.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_resolution.jpg")
    plt.show()

def mass_resolution(snap,
                        disk_radius = None,
                        disk_half_height = None,
                        bins = 500,
                        norm = 'log',
                        cmap = 'viridis',
                        save_path = None):
    
    if (disk_radius is not None) & (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    log_r_eff = np.log10(snap.effective_cell_radius.to_value(u.pc))
    log_mass = np.log10(snap.mass.to_value(u.solMass) + 1e-10)  # add small value to avoid log(0)

    if norm == 'log':
        norm = LogNorm()

    # Plot resolutiom
    fig, ax = plt.subplots(figsize=(10, 8))
    if (disk_radius is not None) & (disk_half_height is not None):
        hist2d = plt.hist2d(log_mass[snap.disk], log_r_eff[snap.disk], cmap=cmap, bins=bins, norm=norm)# norm=mcolors.LogNorm())
    else:
        hist2d = plt.hist2d(log_mass, log_r_eff, cmap=cmap, bins=bins, norm=norm)# norm=mcolors.LogNorm())
    # Add colorbar
    cbar = plt.colorbar(hist2d[3], ax=ax)
    cbar.ax.set_ylabel('Counts')
    # cbar.ax.tick_params(labelsize=14)
    # plt.xticks(fontsize=14)
    # plt.yticks(fontsize=14)
    plt.title(f'Mass Resolution of {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)')
    plt.xlabel(r'Log(Mass) [M$_{\odot}$]')
    plt.ylabel(r'Log(Effective radius) [pc]')
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_mass_resolution.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_mass_resolution.jpg")
    plt.show()


def ndensity_hist(snap,
                   axis = 'z',
                   disk_radius = None,
                   disk_half_height = None,
                   bins = 500,
                   norm = 'log',
                   vmin = None,
                   vmax = None,
                   cmap = 'viridis',
                   save_path = None):
    
    # get disk mask if applicable
    if (disk_radius is not None) & (disk_half_height is not None) & (axis != 'z'):
        raise ValueError('Disk mask can only be applied to z-axis projection.')
    if (disk_radius is not None) & (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    # Get coordinates
    x, y, z = snap.get_centered_coordinates().T.to_value(u.kpc)

    # Get number density weights
    weights = snap.ndensity.to_value(1/u.cm**3)

    # set color limits
    if vmin is None:
        vmin = np.min(weights)
    if vmax is None:
        vmax = np.max(weights)

    if norm == 'log':
        if vmin <= 0:
            warnings.warn("LogNorm requires vmin to be positive. Setting vmin to 1e-10.", UserWarning)
            vmin = 1e-10  # avoid log(0) issues        
        if vmax <= 0:
            raise ValueError("LogNorm requires vmin and vmax to be positive.")
        norm = LogNorm(vmin=vmin, vmax=vmax)

    # Start figure
    fig = plt.figure(figsize=(10,8))
    
    if axis == 'z':
        if (disk_radius is not None) & (disk_half_height is not None): # apply disk mask if applicable
            hist2d = plt.hist2d(x[snap.disk], y[snap.disk], bins=bins, norm=norm, weights=weights[snap.disk], cmap=cmap)#, cmin=10)
        else: # no disk mask
            hist2d = plt.hist2d(x, y, bins=bins, norm=norm, weights=weights, cmap=cmap)#, cmin=10)
        plt.xlabel(r'X [kpc]')
        plt.ylabel(r'Y [kpc]')
    elif axis == 'x':
        hist2d = plt.hist2d(y, z, bins=bins, norm=norm, weights=weights, cmap=cmap)#, cmin=10)
        plt.xlabel(r'Y [kpc]')
        plt.ylabel(r'Z [kpc]')
    elif axis == 'y':
        hist2d = plt.hist2d(x, z, bins=bins, norm=norm, weights=weights, cmap=cmap)
        plt.xlabel(r'X [kpc]')
        plt.ylabel(r'Z [kpc]')
    else:
        raise ValueError('Invalid axis. Choose from x, y, or z.')
           
    # colorbar
    cbar = plt.colorbar(hist2d[3], ax=plt.gca())
    cbar.ax.set_ylabel(r'Number Density [cm$^{-3}$]')
    # title
    plt.title(f'Number Density-Weighted {str(axis).upper()}-Projection of {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)', fontsize=16)
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_{axis}_ndensity_hist.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_{axis}_ndensity_hist.jpg")
    plt.show()


def potential_hist(snap,
                   axis = 'z',
                   disk_radius = None,
                   disk_half_height = None,
                   bins = 500,
                   norm = 'log',
                   vmin = None,
                   vmax = None,
                   cmap = 'viridis',
                   save_path = None):
    
    # check if potential is available
    if not snap.has_potential:
        raise Exception('Snapshot does not have potential data.')
    
    # get disk mask if applicable
    if (disk_radius is not None) & (disk_half_height is not None) & (axis != 'z'):
        raise ValueError('Disk mask can only be applied to z-axis projection.')
    if (disk_radius is not None) & (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    # Get coordinates
    x, y, z = snap.get_centered_coordinates().T.to_value(u.kpc)

    # Get potential weights
    weights = np.abs(snap.potential.to_value((u.pc/u.yr)**2))

    # set color limits
    if vmin is None:
        vmin = np.min(weights)
    if vmax is None:
        vmax = np.max(weights)

    if norm == 'log':
        if vmin <= 0:
            warnings.warn("LogNorm requires vmin to be positive. Setting vmin to 1e-10.", UserWarning)
            vmin = 1e-10  # avoid log(0) issues        
        if vmax <= 0:
            raise ValueError("LogNorm requires vmin and vmax to be positive.")
        norm = LogNorm(vmin=vmin, vmax=vmax)
    elif norm == 'linear':
        norm = Normalize(vmin=vmin, vmax=vmax)

    # Start figure
    fig = plt.figure(figsize=(10,8))
    
    if axis == 'z':
        if (disk_radius is not None) & (disk_half_height is not None): # apply disk mask if applicable
            hist2d = plt.hist2d(x[snap.disk], y[snap.disk], bins=bins, norm=norm, weights=weights[snap.disk], cmap=cmap)#, cmin=10)
        else: # no disk mask
            hist2d = plt.hist2d(x, y, bins=bins, norm=norm, weights=weights, cmap=cmap)#, cmin=10)
        plt.xlabel(r'X [kpc]')
        plt.ylabel(r'Y [kpc]')
    elif axis == 'x':
        hist2d = plt.hist2d(y, z, bins=bins, norm=norm, weights=weights, cmap=cmap)#, cmin=10)
        plt.xlabel(r'Y [kpc]')
        plt.ylabel(r'Z [kpc]')
    elif axis == 'y':
        hist2d = plt.hist2d(x, z, bins=bins, norm=norm, weights=weights, cmap=cmap)
        plt.xlabel(r'X [kpc]')
        plt.ylabel(r'Z [kpc]')
    else:
        raise ValueError('Invalid axis. Choose from x, y, or z.')
    
    # colorbar
    cbar = plt.colorbar(hist2d[3], ax=plt.gca())
    cbar.ax.set_ylabel(r'Potential [pc$^2$ yr$^{-2}$]')

    # title
    plt.title(f'Potential-Weighted {str(axis).upper()}-Projection of {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)', fontsize=16)
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_{axis}_potential_hist.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_{axis}_potential_hist.jpg")
    plt.show()



def velocity_profile(snap,
                          disk_radius = None,
                          disk_half_height = None,
                          bins = 500,
                          norm = 'log',
                          cmap = 'viridis',
                          save_path = None):
    
    # TODO: add more velocity components & add velocity from other particle types

    if (disk_radius is not None) & (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)
                          
    # Get coordinates
    x, y, z = snap.get_centered_coordinates().T.to_value(u.kpc)
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    vx = snap.velocity[:, 0].to_value(u.km/u.s)
    vy = snap.velocity[:, 1].to_value(u.km/u.s)
    vth = vy * np.cos(theta) - vx * np.sin(theta) # Azimuthal velocity

    if norm == 'log':
        norm = LogNorm()

    fig, ax = plt.subplots(figsize=(10, 8))
    if (disk_radius is not None) & (disk_half_height is not None):
        hist2d = plt.hist2d(r[snap.disk], vth[snap.disk], cmap=cmap, bins=bins, norm=norm)# norm=mcolors.LogNorm())
    else:
        hist2d = plt.hist2d(r, vth, cmap=cmap, bins=bins, norm=norm)
    # Add colorbar
    cbar = plt.colorbar(hist2d[3], ax=ax)
    cbar.ax.set_ylabel('Counts')
    plt.title(f'Velocity Profile of {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)')
    plt.xlabel(r'Radius [kpc]')
    plt.ylabel(r'Azimuthal Velocity (V$_{\theta}$) [km/s]')
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_velocity_profile.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_velocity_profile.jpg")
    plt.show()


def ndensity_profile(snap,
                          disk_radius = None,
                          disk_half_height = None,
                          bins = 500,
                          norm = 'log',
                          cmap = 'viridis',
                          save_path = None):
    
    if (disk_radius is not None) & (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)
                          
    # Get coordinates
    x, y, z = snap.get_centered_coordinates().T.to_value(u.kpc)
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)

    if norm == 'log':
        norm = LogNorm()

    ndensity = snap.ndensity.to_value(1/u.cm**3)
    density = snap.density.to_value(u.g/u.cm**3)

    description = stats_description(ndensity)
    print(f"Number Density Statistics:\n{description}")

    description = stats_description(density)
    print(f"Density Statistics:\n{description}")

    fig, ax = plt.subplots(figsize=(10, 8))
    if (disk_radius is not None) & (disk_half_height is not None):
        hist2d = plt.hist2d(r[snap.disk], ndensity[snap.disk], cmap=cmap, bins=bins, norm=norm)# norm=mcolors.LogNorm())
    else:
        hist2d = plt.hist2d(r, ndensity, cmap=cmap, bins=bins, norm=norm)
    # Add colorbar
    cbar = plt.colorbar(hist2d[3], ax=ax)
    cbar.ax.set_ylabel('Counts')
    plt.title(f'Number Density Profile of {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)')
    plt.xlabel(r'Radius [kpc]')
    plt.ylabel(r'Number Density [cm$^{-3}$]')
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_ndensity_profile.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_ndensity_profile.jpg")
    plt.show()


def temperature_profile(snap,
                          disk_radius = None,
                          disk_half_height = None,
                          bins = 500,
                          norm = 'log',
                          cmap = 'viridis',
                          save_path = None):
                          
    if (disk_radius is not None) & (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    # Get coordinates
    x, y, z = snap.get_centered_coordinates().T.to_value(u.kpc)
    r = np.sqrt(x**2 + y**2)

    if norm == 'log':
        norm = LogNorm()

    log_temperature = np.log10(snap.temperature.to_value(u.K) + 1e-10)  # add small value to avoid log(0)

    inf_count = 0
    for temp in log_temperature:
        if np.isinf(temp):
            inf_count += 1

    if inf_count > 0:
        print(f"Percentage of infinite temperature values: {inf_count / len(snap.temperature) * 100:.2f}% \nReplacing inf values with 1e8 K")
        log_temperature[np.isinf(log_temperature)] = 8

    fig, ax = plt.subplots(figsize=(10, 8))
    if (disk_radius is not None) & (disk_half_height is not None):
        hist2d = plt.hist2d(r[snap.disk], log_temperature[snap.disk], cmap=cmap, bins=bins, norm=norm)# norm=mcolors.LogNorm())
    else:
        hist2d = plt.hist2d(r, log_temperature, cmap=cmap, bins=bins, norm=norm)

    # Add colorbar
    cbar = plt.colorbar(hist2d[3], ax=ax)
    cbar.ax.set_ylabel('Counts')
    plt.title(f'Temperature Profile of {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)')
    plt.xlabel(r'Radius [kpc]')
    plt.ylabel(r'log(Temperature) [K]')
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_temperature_profile.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_temperature_profile.jpg")
    plt.show()

def mass_profile(snap,
                disk_radius = None,
                disk_half_height = None,
                bins = 500,
                norm = 'log',
                cmap = 'viridis',
                save_path = None):
                          
    if (disk_radius is not None) & (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    # Get coordinates
    x, y, z = snap.get_centered_coordinates().T.to_value(u.kpc)
    r = np.sqrt(x**2 + y**2)

    if norm == 'log':
        norm = LogNorm()

    log_mass = np.log10(snap.mass.to_value(u.solMass))

    fig, ax = plt.subplots(figsize=(10, 8))
    if (disk_radius is not None) & (disk_half_height is not None):
        hist2d = plt.hist2d(r[snap.disk], log_mass[snap.disk], cmap=cmap, bins=bins, norm=norm)# norm=mcolors.LogNorm())
    else:
        hist2d = plt.hist2d(r, log_mass, cmap=cmap, bins=bins, norm=norm)
    # Add colorbar
    cbar = plt.colorbar(hist2d[3], ax=ax)
    cbar.ax.set_ylabel('Counts')
    plt.title(f'Mass Profile of {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)')
    plt.xlabel(r'Radius [kpc]')
    plt.ylabel(r'log(Mass) [M$_{\odot}$]')
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_mass_profile.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_mass_profile.jpg")
    plt.show()


def phase_space(snap,
                    disk_radius = None,
                    disk_half_height = None,
                    bins = 500,
                    norm = 'log',
                    cmap = 'viridis',
                    save_path = None):
    
    if (disk_radius is not None) & (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    temperature = snap.temperature.to_value(u.K)

    inf_count = 0
    for temp in temperature:
        if np.isinf(temp):
            inf_count += 1
    print(f"Percentage of infinite temperature values: {inf_count / len(snap.temperature) * 100:.2f}% \nReplacing with 1e8 K")

    temperature[np.isinf(temperature)] = 1e8

    log_temp = np.log10(temperature + 1e-10)  # add small value to avoid log(0)
    log_ndens = np.log10(snap.ndensity.to_value(u.cm**(-3)) + 1e-10)  # add small value to avoid log(0)

    if norm == 'log':
        norm = LogNorm()

    # Plot phase space
    fig, ax = plt.subplots(figsize=(10, 8))
    if (disk_radius is not None) & (disk_half_height is not None):
        hist2d = plt.hist2d(log_ndens[snap.disk], log_temp[snap.disk], cmap=cmap, bins=bins, norm=norm)# norm=mcolors.LogNorm())
    else:
        hist2d = plt.hist2d(log_ndens, log_temp, cmap=cmap, bins=bins, norm=norm)# norm=mcolors.LogNorm())
    # Add colorbar
    cbar = plt.colorbar(hist2d[3], ax=ax)
    cbar.ax.set_ylabel('Counts')
    plt.title(f'Phase Space of {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)')
    plt.xlabel(r'Log(Number Density) [cm$^{-3}$]')
    plt.ylabel(r'Log(Temperature) [K]')
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_phase_space.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_phase_space.jpg")
    plt.show()


def jeans_mass_hist(snap,
                    disk_radius = None,
                    disk_half_height = None,
                    gamma = 5/3,
                    bins = 500,
                    norm = 'log',
                    vmin = None,
                    vmax = None,
                    cmap = 'viridis',
                    save_path = None):
    
    # get disk mask if applicable
    if (disk_radius is not None) & (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)
    
    # Get coordinates
    x, y, z = snap.get_centered_coordinates().T.to_value(u.kpc)
    w = 5 * snap.internal_energy.to_value(u.J/u.kg) * (gamma - 1) / constants.G.si.value
    # internal energy [m**2/s**2], G [m**3/(kg s**2)], therefore w [kg/m]
    conversion_factor = (1*u.kg/u.m).to_value(u.solMass/u.kpc)
    w = w * conversion_factor  # convert to M_sun/kpc
    v = (3/(4 * np.pi * snap.density)).to_value(u.kpc**3/u.solMass)
    jeans_mass = w ** (3/2) * v ** (1/2)

    # set color limits
    if vmin is None:
        vmin = np.min(jeans_mass)
    if vmax is None:
        vmax = np.max(jeans_mass)

    if norm == 'log':
        if vmin <= 0:
            warnings.warn("LogNorm requires vmin to be positive. Setting vmin to 1e-10.", UserWarning)
            vmin = 1e-10  # avoid log(0) issues        
        if vmax <= 0:
            raise ValueError("LogNorm requires vmin and vmax to be positive.")
        norm = LogNorm(vmin=vmin, vmax=vmax)

    # Start figure
    fig = plt.figure(figsize=(10,8))
    
    if (disk_radius is not None) & (disk_half_height is not None): # apply disk mask if applicable
        hist2d = plt.hist2d(x[snap.disk], y[snap.disk], bins=bins, norm=norm, weights=jeans_mass[snap.disk], cmap=cmap)#, cmin=10)
    else: # no disk mask
        hist2d = plt.hist2d(x, y, bins=bins, norm=norm, weights=jeans_mass, cmap=cmap)#, cmin=10)
    plt.xlabel(r'X [kpc]')
    plt.ylabel(r'Y [kpc]')

    # colorbar
    cbar = plt.colorbar(hist2d[3], ax=plt.gca())
    cbar.ax.set_ylabel(r'Jeans Mass [M$_{\odot}$]')
    # title
    plt.title(f'Jeans Mass-Weighted Z-Projection of {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)', fontsize=16)
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_jeans_mass_hist.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_jeans_mass_hist.jpg")
    plt.show()


def sfr_hist(snap,
             disk_radius = None,
             disk_half_height = None,
             bins = 100,
             norm = 'linear',
             vmin = None,
             vmax = None,
             cmap = 'viridis',
             save_path = None):
    
    if '_000' in snap.name:
        raise Exception('Cannot compute SFR for the initial snapshot.')
    if snap.has_type_4 is False:
        raise Exception('Snapshot does not contain star particles (type 4). SFR cannot be computed.')

    # Get previous snapshot
    snap.name = snap.filepath.split('/')[-1][0:-5]
    previous_n = (int(snap.name[-1]) - 1) % 10
    previous_snap_name = snap.name[:-1] + str(previous_n)
    previous_filepath = snap.filepath.replace(snap.name, previous_snap_name)

    previous_snap = Snapshot(previous_filepath)

    # Get star formation rates
    if previous_snap.ntype4 > snap.ntype4:
        warnings.warn("Previous snapshot has more star particles. Only SFR of star particles in current snapshot will be plotted.")
        delta_star_mass = snap.mass_type4.to_value(u.solMass) - previous_snap.mass_type4[snap.ntype4].to_value(u.solMass)
        new_star_mass = -1 * previous_snap.mass_type4[snap.ntype4:].to_value(u.solMass)
    else:
        # TODO: check assumption that star particles retain the same index from one snapshot to the next
        delta_star_mass = snap.mass_type4[0:previous_snap.ntype4].to_value(u.solMass) - previous_snap.mass_type4.to_value(u.solMass)
        new_star_mass = snap.mass_type4[previous_snap.ntype4:].to_value(u.solMass)
    delta_star_mass = np.concatenate((delta_star_mass, new_star_mass))
    delta_time = (snap.time - previous_snap.time).to_value(u.yr)
    sfr = delta_star_mass/delta_time

    # Get star particle coordinates
    xs, ys, zs = snap.get_centered_coordinates(parttype=4).T.to_value(u.kpc)
    if (disk_radius is not None):
        range = disk_radius.to_value(u.kpc)
    else:
        range = (snap.boxsize * snap.arepo_length).to_value(u.kpc) / 2
    
    # Get sfr histogram
    hist, xedges, yedges = np.histogram2d(xs, ys, bins=bins, range=[[-range, range], [-range, range]], weights=sfr[0:snap.ntype4], density=False)

    # set color limits
    if vmin is None:
        vmin = np.min(sfr)
    if vmax is None:
        vmax = np.max(sfr)

    if norm == 'linear':
        norm = Normalize(vmin=vmin, vmax=vmax)

    elif norm == 'log':
        if vmin <= 0:
            warnings.warn("LogNorm requires vmin to be positive. Setting vmin to 1e-10.", UserWarning)
            vmin = 1e-10  # avoid log(0) issues        
        if vmax <= 0:
            raise ValueError("LogNorm requires vmin and vmax to be positive.")
        norm = LogNorm(vmin=vmin, vmax=vmax)

    # # Create an alpha mask: 0 where data < threshold, 1 otherwise
    # alpha = np.where(hist.T >= 0.0, 1, 0.0)

    # Plot SFR histogram
    fig, ax = plt.subplots(figsize=(10, 8))
    im = plt.imshow(hist.T, interpolation='nearest', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cmap, norm=norm)#, alpha=alpha)
    cbar = plt.colorbar(im, ax=ax)
    cbar.ax.set_ylabel(r'SFR [M$_{\odot}$ yr$^{-1}$]')
    plt.title(f'SFRs of {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)')
    plt.xlabel(r'X [kpc]')
    plt.ylabel(r'Y [kpc]')
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_sfr_hist.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_sfr_hist.jpg")
    plt.show()

def longitude_velocity(snap,
                        observer_coords_in_kpc = [8, 0, 0],
                        disk_radius = None,
                        disk_half_height = None,
                        bins = 500,
                        norm = 'linear',
                        vmin = None,
                        vmax = None,
                        cmap = 'viridis',
                        save_path = None):
    
    if (disk_radius is not None) & (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    coordinates = snap.get_centered_coordinates().to_value(u.kpc)
    velocities = snap.velocity.to_value(u.km/u.s)
    obs = np.asarray(observer_coords_in_kpc)

    lons = compute_galactic_longitudes(coordinates, obs)
    los_velocities = compute_los_velocity(coordinates, velocities, obs)

    if norm == 'log':
        norm = LogNorm()

    # Plot resolutiom
    fig, ax = plt.subplots(figsize=(10, 8))
    if (disk_radius is not None) & (disk_half_height is not None):
        hist2d = plt.hist2d(lons[snap.disk], los_velocities[snap.disk], cmap=cmap, bins=bins, norm=norm)
    else:
        hist2d = plt.hist2d(lons, los_velocities, cmap=cmap, bins=bins, norm=norm)

    # Add colorbar
    cbar = plt.colorbar(hist2d[3], ax=ax)
    cbar.ax.set_ylabel('Counts')
    plt.title(f'Longitude vs Velocity of {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)')
    plt.xlabel(r'Galactic Longitude [deg]')
    plt.ylabel(r'LoS Velocity [km/s]')
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_LV.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_LV.jpg")
    plt.show()