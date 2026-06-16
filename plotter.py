# -*- coding: utf-8 -*-
"""
Created on Fri May 30 11:43:00 2025

@author: zoefaes
"""

# Imports
import numpy as np
import astropy.units as u
from astropy import constants
from scipy.stats import binned_statistic_2d
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, LinearSegmentedColormap
import matplotlib.animation as animation
import h5py
import warnings
from ArepoISM_tools.snapshot_tools import *
from matplotlib import rcParams
from pathlib import Path
import re

rcParams['font.size'] = 20
# plt.style.use('dark_background')

def truncate_colormap(cmap_name, minval=0.0, maxval=1.0, n=100):
    cmap = plt.get_cmap(cmap_name)
    new_cmap = LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def scale_array(arr, min_val, max_val):
    arr_min = np.min(arr)
    arr_max = np.max(arr)
    return min_val + (arr - arr_min) * (max_val - min_val) / (arr_max - arr_min)

def get_previous_snapshot(snap, n=1):
    """
    Given a snap object with a `.filepath` attribute,
    return the previous snapshot name and filepath.

    Example:
        snapshot600.hdf5 -> snapshot599.hdf5
    """

    path = Path(snap.filepath)

    # Extract trailing number
    match = re.search(r'(\d+)$', snap.name)
    if not match:
        raise ValueError(f"No numeric suffix found in snapshot name: {snap.name}")

    number_str = match.group(1)
    number = int(number_str)
    previous_number = number - n

    if (previous_number) < 0:
        raise ValueError(f"Cannot get {n}th previous snapshot for snapshot: {snap.name}")

    # Preserve zero padding
    previous_number_str = str(previous_number).zfill(3)

    # Rebuild snapshot name
    previous_snap_name = (snap.name[:match.start(1)] + previous_number_str)

    # Rebuild filepath
    previous_filepath = str(path.with_name(previous_snap_name + path.suffix))

    return previous_snap_name, previous_filepath


def interpolate_to_grid(x, y, values, grid_size=(100, 100), method='linear'):
    """
    Interpolate scattered data onto a regular grid.

    Parameters:
    -----------
    x : np.ndarray
        Array of x coordinates of the data points.
    y : np.ndarray
        Array of y coordinates of the data points.
    values : np.ndarray
        Array of shape (N,) with values at the data points.
    grid_size : tuple
        Number of grid points in x and y directions (nx, ny).
    method : str
        Interpolation method: 'nearest', 'linear', or 'cubic'.

    Returns:
    --------
    grid_x, grid_y : np.ndarray
        2D arrays with the x and y coordinates of the grid.
    grid_values : np.ndarray
        2D array with interpolated values on the grid.
    """

    # Define grid
    xi = np.linspace(x.min(), x.max(), grid_size[0])
    yi = np.linspace(y.min(), y.max(), grid_size[1])
    grid_x, grid_y = np.meshgrid(xi, yi)

    # Interpolate values onto grid
    grid_values = griddata(np.column_stack((x, y)), values, (grid_x, grid_y), method=method)

    return grid_x, grid_y, grid_values


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


def compute_epicyclic_frequency_binned(coordinates, velocities, nbins=50, rmin=None, rmax=None, smoothing=0):
    """
    Compute epicyclic frequency κ(R) in radial bins using binned tangential velocity.
    
    Parameters:
        positions: (N, 3) array of particle positions [x, y, z]
        velocities: (N, 3) array of particle velocities [vx, vy, vz]
        nbins: number of radial bins
        rmin, rmax: min and max radius to bin (auto if None)
        smoothing: optional Gaussian smoothing sigma in bins (0 = no smoothing)
        
    Returns:
        R_bin_centers: array of bin centers
        kappa: epicyclic frequency at bin centers
    """
    x, y, z = coordinates
    vx, vy, vz = velocities
    
    # Radial distance
    R = np.sqrt(x**2 + y**2)
    
    # Tangential velocity in cylindrical coordinates
    v_phi = (x * vy - y * vx) / (R + 1e-10)
    
    # Define radial bins
    if rmin is None:
        rmin = R.min()
    if rmax is None:
        rmax = R.max()
        
    bins = np.linspace(rmin, rmax, nbins + 1)
    R_bin_centers = 0.5 * (bins[:-1] + bins[1:])
    
    # Bin statistics
    vphi_mean = np.zeros(nbins)
    counts = np.zeros(nbins)
    
    # Digitize into bins
    inds = np.digitize(R, bins) - 1
    valid = (inds >= 0) & (inds < nbins)
    
    for i in range(nbins):
        in_bin = inds == i
        if np.any(in_bin):
            vphi_mean[i] = np.mean(v_phi[in_bin])
            counts[i] = np.sum(in_bin)
    
    # Angular velocity Omega = v_phi / R
    Omega = vphi_mean / (R_bin_centers + 1e-10)
    
    # Smooth if requested
    if smoothing > 0:
        from scipy.ndimage import gaussian_filter1d
        Omega = gaussian_filter1d(Omega, smoothing)
    
    # Compute derivative of Omega^2
    Omega2 = Omega**2
    dOmega2_dR = np.gradient(Omega2, R_bin_centers)
    
    # Epicyclic frequency
    kappa2 = R_bin_centers * dOmega2_dR + 4 * Omega2
    kappa = np.sqrt(np.maximum(kappa2, 0))  # avoid complex numbers
    
    return R_bin_centers, kappa




##################################################################################################
#                                                                                                #
#                                             PLOTS                                              #
#                                                                                                #
##################################################################################################

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
    if ((disk_radius is not None) | (disk_half_height is not None)) & (axis != 'z'):
        raise ValueError('Disk mask can only be applied to z-axis projection.')
    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)
    
    # set color limits
    if vmin is None:
        vmin = (np.min(snap.mass) << snap.arepo_mass).to_value(u.solMass)
    if vmax is None:
        vmax = (np.max(snap.mass) << snap.arepo_mass).to_value(u.solMass)

    if norm == 'log':
        if vmin <= 0:
            warnings.warn("LogNorm requires vmin to be positive. Setting vmin to 1e-10.", UserWarning)
            vmin = 1e-10  # avoid log(0) issues        
        if vmax <= 0:
            raise ValueError("LogNorm requires vmin and vmax to be positive.")
        norm = LogNorm(vmin=vmin, vmax=vmax)

    # Get coordinates
    x, y, z = snap.get_centered_coordinates().T

    if show_stars:
        if snap.has_type_4:
            xs, ys, zs = snap.get_centered_coordinates(part_type=4).T * snap.arepo_length.to(u.kpc)
        else:
            warnings.warn(f"Snapshot {snap.name} does not have star particles.")

    # Start figure
    fig = plt.figure(figsize=(10,8))
    
    if axis == 'z':
        if (disk_radius is not None) | (disk_half_height is not None): # apply disk mask if applicable
            x = (x[snap.disk] << snap.arepo_length).to_value(u.kpc)
            y = (y[snap.disk] << snap.arepo_length).to_value(u.kpc)
            w = (snap.mass[snap.disk] << snap.arepo_mass).to_value(u.solMass)
            hist2d = plt.hist2d(x, y, bins=bins, norm=norm, weights=w, cmap=cmap)#, cmin=10)
        else: # no disk mask
            x = (x << snap.arepo_length).to_value(u.kpc)
            y = (y << snap.arepo_length).to_value(u.kpc)
            w = (snap.mass << snap.arepo_mass).to_value(u.solMass)
            hist2d = plt.hist2d(x, y, bins=bins, norm=norm, weights=w, cmap=cmap)#, cmin=10)
        if show_stars and snap.has_type_4:
            plt.scatter(xs, ys, s=3, c='r', marker='o')
        plt.xlabel(r'X [kpc]')
        plt.ylabel(r'Y [kpc]')
    elif axis == 'x':
        y = (y << snap.arepo_length).to_value(u.kpc)
        z = (z << snap.arepo_length).to_value(u.kpc)
        w = (snap.mass << snap.arepo_mass).to_value(u.solMass)
        hist2d = plt.hist2d(y, z, bins=bins, norm=norm, weights=w, cmap=cmap)#, cmin=10)
        if show_stars and snap.has_type_4:
            ys = (ys << snap.arepo_length).to_value(u.kpc)
            zs = (zs << snap.arepo_length).to_value(u.kpc)
            plt.scatter(ys, zs, s=3, c='r', marker='o')
        plt.xlabel(r'Y [kpc]')
        plt.ylabel(r'Z [kpc]')
    elif axis == 'y':
        x = (x << snap.arepo_length).to_value(u.kpc)
        z = (z << snap.arepo_length).to_value(u.kpc)
        w = (snap.mass << snap.arepo_mass).to_value(u.solMass)
        hist2d = plt.hist2d(x, z, bins=bins, norm=norm, weights=w, cmap=cmap)
        if show_stars and snap.has_type_4:
            xs = (xs << snap.arepo_length).to_value(u.kpc)
            zs = (zs << snap.arepo_length).to_value(u.kpc)
            plt.scatter(xs, zs, s=3, c='r', marker='o')
        plt.xlabel(r'X [kpc]')
        plt.ylabel(r'Z [kpc]')
    else:
        raise ValueError('Invalid axis. Choose from x, y, or z.')
       
    # colorbar
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'Mass [M$_{\odot}$]')
    # title
    plt.title(f'Mass-Weighted {str(axis).upper()}-Projection of {snap.name} ({snap.time:.2f})', fontsize=16)
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_{axis}_mass_hist.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_{axis}_mass_hist.jpg")
    plt.show()

def mass_3D(snap,
            cmap = 'viridis_r',
            s = 10,
            save_path = None):

    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(projection='3d')

    x, y, z = (snap.get_centered_coordinates().T << snap.arepo_length).to_value(u.kpc)

    if len(x) > 10000:
        warnings.warn('Number of particles is large. Plot may take a long time to render.', UserWarning)

    sizes = scale_array(snap.mass, 10, 30)
    cell_mass = ax.scatter(x, y, z, c=(snap.mass << snap.arepo_mass).to_value(u.solMass), cmap=cmap, s=sizes)
    cbar = plt.colorbar(cell_mass, ax=ax)
    cbar.ax.set_ylabel(r'Mass [M$_{\odot}$]')
    cbar.ax.tick_params(labelsize=14)
    # plt.xticks(fontsize=14)
    # plt.yticks(fontsize=14)
    plt.title(f'Mass of cells for {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)')
    ax.set_xlabel(r'X [kpc]')
    ax.set_ylabel(r'Y [kpc]')
    ax.set_zlabel(r'Z [kpc]')
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_mass_3D.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_mass_3D.jpg")
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
                        species = 'all',
                        save_path = None):
    
    if not isinstance(species, str):
        raise ValueError("Species must be a string. Choose from 'all', 'H2', 'H+', or 'CO'.")

    # get disk mask if applicable
    if ((disk_radius is not None) | (disk_half_height is not None)) & (axis != 'z'):
        raise ValueError('Disk mask can only be applied to z-axis projection.')
    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    # Get coordinates
    x, y, z = snap.get_centered_coordinates().T
    if show_stars:
        if snap.has_type_4:
            xs, ys, zs = snap.get_centered_coordinates(part_type=4).T
        else:
            warnings.warn(f"Snapshot {snap.name} does not have star particles.")

    # Get area per bin
    if (axis == 'z') & (disk_radius is not None) & (disk_half_height is not None):
        length = 2 * disk_half_height.to_value(u.cm)
    else:
        length = (snap.boxsize * snap.arepo_length).to_value(u.cm)

    # Calculate mass per area
    weights = (snap.mass << snap.arepo_mass).to_value(u.solMass)

    if species == 'all':
        weights = weights
        species_label = 'Gas'

    elif species == 'H2':
        weights = weights * snap.xH2
        species_label = 'H$_2$'

    elif species in ['H+', 'HII', 'Hp']:
        weights = weights * snap.xHp
        species_label = 'H$^+$'

    elif species in ['H', 'HI']:
        weights = weights * (1.0 - snap.xHp - 2.0*snap.xH2)
        species_label = 'H'

    elif species == 'CO':
        weights = weights * snap.xCO
        species_label = 'CO'

    else:
        raise ValueError("Invalid species. Choose from 'all', 'H2', 'H+', 'H', or 'CO'.")

    # Start figure
    fig, ax = plt.subplots(figsize=(10,8))
    
    if axis == 'z':
        if (disk_radius is not None) | (disk_half_height is not None): # apply disk mask if applicable
            x = (x[snap.disk] << snap.arepo_length).to_value(u.kpc)
            y = (y[snap.disk] << snap.arepo_length).to_value(u.kpc)
            counts, xedges, yedges = np.histogram2d(x, y, bins=bins, weights=weights[snap.disk])
        else: # no disk mask
            x = (x << snap.arepo_length).to_value(u.kpc)
            y = (y << snap.arepo_length).to_value(u.kpc)
            counts, xedges, yedges = np.histogram2d(x, y, bins=bins, weights=weights)

        if show_stars and snap.has_type_4:
            xs = (xs << snap.arepo_length).to_value(u.kpc)
            ys = (ys << snap.arepo_length).to_value(u.kpc)
        plt.xlabel(r'X [kpc]')
        plt.ylabel(r'Y [kpc]')

    elif axis == 'x':
        y = (y << snap.arepo_length).to_value(u.kpc)
        z = (z << snap.arepo_length).to_value(u.kpc)
        counts, xedges, yedges = np.histogram2d(y, z, bins=bins, weights=weights)

        if show_stars and snap.has_type_4:
            xs = (ys << snap.arepo_length).to_value(u.kpc)
            ys = (zs << snap.arepo_length).to_value(u.kpc)
        plt.xlabel(r'Y [kpc]')
        plt.ylabel(r'Z [kpc]')

    elif axis == 'y':
        x = (x << snap.arepo_length).to_value(u.kpc)
        z = (z << snap.arepo_length).to_value(u.kpc)
        counts, xedges, yedges = np.histogram2d(x, z, bins=bins, weights=weights)
        if show_stars and snap.has_type_4:
            xs = (xs << snap.arepo_length).to_value(u.kpc)
            ys = (zs << snap.arepo_length).to_value(u.kpc)
        plt.xlabel(r'X [kpc]')
        plt.ylabel(r'Z [kpc]')
    else:
        raise ValueError('Invalid axis. Choose from x, y, or z.')

    # convert bin size to area
    dx = np.diff(xedges)[0]
    dy = np.diff(yedges)[0]
    pixel_area = dx * dy  # kpc^2

    # mass surface density (solMass/kpc^2)
    counts = counts / pixel_area
    # convert to g/cm^2
    counts = counts * (constants.M_sun.to_value(u.g) / (u.kpc.to(u.cm))**2)

    # convert to number column density
    mu = 1.4 # mean molecular weight
    m_p = 1.6726219e-24  # g

    counts = counts / (mu * m_p)
    
    # Replace zeros with a very small number to show the lowest color
    # Or alternatively, mask them and set 'bad' color in cmap
    masked_counts = np.ma.masked_where(counts == 0, counts)

    # Setup colormap to show 'bad' (masked) values using lowest color
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad(cmap_obj(0.0))

    # set color limits
    if vmin is None:
        vmin = np.min(counts)
    if vmax is None:
        vmax = np.max(counts)

    if norm == 'log':
        if vmin <= 0:
            warnings.warn("LogNorm requires vmin to be positive. Setting vmin to 1e-10.", UserWarning)
            vmin = 1e-10  # avoid log(0) issues        
        if vmax <= 0:
            raise ValueError("LogNorm requires vmin and vmax to be positive.")
        norm = LogNorm(vmin=vmin, vmax=vmax)

    # Plot using pcolormesh
    h = ax.pcolormesh(xedges, yedges, masked_counts.T, norm=norm, cmap=cmap_obj, shading='auto')

    # plot stars
    if show_stars and snap.has_type_4:
        ax.scatter(xs, ys, s=3, c='r', marker='.')
    
    # t = np.linspace(0, 2*np.pi, 100)
    # r = 0.3
    # x0 = 1.5
    # y0 = -1
    # xc = x0 + r*np.cos(t)
    # yc = y0 + r*np.sin(t)

    # plt.plot(xc, yc, c='firebrick', linewidth=2)
       
    # colorbar
    cbar = plt.colorbar(h, ax=ax)
    cbar.ax.set_ylabel(species_label + r' Column Density [cm$^{-2}$]')
    # title
    plt.title(f'{species_label} Column Density {str(axis).upper()}-Projection {snap.name} ({snap.time:.2f})')
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_{axis}_{species}_column_density.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_{axis}_{species}_column_density.jpg")
    plt.show()


def zoom_column_density(snap,
                        axis = 'z',
                        margin = 0.1, # kpc
                        bins = 1000,
                        norm = 'log',
                        vmin = None,
                        vmax = None,
                        cmap = 'viridis',
                        show_stars = False,
                        save_path = None):

    if snap.has_tracer_field is False:
        raise Exception(f"Snapshot {snap.name} does not have a tracer field.")
    
    if isinstance(margin, u.Quantity):
        margin = margin.to_value(u.kpc)
    else:
        print("margin is not an astropy Quantity. Assuming it is in kpc.")

    coords = snap.get_centered_coordinates() * snap.arepo_length.to(u.kpc)
    xz, yz, zz = coords[(snap.tracer_field > 1e-5)].T
    coords = coords.T # don't transpose before tracer field masking
    xzmin = xz.min()
    xzmax = xz.max()
    yzmin = yz.min()
    yzmax = yz.max()
    zzmin = zz.min()
    zzmax = zz.max()

    cube = (coords[0] > (xzmin - margin)) & (coords[0] < (xzmax + margin)) & (coords[1] > (yzmin - margin)) & (coords[1] < (yzmax + margin)) & (coords[2] > (zzmin - margin)) & (coords[2] < (zzmax + margin))
    xcube = coords[0][cube]
    ycube = coords[1][cube]
    zcube = coords[2][cube]

    if show_stars:
        if snap.has_type_4:
            star_coords = snap.get_centered_coordinates(part_type=4).T * snap.arepo_length.to(u.kpc)
        else:
            warnings.warn(f"Snapshot {snap.name} does not have star particles.")
        star_cube = (star_coords[0] > (xzmin - margin)) & (star_coords[0] < (xzmax + margin)) & (star_coords[1] > (yzmin - margin)) & (star_coords[1] < (yzmax + margin)) & (star_coords[2] > (zzmin - margin)) & (star_coords[2] < (zzmax + margin))
        xs = star_coords[0][star_cube]
        ys = star_coords[1][star_cube]
        zs = star_coords[2][star_cube]

    # Get area per bin
    if (axis == 'z'):
        length = (zzmax - zzmin) * u.kpc.to(u.cm)
    elif (axis == 'y'):
        length = (yzmax - yzmin) * u.kpc.to(u.cm)
    elif (axis == 'x'):
        length = (xzmax - xzmin) * u.kpc.to(u.cm)
    else:
        raise Exception("Invalid axis. Choose 'x', 'y' or 'z'.")

    # Calculate mass per area
    ndensity = snap.ndensity[cube] * (1/snap.arepo_length**3).to_value(1/u.cm**3)
    weights = ndensity * length

    weights = snap.mass[cube] * snap.arepo_mass.to(u.solMass)

    fig, ax = plt.subplots(figsize=(13,10))
    
    if axis == 'z':
        counts, xedges, yedges = np.histogram2d(xcube, ycube, bins=bins, weights=weights)
        if show_stars and snap.has_type_4:
            ax.scatter(xs, ys, s=3, c='red', marker='*', zorder=100)
        ax.set_xlabel(r'X [kpc]')
        ax.set_ylabel(r'Y [kpc]')


    elif axis == 'x':
        counts, xedges, yedges = np.histogram2d(ycube, zcube, bins=bins, weights=weights)
        if show_stars and snap.has_type_4:
            ax.scatter(ys, zs, s=3, c='red', marker='*', zorder=100)
        ax.set_xlabel(r'Y [kpc]')
        ax.set_ylabel(r'Z [kpc]')

    elif axis == 'y':
        counts, xedges, yedges = np.histogram2d(xcube, zcube, bins=bins, weights=weights)
        if show_stars and snap.has_type_4:
            ax.scatter(xs, zs, s=3, c='red', marker='*', zorder=100)
        ax.set_xlabel(r'X [kpc]')
        ax.set_ylabel(r'Z [kpc]')
    else:
        raise ValueError('Invalid axis. Choose from x, y, or z.')
    
    
    def convert_to_column_density(counts, xedges, yedges):
        dx = np.diff(xedges)[0]
        dy = np.diff(yedges)[0]
        pixel_area = dx * dy  # kpc^2

        counts = counts / pixel_area
        counts = counts * (constants.M_sun.to_value(u.g) / (u.kpc.to(u.cm))**2)

        mu = 1.4
        m_p = 1.6726219e-24
        return counts / (mu * m_p)

    counts = convert_to_column_density(counts, xedges, yedges)

    # set color limits
    if vmin is None:
        vmin = np.min(counts)
    if vmax is None:
        vmax = np.max(counts)

    if norm == 'log':
        if vmin <= 0:
            warnings.warn("LogNorm requires vmin to be positive. Setting vmin to 1e-10.", UserWarning)
            vmin = 1e-10  # avoid log(0) issues        
        if vmax <= 0:
            raise ValueError("LogNorm requires vmin and vmax to be positive.")
        norm = LogNorm(vmin=vmin, vmax=vmax)


    # Replace zeros with a very small number to show the lowest color
    # Or alternatively, mask them and set 'bad' color in cmap
    masked_counts = np.ma.masked_where(counts == 0, counts)

    # Setup colormap to show 'bad' (masked) values using lowest color
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad(cmap_obj(0.0))

    # Plot using pcolormesh
    h = ax.pcolormesh(xedges, yedges, masked_counts.T, norm=norm, cmap=cmap_obj, shading='auto')
       
    # colorbar
    cbar = plt.colorbar(h, ax=ax)
    cbar.ax.set_ylabel(r'Column Density [cm$^{-2}$]')
    # title
    # plt.title(f'Column Density {str(axis).upper()}-Projection {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)')
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_{axis}_column_density_proj_zoom.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_{axis}_column_density_proj_zoom.jpg")
    plt.show()


def column_density_projection_comparison(zoom_snap,
                                        base_snap,
                                        axis = 'z',
                                        margin = 0.1, # kpc
                                        zoom_bins = 1000,
                                        base_bins = 300,
                                        norm = 'log',
                                        vmin = None,
                                        vmax = None,
                                        cmap = 'viridis',
                                        show_stars = False,
                                        save_path = None):
    
    if zoom_snap.has_tracer_field is False:
        raise Exception(f"Snapshot {zoom_snap.name} does not have a tracer field.")

    if isinstance(margin, u.Quantity):
        margin = margin.to_value(base_snap.arepo_length)
    else:
        print("margin is not an astropy Quantity. Assuming it is in kpc.")
        margin = (margin << u.kpc).to_value(base_snap.arepo_length)
    
    # Get coordinates
    if show_stars:
        if zoom_snap.has_type_4:
            xsz, ysz, zsz = zoom_snap.get_centered_coordinates(part_type=4).T
        else:
            warnings.warn(f"Snapshot {zoom_snap.name} does not have star particles.")
            show_stars = False
        if base_snap.has_type_4:
            xsb, ysb, zsb = base_snap.get_centered_coordinates(part_type=4).T
        else:
            warnings.warn(f"Snapshot {base_snap.name} does not have star particles.")
            show_stars = False

    bcoords = base_snap.get_centered_coordinates().T
    zcoords = zoom_snap.get_centered_coordinates()
    xz, yz, zz = zcoords[(zoom_snap.tracer_field > 1e-5)].T
    zcoords = zcoords.T

    xzmin, xzmax = xz.min(), xz.max()
    yzmin, yzmax = yz.min(), yz.max()
    zzmin, zzmax = zz.min(), zz.max()

    zcube = (zcoords[0] > (xzmin - margin)) & (zcoords[0] < (xzmax + margin)) & \
            (zcoords[1] > (yzmin - margin)) & (zcoords[1] < (yzmax + margin)) & \
            (zcoords[2] > (zzmin - margin)) & (zcoords[2] < (zzmax + margin))

    bcube = (bcoords[0] > (xzmin - margin)) & (bcoords[0] < (xzmax + margin)) & \
            (bcoords[1] > (yzmin - margin)) & (bcoords[1] < (yzmax + margin)) & \
            (bcoords[2] > (zzmin - margin)) & (bcoords[2] < (zzmax + margin))

    xzmin *= zoom_snap.arepo_length.to(u.kpc)
    yzmin *= zoom_snap.arepo_length.to(u.kpc)
    zzmin *= zoom_snap.arepo_length.to(u.kpc)
    xzmax *= zoom_snap.arepo_length.to(u.kpc)
    yzmax *= zoom_snap.arepo_length.to(u.kpc)
    zzmax *= zoom_snap.arepo_length.to(u.kpc)

    # -------------------------
    # MASS weights (FIXED)
    # -------------------------
    zweights = (zoom_snap.mass[zcube] << zoom_snap.arepo_mass).to_value(u.solMass)
    bweights = (base_snap.mass[bcube] << base_snap.arepo_mass).to_value(u.solMass)

    fig, (axz, axb) = plt.subplots(1, 2, figsize=(18,8), layout='constrained')
    
    if axis == 'z':
        xz = (zcoords[0][zcube] << zoom_snap.arepo_length).to_value(u.kpc)
        yz = (zcoords[1][zcube] << zoom_snap.arepo_length).to_value(u.kpc)
        xb = (bcoords[0][bcube] << base_snap.arepo_length).to_value(u.kpc)
        yb = (bcoords[1][bcube] << base_snap.arepo_length).to_value(u.kpc)

        zcounts, zxedges, zyedges = np.histogram2d(xz, yz, bins=zoom_bins, weights=zweights)
        bcounts, bxedges, byedges = np.histogram2d(xb, yb, bins=base_bins, weights=bweights)

        if show_stars:
            xsz = xsz * zoom_snap.arepo_length.to(u.kpc)
            ysz = ysz * zoom_snap.arepo_length.to(u.kpc)
            xsb = xsb * zoom_snap.arepo_length.to(u.kpc)
            ysb = ysb * zoom_snap.arepo_length.to(u.kpc)

        axz.set_xlabel(r'X [kpc]')
        axz.set_ylabel(r'Y [kpc]')
        axb.set_xlabel(r'X [kpc]')
        axz.set_xlim(xzmin, xzmax)
        axz.set_ylim(yzmin, yzmax)
        axb.set_xlim(xzmin, xzmax)
        axb.set_ylim(yzmin, yzmax)

    elif axis == 'x':
        yz = (zcoords[1][zcube] << zoom_snap.arepo_length).to_value(u.kpc)
        zz = (zcoords[2][zcube] << zoom_snap.arepo_length).to_value(u.kpc)
        yb = (bcoords[1][bcube] << base_snap.arepo_length).to_value(u.kpc)
        zb = (bcoords[2][bcube] << base_snap.arepo_length).to_value(u.kpc)

        zcounts, zxedges, zyedges = np.histogram2d(yz, zz, bins=zoom_bins, weights=zweights)
        bcounts, bxedges, byedges = np.histogram2d(yb, zb, bins=base_bins, weights=bweights)

        if show_stars:
            xsz = ysz * zoom_snap.arepo_length.to(u.kpc)
            ysz = zsz * zoom_snap.arepo_length.to(u.kpc)
            xsb = ysb * zoom_snap.arepo_length.to(u.kpc)
            ysb = zsb * zoom_snap.arepo_length.to(u.kpc)

        axz.set_xlabel(r'Y [kpc]')
        axz.set_ylabel(r'Z [kpc]')
        axb.set_xlabel(r'Y [kpc]')
        axz.set_xlim(yzmin, yzmax)
        axz.set_ylim(zzmin, zzmax)
        axb.set_xlim(yzmin, yzmax)
        axb.set_ylim(zzmin, zzmax)

    elif axis == 'y':
        xz = (zcoords[0][zcube] << zoom_snap.arepo_length).to_value(u.kpc)
        zz = (zcoords[2][zcube] << zoom_snap.arepo_length).to_value(u.kpc)
        xb = (bcoords[0][bcube] << base_snap.arepo_length).to_value(u.kpc)
        zb = (bcoords[2][bcube] << base_snap.arepo_length).to_value(u.kpc)

        zcounts, zxedges, zyedges = np.histogram2d(xz, zz, bins=zoom_bins, weights=zweights)
        bcounts, bxedges, byedges = np.histogram2d(xb, zb, bins=base_bins, weights=bweights)

        if show_stars:
            xsz = xsz * zoom_snap.arepo_length.to(u.kpc)
            zsz = zsz * zoom_snap.arepo_length.to(u.kpc)
            xsb = xsb * zoom_snap.arepo_length.to(u.kpc)
            zsb = zsb * zoom_snap.arepo_length.to(u.kpc)

        axz.set_xlabel(r'X [kpc]')
        axz.set_ylabel(r'Z [kpc]')
        axb.set_xlabel(r'X [kpc]')
        axz.set_xlim(xzmin, xzmax)
        axz.set_ylim(zzmin, zzmax)
        axb.set_xlim(xzmin, xzmax)
        axb.set_ylim(zzmin, zzmax)

    else:
        raise ValueError('Invalid axis.')

    # -------------------------
    # convert to column density
    # -------------------------
    def convert_to_column_density(counts, xedges, yedges):
        dx = np.diff(xedges)[0]
        dy = np.diff(yedges)[0]
        pixel_area = dx * dy  # kpc^2

        counts = counts / pixel_area
        counts = counts * (constants.M_sun.to_value(u.g) / (u.kpc.to(u.cm))**2)

        mu = 1.4
        m_p = 1.6726219e-24
        return counts / (mu * m_p)

    zcounts = convert_to_column_density(zcounts, zxedges, zyedges)
    bcounts = convert_to_column_density(bcounts, bxedges, byedges)

    # -------------------------
    # normalization (FIXED)
    # -------------------------
    if vmin is None:
        vmin = np.min(zcounts[zcounts > 0])
    if vmax is None:
        vmax = np.max(zcounts)

    if norm == 'log':
        norm = LogNorm(vmin=vmin, vmax=vmax)

    # -------------------------
    # masking + plotting
    # -------------------------
    masked_zcounts = np.ma.masked_where(zcounts == 0, zcounts)
    masked_bcounts = np.ma.masked_where(bcounts == 0, bcounts)

    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad(cmap_obj(0.0))

    zh = axz.pcolormesh(zxedges, zyedges, masked_zcounts.T,
                        norm=norm, cmap=cmap_obj, shading='auto')
    bh = axb.pcolormesh(bxedges, byedges, masked_bcounts.T,
                        norm=norm, cmap=cmap_obj, shading='auto')
    
    if show_stars:
            axz.scatter(xsz, ysz, s=3, c='red', marker='*')
            axb.scatter(xsb, ysb, s=3, c='red', marker='*')

    cbar = fig.colorbar(zh, ax=axb)
    cbar.ax.set_ylabel(r'Column Density [cm$^{-2}$]')

    axz.set_title(f'Zoom ({zoom_snap.time:.2f})')
    axb.set_title(f'Base ({base_snap.time:.2f})')

    if save_path is not None:
        plt.savefig(f'{save_path}{zoom_snap.name}_{axis}_column_density_compare_zoom.jpg',
                    dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{zoom_snap.name}_{axis}_column_density_compare_zoom.jpg")

    plt.show()


def zoom_levels(snap,
                axis='z',
                margin=0.5,  # kpc
                bins=1000,
                norm='log',
                vmin=None,
                vmax=None,
                cmaps=['Greens', 'Reds', 'Blues'],
                show_stars=False,
                save_path=None):

    if snap.has_tracer_field is False:
        raise Exception(f"Snapshot {snap.name} does not have a tracer field.")
    elif snap.zoom_level is None:
        raise Exception(f"Snapshot {snap.name} does not have zoom levels.")
    
    plt.rcParams.update({'font.size': 20})

    # -------------------------
    # coordinates
    # -------------------------
    coords = snap.get_centered_coordinates()
    x0, y0, z0 = coords[(snap.zoom_level == 0)].T
    x1, y1, z1 = coords[(snap.zoom_level == 1)].T
    x2, y2, z2 = coords[(snap.zoom_level == 2)].T

    if show_stars and snap.has_type_4:
        xs, ys, zs = snap.get_centered_coordinates(part_type=4).T

    # -------------------------
    # mass weights (SOLAR MASSES)
    # -------------------------
    mass0 = (snap.mass[snap.zoom_level == 0] << snap.arepo_mass).to_value(u.solMass)
    mass1 = (snap.mass[snap.zoom_level == 1] << snap.arepo_mass).to_value(u.solMass)
    mass2 = (snap.mass[snap.zoom_level == 2] << snap.arepo_mass).to_value(u.solMass)

    # -------------------------
    # projection axes
    # -------------------------
    if axis == 'z':
        x0p, y0p = x0, y0
        x1p, y1p = x1, y1
        x2p, y2p = x2, y2
        xlabel, ylabel = 'X [kpc]', 'Y [kpc]'

    elif axis == 'x':
        x0p, y0p = y0, z0
        x1p, y1p = y1, z1
        x2p, y2p = y2, z2
        xlabel, ylabel = 'Y [kpc]', 'Z [kpc]'

    elif axis == 'y':
        x0p, y0p = x0, z0
        x1p, y1p = x1, z1
        x2p, y2p = x2, z2
        xlabel, ylabel = 'X [kpc]', 'Z [kpc]'
    else:
        raise Exception("Invalid axis.")

    # -------------------------
    # convert to kpc
    # -------------------------
    x0p = (x0p << snap.arepo_length).to_value(u.kpc)
    y0p = (y0p << snap.arepo_length).to_value(u.kpc)
    x1p = (x1p << snap.arepo_length).to_value(u.kpc)
    y1p = (y1p << snap.arepo_length).to_value(u.kpc)
    x2p = (x2p << snap.arepo_length).to_value(u.kpc)
    y2p = (y2p << snap.arepo_length).to_value(u.kpc)

    # -------------------------
    # shared bin edges
    # -------------------------
    xmin = min(x0p.min(), x1p.min(), x2p.min())
    xmax = max(x0p.max(), x1p.max(), x2p.max())
    ymin = min(y0p.min(), y1p.min(), y2p.min())
    ymax = max(y0p.max(), y1p.max(), y2p.max())

    xedges = np.linspace(xmin, xmax, bins + 1)
    yedges = np.linspace(ymin, ymax, bins + 1)

    # -------------------------
    # histograms (mass per pixel)
    # -------------------------
    counts0 = np.histogram2d(x0p, y0p, bins=[xedges, yedges], weights=mass0)[0]
    counts1 = np.histogram2d(x1p, y1p, bins=[xedges, yedges], weights=mass1)[0]
    counts2 = np.histogram2d(x2p, y2p, bins=[xedges, yedges], weights=mass2)[0]

    # -------------------------
    # convert to column density
    # -------------------------
    dx = np.diff(xedges)[0]
    dy = np.diff(yedges)[0]
    pixel_area = dx * dy  # kpc^2

    def to_column_density(counts):
        counts = counts / pixel_area  # solMass / kpc^2
        counts = counts * (constants.M_sun.to_value(u.g) / (u.kpc.to(u.cm))**2)  # g/cm^2
        mu = 1.4
        m_p = 1.6726219e-24
        return counts / (mu * m_p)  # cm^-2

    counts0 = to_column_density(counts0)
    counts1 = to_column_density(counts1)
    counts2 = to_column_density(counts2)

    # -------------------------
    # normalization
    # -------------------------
    if vmin is None:
        vmin = np.min([counts0[counts0 > 0].min(),
                       counts1[counts1 > 0].min(),
                       counts2[counts2 > 0].min()])

    if vmax is None:
        vmax = np.max([counts0.max(), counts1.max(), counts2.max()])

    if norm == 'log':
        norm = LogNorm(vmin=vmin, vmax=vmax)

    # -------------------------
    # plotting
    # -------------------------
    fig, ax = plt.subplots(figsize=(13, 10))

    ax.pcolormesh(xedges, yedges, counts0.T, norm=norm, cmap=cmaps[0], shading='auto')
    ax.pcolormesh(xedges, yedges, counts1.T, norm=norm, cmap=cmaps[1], shading='auto')
    ax.pcolormesh(xedges, yedges, counts2.T, norm=norm, cmap=cmaps[2], shading='auto')

    if show_stars and snap.has_type_4:
        if axis == 'z':
            xs = (xs << snap.arepo_length).to_value(u.kpc)
            ys = (ys << snap.arepo_length).to_value(u.kpc)
            ax.scatter(xs, ys, s=3, c='white', marker='*')

    ax.set_xlim(xedges[0], xedges[-1])
    ax.set_ylim(yedges[0], yedges[-1])

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    cbar = plt.colorbar(ax.collections[0], ax=ax)
    cbar.ax.set_ylabel(r'Column Density [cm$^{-2}$]')

    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_{axis}_zoom_levels.jpg',
                    dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_{axis}_zoom_levels.jpg")

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
    if ((disk_radius is not None) | (disk_half_height is not None)) & (axis != 'z'):
        raise ValueError('Disk mask can only be applied to z-axis projection.')
    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    inf_mask = np.isinf(snap.temperature)
    inf_count = inf_mask.sum()
    if inf_count > 0:
        print(f"Percentage of infinite temperature values: {inf_count / len(snap.temperature) * 100:.2f}% \nReplacing with 1e8 K")
        snap.temperature[inf_mask] = 1e8

    # Get coordinates
    x, y, z = snap.get_centered_coordinates().T

    # Start figure
    fig, ax = plt.subplots(figsize=(10,8))
    
    if axis == 'z':
        if (disk_radius is not None) | (disk_half_height is not None): # apply disk mask if applicable
            x = (x[snap.disk] << snap.arepo_length).to_value(u.kpc)
            y = (y[snap.disk] << snap.arepo_length).to_value(u.kpc)
            counts, xedges, yedges = np.histogram2d(x, y, bins=bins, weights=snap.temperature[snap.disk])
        else: # no disk mask
            x = (x << snap.arepo_length).to_value(u.kpc)
            y = (y << snap.arepo_length).to_value(u.kpc)
            counts, xedges, yedges = np.histogram2d(x, y, bins=bins, weights=snap.temperature)
        plt.xlabel(r'X [kpc]')
        plt.ylabel(r'Y [kpc]')
    elif axis == 'x':
        y = (y << snap.arepo_length).to_value(u.kpc)
        z = (z << snap.arepo_length).to_value(u.kpc)
        counts, xedges, yedges = np.histogram2d(y, z, bins=bins, weights=snap.temperature)
        plt.xlabel(r'Y [kpc]')
        plt.ylabel(r'Z [kpc]')
    elif axis == 'y':
        x = (x << snap.arepo_length).to_value(u.kpc)
        z = (z << snap.arepo_length).to_value(u.kpc)
        counts, xedges, yedges = np.histogram2d(x, z, bins=bins, weights=snap.temperature)
        plt.xlabel(r'X [kpc]')
        plt.ylabel(r'Z [kpc]')
    else:
        raise ValueError('Invalid axis. Choose from x, y, or z.')

    # Replace zeros with a very small number to show the lowest color
    # Or alternatively, mask them and set 'bad' color in cmap
    masked_counts = np.ma.masked_where(counts == 0, counts)

    # Setup colormap to show 'bad' (masked) values using mid color
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad(cmap_obj(0.4))

    # set color limits
    if vmin is None:
        vmin = np.min(snap.temperature)
    if vmax is None:
        vmax = np.max(snap.temperature)

    if norm == 'log':
        if vmin <= 0:
            warnings.warn("LogNorm requires vmin to be positive. Setting vmin to 1e-10.", UserWarning)
            vmin = 1e-10  # avoid log(0) issues        
        if vmax <= 0:
            raise ValueError("LogNorm requires vmin and vmax to be positive.")
        norm = LogNorm(vmin=vmin, vmax=vmax)

    # Plot using pcolormesh
    h = ax.pcolormesh(xedges, yedges, masked_counts.T, norm=norm, cmap=cmap_obj, shading='auto')
       
    # colorbar
    cbar = plt.colorbar(h, ax=ax)
    cbar.ax.set_ylabel(r'Temperature [K]')
    # title
    plt.title(f'Temperature {str(axis).upper()}-Projection of {snap.name} ({snap.time:.2f})', fontsize=16)
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
    
    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    log_r_eff = np.log10((np.abs(snap.effective_cell_radius) << snap.arepo_length).to_value(u.pc))
    log_ndens = np.log10((np.abs(snap.ndensity) << (snap.arepo_length ** (-3))).to_value(u.cm**(-3)) + 1e-10)  # add small value to avoid log(0)

    if norm == 'log':
        norm = LogNorm()

    # Plot resolutiom
    fig, ax = plt.subplots(figsize=(10, 8))
    if (disk_radius is not None) | (disk_half_height is not None):
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
    
    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    log_r_eff = np.log10((snap.effective_cell_radius << snap.arepo_length).to_value(u.pc))
    log_mass = np.log10((snap.mass << snap.arepo_mass).to_value(u.solMass) + 1e-10)  # add small value to avoid log(0)

    if norm == 'log':
        norm = LogNorm()

    # Plot resolutiom
    fig, ax = plt.subplots(figsize=(10, 8))
    if (disk_radius is not None) | (disk_half_height is not None):
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
    if ((disk_radius is not None) | (disk_half_height is not None)) & (axis != 'z'):
        raise ValueError('Disk mask can only be applied to z-axis projection.')
    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    # Get coordinates
    x, y, z = (snap.get_centered_coordinates().T << snap.arepo_length).to_value(u.kpc)

    # Get number density weights
    weights = (snap.ndensity << (snap.arepo_length ** (-3))).to_value(u.cm**(-3))

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
        if (disk_radius is not None) | (disk_half_height is not None): # apply disk mask if applicable
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
    if ((disk_radius is not None) | (disk_half_height is not None)) & (axis != 'z'):
        raise ValueError('Disk mask can only be applied to z-axis projection.')
    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    # Get coordinates
    x, y, z = (snap.get_centered_coordinates().T << snap.arepo_length).to_value(u.kpc)

    # Get potential weights
    weights = (np.abs(snap.potential) << (snap.arepo_length/snap.arepo_time)**2).to_value((u.pc/u.yr)**2)

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
        if (disk_radius is not None) | (disk_half_height is not None): # apply disk mask if applicable
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

    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)
                          
    # Get coordinates
    x, y, z = (snap.get_centered_coordinates().T << snap.arepo_length).to_value(u.kpc)
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    vx = (snap.velocity[:, 0] << snap.arepo_velocity).to_value(u.km/u.s)
    vy = (snap.velocity[:, 1] << snap.arepo_velocity).to_value(u.km/u.s)
    vth = vy * np.cos(theta) - vx * np.sin(theta) # Azimuthal velocity

    if norm == 'log':
        norm = LogNorm()

    fig, ax = plt.subplots(figsize=(10, 8))
    if (disk_radius is not None) | (disk_half_height is not None):
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


def mag_field_profile(snap,
                          component = None,
                          disk_radius = None,
                          disk_half_height = None,
                          bins = 500,
                          norm = 'log',
                          cmap = 'viridis',
                          save_path = None):
    
    # Check if magnetic field is available
    if not snap.has_mag_field:
        raise Exception('Snapshot does not have magnetic field data.')
    
    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)
                          
    # Get coordinates
    x, y, z = (snap.get_centered_coordinates().T << snap.arepo_length).to_value(u.kpc)
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    bx, by, bz = snap.mag_field.T * 1/(snap.arepo_mag.cgs.scale * 10**9) # convert to nG
    if component is None:
        b = np.sqrt(bx**2 + by**2 + bz**2) # |B|
        comp = '|B|'
    elif component.lower().strip() in ['r', 'br', 'b_r', 'radial']:
        b = bx*np.cos(theta) + by*np.sin(theta)
        comp = 'B$_r$'
    elif component.lower().strip() in ['theta', 't', 'btheta', 'b_theta', 'azimuthal']:
        b = -bx*np.sin(theta) + by*np.cos(theta)
        comp = 'B$_{\theta}$'
    elif component.lower().strip() in ['z', 'bz', 'b_z', 'vertical']:
        b = bz
        comp = 'B$_z$'
    else:
        raise ValueError('Invalid component. Choose from "r", "theta", "z", or leave as None for |B|.')
    
    if norm == 'log':
        norm = LogNorm()

    fig, ax = plt.subplots(figsize=(10, 8))
    if (disk_radius is not None) | (disk_half_height is not None):
        hist2d = plt.hist2d(r[snap.disk], b[snap.disk], cmap=cmap, bins=bins, norm=norm)# norm=mcolors.LogNorm())
    else:
        hist2d = plt.hist2d(r, b, cmap=cmap, bins=bins, norm=norm)
    # Add colorbar
    cbar = plt.colorbar(hist2d[3], ax=ax)
    cbar.ax.set_ylabel('Counts')
    plt.title(f'{comp} Profile of {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)')
    plt.xlabel(r'Radius [kpc]')
    plt.ylabel(comp + r' [nG]')
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_{component}_mag_profile.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_{component}_mag_profile.jpg")
    plt.show()


def ndensity_profile(snap,
                          disk_radius = None,
                          disk_half_height = None,
                          bins = 500,
                          norm = 'log',
                          cmap = 'viridis',
                          save_path = None):
    
    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)
                          
    # Get coordinates
    x, y, z = (snap.get_centered_coordinates().T << snap.arepo_length).to_value(u.kpc)
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)

    if norm == 'log':
        norm = LogNorm()

    ndensity = (snap.ndensity << (snap.arepo_length**(-3))).to_value(1/u.cm**3)
    density = (snap.density << (snap.arepo_mass / snap.arepo_length**3)).to_value(u.g/u.cm**3)

    fig, ax = plt.subplots(figsize=(10, 8))
    if (disk_radius is not None) | (disk_half_height is not None):
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
                          
    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    # Get coordinates
    x, y, z = (snap.get_centered_coordinates().T << snap.arepo_length).to_value(u.kpc)
    r = np.sqrt(x**2 + y**2)

    if norm == 'log':
        norm = LogNorm()

    log_temperature = np.log10(snap.temperature + 1e-10)  # add small value to avoid log(0)

    inf_count = 0
    for temp in log_temperature:
        if np.isinf(temp):
            inf_count += 1

    if inf_count > 0:
        print(f"Percentage of infinite temperature values: {inf_count / len(snap.temperature) * 100:.2f}% \nReplacing inf values with 1e8 K")
        log_temperature[np.isinf(log_temperature)] = 8

    fig, ax = plt.subplots(figsize=(10, 8))
    if (disk_radius is not None) | (disk_half_height is not None):
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
                          
    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    # Get coordinates
    x, y, z = (snap.get_centered_coordinates().T << snap.arepo_length).to_value(u.kpc)
    r = np.sqrt(x**2 + y**2)

    if norm == 'log':
        norm = LogNorm()

    log_mass = np.log10((snap.mass << snap.arepo_mass).to_value(u.solMass))

    fig, ax = plt.subplots(figsize=(10, 8))
    if (disk_radius is not None) | (disk_half_height is not None):
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


def gas_toomre_profile(snap,
                        disk_radius = None,
                        disk_half_height = None,
                        bins = 100,
                        norm = 'log',
                        cmap = 'viridis',
                        save_path = None):
    
    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)
        x, y, z = (snap.get_centered_coordinates()[snap.disk]).T * snap.arepo_length.to(u.kpc)
        vx, vy, vz = (snap.velocity[snap.disk]).T * snap.arepo_velocity.to(u.km/u.s)
    else:
        x, y, z = snap.get_centered_coordinates().T * snap.arepo_length.to(u.kpc)
        vx, vy, vz = snap.velocity.T * snap.arepo_velocity.to(u.km/u.s)
    if disk_radius is None:
        disk_radius = max(np.sqrt(x**2 + y**2)) * u.kpc
    r, kappa = compute_epicyclic_frequency_binned([x,y,z], [vx, vy, vz], nbins=bins, rmin=0, rmax=disk_radius.to_value(u.kpc))

    if norm == 'log':
        norm = LogNorm()

    fig, ax = plt.subplots(figsize=(10, 8))
    if (disk_radius is not None) | (disk_half_height is not None):
        hist2d = plt.hist2d(r, kappa, cmap=cmap, bins=bins, norm=norm)# norm=mcolors.LogNorm())
    else:
        hist2d = plt.hist2d(r, kappa, cmap=cmap, bins=bins, norm=norm)
    # Add colorbar
    cbar = plt.colorbar(hist2d[3], ax=ax)
    cbar.ax.set_ylabel('Counts')
    plt.title(f'Toomre Profile of {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)')
    plt.xlabel(r'Radius [kpc]')
    plt.ylabel(r'Toomre Parameter')
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_toomre_profile.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_toomre_profile.jpg")
    plt.show()



def phase_space(snap,
                    disk_radius = None,
                    disk_half_height = None,
                    bins = 500,
                    weights = 'mass',
                    vmin = None,
                    vmax = None,
                    norm = 'log',
                    cmap = 'viridis',
                    save_path = None):
    
    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    temperature = snap.temperature

    inf_count = 0
    for temp in temperature:
        if np.isinf(temp):
            inf_count += 1
    if inf_count > 0:
        print(f"Percentage of infinite temperature values: {inf_count / len(snap.temperature) * 100:.2f}% \nReplacing with 1e8 K")

    temperature[np.isinf(temperature)] = 1e8

    log_temp = np.log10(temperature + 1e-10)  # add small value to avoid log(0)
    log_ndens = np.log10((snap.ndensity << (snap.arepo_length ** (-3))).to_value(u.cm**(-3)) + 1e-10)  # add small value to avoid log(0)

    if norm == 'log':
        norm = LogNorm(vmin=vmin, vmax=vmax)

    weights = (snap.mass << snap.arepo_mass).to_value(u.solMass) if (weights == 'mass') else None

    # Plot phase space
    fig, ax = plt.subplots(figsize=(10, 8))
    if (disk_radius is not None) | (disk_half_height is not None):
        if weights is None:
            hist2d = plt.hist2d(log_ndens[snap.disk], log_temp[snap.disk], weights=weights, cmap=cmap, bins=bins, norm=norm)# norm=mcolors.LogNorm())
        else:
            hist2d = plt.hist2d(log_ndens[snap.disk], log_temp[snap.disk], weights=weights[snap.disk], cmap=cmap, bins=bins, norm=norm)# norm=mcolors.LogNorm())
    else:
        hist2d = plt.hist2d(log_ndens, log_temp, weights=weights, cmap=cmap, bins=bins, norm=norm)# norm=mcolors.LogNorm())
    # Add colorbar
    cbar = plt.colorbar(hist2d[3], ax=ax)
    if weights is not None:
        cbar.ax.set_ylabel(r'Mass [M$_{\odot}$]')
    else:
        cbar.ax.set_ylabel('Counts')
    plt.title(f'Phase Space of {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)')
    plt.xlabel(r'Log(Number Density) [cm$^{-3}$]')
    plt.ylabel(r'Log(Temperature) [K]')
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_phase_space.jpg', dpi=200, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_phase_space.jpg")
    plt.show()


def mag_ndensity(snap,
                      disk_radius = None,
                      disk_half_height = None,
                      bins = 500,
                      norm = 'log',
                      cmap = 'viridis',
                      save_path = None):
    # Check if magnetic field is available
    if not snap.has_mag_field:
        raise Exception('Snapshot does not have magnetic field data.')
    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    bx, by, bz = snap.mag_field.T * 1/(snap.arepo_mag.cgs.scale * 10**9) # convert to nG
    b = np.sqrt(bx**2 + by**2 + bz**2) # |B|

    log_ndens = np.log10((snap.ndensity << (snap.arepo_length ** (-3))).to_value(u.cm**(-3)) + 1e-10)  # add small value to avoid log(0)

    if norm == 'log':
        norm = LogNorm()

    # Plot resolutiom
    fig, ax = plt.subplots(figsize=(10, 8))
    if (disk_radius is not None) | (disk_half_height is not None):
        hist2d = plt.hist2d(log_ndens[snap.disk], b[snap.disk], cmap=cmap, bins=bins, norm=norm)# norm=mcolors.LogNorm())
    else:
        hist2d = plt.hist2d(log_ndens, b, cmap=cmap, bins=bins, norm=norm)# norm=mcolors.LogNorm())
    # Add colorbar
    cbar = plt.colorbar(hist2d[3], ax=ax)
    cbar.ax.set_ylabel('Counts')
    plt.title(f'Magnetic Field Strength vs Number Density for {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)')
    plt.xlabel(r'Log(Number Density) [cm$^{-3}$]')
    plt.ylabel(r'|B| [$\mu$G]')
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_B_ndensity.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_B_ndensity.jpg")
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
    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)
    
    # Get coordinates
    x, y, z = (snap.get_centered_coordinates() << snap.arepo_length).T.to_value(u.kpc)
    w = 5 * (snap.internal_energy << (snap.arepo_energy/snap.arepo_mass)).to_value(u.J/u.kg) * (gamma - 1) / constants.G.si.value
    # internal energy [m**2/s**2], G [m**3/(kg s**2)], therefore w [kg/m]
    conversion_factor = (1*u.kg/u.m).to_value(u.solMass/u.kpc)
    w = w * conversion_factor  # convert to M_sun/kpc
    v = ((3/(4 * np.pi * snap.density)) << (1/snap.arepo_density)).to_value(u.kpc**3/u.solMass)
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
    
    if (disk_radius is not None) | (disk_half_height is not None): # apply disk mask if applicable
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
             dt = 1,
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
    previous_snap_name, previous_filepath = get_previous_snapshot(snap, n=dt)

    previous_snap = Snapshot(previous_filepath)

    # Get star formation rates
    if previous_snap.n_type_4 > snap.n_type_4:
        warnings.warn("Previous snapshot has more star particles. Only SFR of star particles in current snapshot will be plotted.")
        delta_star_mass = (snap.mass_type_4 - previous_snap.mass_type_4[snap.n_type_4])  * snap.arepo_mass.to(u.solMass)
        new_star_mass = -previous_snap.mass_type_4[snap.n_type_4:]  * snap.arepo_mass.to(u.solMass)
    else:
        # TODO: check assumption that star particles retain the same index from one snapshot to the next
        delta_star_mass = ((snap.mass_type_4[0:previous_snap.n_type_4] - previous_snap.mass_type_4) << snap.arepo_mass).to_value(u.solMass)
        new_star_mass = snap.mass_type_4[previous_snap.n_type_4:] * snap.arepo_mass.to(u.solMass)
    delta_star_mass = np.concatenate((delta_star_mass, new_star_mass))
    delta_time = (snap.time - previous_snap.time).to_value(u.yr)
    sfr = delta_star_mass/delta_time

    # Get star particle coordinates
    xs, ys, zs = snap.get_centered_coordinates(part_type=4).T * snap.arepo_length.to(u.kpc)
    if (disk_radius is not None):
        range = disk_radius.to_value(u.kpc)
    else:
        range = snap.boxsize * snap.arepo_length.to(u.kpc) / 2

    # Get sfr histogram
    hist, xedges, yedges = np.histogram2d(xs, ys, bins=bins, range=[[-range, range], [-range, range]], weights=sfr[0:snap.n_type_4], density=False)

    # Ignore negative SFRs (star particle mass loss due to SN feedback)
    hist[hist < 0] = 0

    # set color limits
    if vmin is None:
        vmin = np.min(hist)
    if vmax is None:
        vmax = np.max(hist)

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
    plt.title(f'SFRs of {snap.name} ({snap.time:.2f})')
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
                        mass_weighted = True,
                        vmin = None,
                        vmax = None,
                        cmap = 'viridis',
                        save_path = None):
    
    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_disk(radius = disk_radius, half_height = disk_half_height)

    coordinates = (snap.get_centered_coordinates() << snap.arepo_length).to_value(u.kpc)
    velocities = (snap.velocity << snap.arepo_velocity).to_value(u.km/u.s)
    obs = np.asarray(observer_coords_in_kpc)

    lons = compute_galactic_longitudes(coordinates, obs)
    los_velocities = compute_los_velocity(coordinates, velocities, obs)

    if mass_weighted is True:
        weights = snap.mass * snap.arepo_mass.to(u.Msun)
    else:
        weights = None

    # Plot resolutiom
    fig, ax = plt.subplots(figsize=(10, 8))
    if (disk_radius is not None) | (disk_half_height is not None):
        weights = weights[snap.disk] if mass_weighted is True else None
        hist, xedges, yedges = np.histogram2d(lons[snap.disk], los_velocities[snap.disk], bins=bins, weights=weights)
    else:
        hist, xedges, yedges = np.histogram2d(lons, los_velocities, bins=bins, weights=weights)

    masked_hist = np.ma.masked_where(hist == 0, hist)

    # Setup colormap to show 'bad' (masked) values using lowest color
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad(cmap_obj(0.0))

    # set color limits
    if vmin is None:
        vmin = np.min(hist)
    if vmax is None:
        vmax = np.max(hist)

    if norm == 'log':
        if vmin <= 0:
            warnings.warn("LogNorm requires vmin to be positive. Setting vmin to 1e-10.", UserWarning)
            vmin = 1e-10  # avoid log(0) issues        
        if vmax <= 0:
            raise ValueError("LogNorm requires vmin and vmax to be positive.")
        norm = LogNorm(vmin=vmin, vmax=vmax)

    # Plot using pcolormesh
    h = ax.pcolormesh(xedges, yedges, masked_hist.T, norm=norm, cmap=cmap_obj, shading='auto')

    # Add colorbar
    cbar_label = r'Mass [M$_\odot$]' if mass_weighted is True else 'Counts'
    cbar = plt.colorbar(h, ax=ax)
    cbar.ax.set_ylabel(cbar_label)
    plt.title(f'Longitude vs Velocity of {snap.model} ({snap.time:.2f})')
    plt.xlabel(r'Galactic Longitude [deg]')
    plt.ylabel(r'Line of Sight Velocity [km/s]')
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_LV.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_LV.jpg")
    plt.show()

def stellar_mass_evolution(filepaths):
    snap = Snapshot(filepaths[0])
    times = []
    stellar_masses = []
    for filepath in filepaths:
        with h5py.File(filepath, 'r') as file:
            stellar_masses.append(np.sum(file['part_type4']['Masses']) * snap.arepo_mass.to(u.Msun))
            times.append((file['Header'].attrs['Time']) * snap.arepo_time.to(u.Myr))
    
    plt.figure(figsize=(10, 6))
    plt.plot(times, stellar_masses)
    plt.xlabel('Time [Myr]')
    plt.ylabel(r'Star Particle Mass [M$_{\odot}$]')
    plt.title('Stellar Mass Evolution')
    plt.grid()
    plt.show()

def star_mass_hist(snap,
                    axis='z',
                    disk_radius=None,
                    disk_half_height=None,
                    bins=500,
                    norm='log',
                    vmin=None,
                    vmax=None,
                    cmap='viridis',
                    show_tracer_field=False,
                    ax=None,
                    save_path=None):

    if ((disk_radius is not None) | (disk_half_height is not None)) & (axis != 'z'):
        raise ValueError('Disk mask can only be applied to z-axis projection.')
    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_disk(radius=disk_radius, half_height=disk_half_height, part_type=4)

    # Get coordinates
    if snap.has_type_4:
        xs, ys, zs = snap.get_centered_coordinates(part_type=4).T * snap.arepo_length.to(u.kpc)
    else:
        raise Exception(f"Snapshot {snap.name} does not have star particles.")

    # Calculate mass per area
    weights = snap.mass_type_4 * snap.arepo_mass.to(u.solMass)

    # Coordinate projection
    if axis == 'z':
        coords = (xs[snap.disk], ys[snap.disk]) if (disk_radius is not None) else (xs, ys)
        weights = weights[snap.disk] if (disk_radius is not None) else weights
        xlabel, ylabel = 'X [kpc]', 'Y [kpc]'
    elif axis == 'x':
        coords = (ys, zs)
        xlabel, ylabel = 'Y [kpc]', 'Z [kpc]'
    elif axis == 'y':
        coords = (xs, zs)
        xlabel, ylabel = 'X [kpc]', 'Z [kpc]'
    else:
        raise ValueError('Invalid axis. Choose from x, y, or z.')

    # Histogram range
    if disk_radius is not None:
        x_extent = [-disk_radius.to_value(u.kpc), disk_radius.to_value(u.kpc)]
        y_extent = [-disk_radius.to_value(u.kpc), disk_radius.to_value(u.kpc)]
    else:
        x_extent = [np.min(coords[0]), np.max(coords[0])]
        y_extent = [np.min(coords[1]), np.max(coords[1])]
    extent = x_extent + y_extent

    # Compute histogram
    H, xedges, yedges = np.histogram2d(coords[0], coords[1], bins=bins, weights=weights,
                                       range=[[x_extent[0], x_extent[1]], [y_extent[0], y_extent[1]]])

    # Fill empty bins with vmin to get uniform background
    # Set vmin/vmax based on non-zero entries
    nonzero = H[H > 0]
    if vmin is None:
        vmin = np.min(nonzero) if nonzero.size > 0 else 1e-10
    if vmax is None:
        vmax = np.max(nonzero) if nonzero.size > 0 else vmin * 10

    if norm == 'log':
        if vmin <= 0:
            warnings.warn("LogNorm requires vmin to be positive. Setting vmin to 1e-10.", UserWarning)
            vmin = 1e-10
        if vmax <= 0:
            raise ValueError("LogNorm requires vmin and vmax to be positive.")
        norm_obj = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm_obj = None

    H_filled = np.nan_to_num(H, nan=vmin)

    # Enforce that all bins are ≥ vmin
    H_filled = np.where(H_filled > 0, H_filled, vmin)

    # Show with imshow
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))

    img = ax.imshow(H_filled.T, origin='lower', cmap=cmap, norm=norm_obj,
                    extent=extent, aspect='auto')

    if show_tracer_field:
        if snap.has_tracer_field:
            xf, yf, zf = (snap.get_centered_coordinates()[snap.tracer_field > 1e-3] * snap.arepo_length.to(u.kpc)).T
            vf = snap.tracer_field[snap.tracer_field > 1e-3]
            vf_norm = (vf - np.min(vf)) / (np.max(vf) - np.min(vf) + 1e-10)
            if axis == 'z':
                scatter = ax.scatter(xf, yf, s=3, c='r', alpha=vf_norm/5)
            elif axis == 'x':
                scatter = ax.scatter(yf, zf, s=3, c='r', alpha=vf_norm/5)
            elif axis == 'y':
                scatter = ax.scatter(xf, zf, s=3, c='r', alpha=vf_norm/5)
        else:
            warnings.warn(f"Snapshot {snap.name} does not have a tracer field.")
            show_tracer_field = False

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    cbar = plt.colorbar(img, ax=ax)
    cbar.ax.set_ylabel(r'Stellar Mass [M$_{\odot}$]')
    ax.set_title(f'Stellar Mass {axis.upper()}-Projection {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)')

    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_{axis}_stellar_mass.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_{axis}_stellar_mass.jpg")


def velocity_streamlines(snap, 
                         axis='z', 
                         bins=200,
                         interp_method='linear',
                         linedensity=1.5,
                         linewidth=1,
                         disk_radius=None,
                         disk_half_height=None,
                         cmap='viridis',
                         vmin=None,
                         vmax=None,
                         norm=None,
                         save_path=None):
    """Plot velocity streamlines in the snapshot.

    Parameters
    ----------
    snap : Snapshot object
        The snapshot for which to plot velocity streamlines.
    axis : str, optional
        The axis along which to plot the streamlines ('x', 'y', or 'z'), by default 'z'.
    save_path : str, optional
        The path where the plot should be saved, by default None.
    """

    # get disk mask if applicable
    if ((disk_radius is not None) | (disk_half_height is not None)) & (axis != 'z'):
        raise ValueError('Disk mask can only be applied to z-axis projection.')
    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_box(xmax = disk_radius, ymax = disk_radius, zmax = disk_half_height)

    # Get coordinates and velocities
    x, y, z = snap.get_centered_coordinates().T * snap.arepo_length.to(u.kpc)
    vx, vy, vz = snap.velocity.T * snap.arepo_velocity.to(u.km/u.s)
    # speed = np.linalg.norm(snap.velocity, axis=1) * snap.arepo_velocity.to(u.km/u.s)

    if axis == 'z':
        speed = np.sqrt(vx**2 + vy**2)
        if (disk_radius is not None) | (disk_half_height is not None): # apply disk mask if applicable
            X, Y, U = interpolate_to_grid(x[snap.box], y[snap.box], vx[snap.box], grid_size=(bins, bins), method=interp_method)
            X, Y, V = interpolate_to_grid(x[snap.box], y[snap.box], vy[snap.box], grid_size=(bins, bins), method=interp_method)
            speeds, xedges, yedges = np.histogram2d(x[snap.box], y[snap.box], bins=bins, weights=speed[snap.box])
            ncounts = np.histogram2d(x[snap.box], y[snap.box], bins=bins)[0]
        else: # no disk mask
            X, Y, U = interpolate_to_grid(x, y, vx, grid_size=(bins, bins), method=interp_method)
            X, Y, V = interpolate_to_grid(x, y, vy, grid_size=(bins, bins), method=interp_method)
            speeds, xedges, yedges = np.histogram2d(x, y, bins=bins, weights=speed)
            ncounts = np.histogram2d(x, y, bins=bins)[0]
        xlabel, ylabel = 'X [kpc]', 'Y [kpc]'
    elif axis == 'x':
        speed = np.sqrt(vy**2 + vz**2)
        X, Y, U = interpolate_to_grid(y, z, vy, grid_size=(bins, bins), method=interp_method)
        X, Y, V = interpolate_to_grid(y, z, vz, grid_size=(bins, bins), method=interp_method)
        speeds, xedges, yedges = np.histogram2d(y, z, bins=bins, weights=speed)
        ncounts = np.histogram2d(y, z, bins=bins)[0]
        xlabel, ylabel = 'Y [kpc]', 'Z [kpc]'
    elif axis == 'y':
        speed = np.sqrt(vx**2 + vz**2)
        X, Y, U = interpolate_to_grid(x, z, vx, grid_size=(bins, bins), method=interp_method)
        X, Y, V = interpolate_to_grid(x, z, vz, grid_size=(bins, bins), method=interp_method)
        speeds, xedges, yedges = np.histogram2d(x, z, bins=bins, weights=speed)
        ncounts = np.histogram2d(x, z, bins=bins)[0]
        xlabel, ylabel = 'X [kpc]', 'Z [kpc]'
    else:
        raise ValueError('Invalid axis. Choose from x, y, or z.')

    fig, ax = plt.subplots(figsize=(10, 8))

    # mean speed
    counts = speeds / ncounts
    counts[ncounts == 0] = 0 # avoid divide-by-zero

    # Replace zeros with a very small number to show the lowest color
    # Or alternatively, mask them and set 'bad' color in cmap
    masked_counts = np.ma.masked_where(counts == 0, counts)

    # Setup colormap to show 'bad' (masked) values using highest color
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad(cmap_obj(0.8))

    # set color limits
    if vmin is None:
        vmin = np.min(counts)
    if vmax is None:
        vmax = np.max(counts)

    if norm == 'log':
        if vmin <= 0:
            warnings.warn("LogNorm requires vmin to be positive. Setting vmin to 1e-10.", UserWarning)
            vmin = 1e-10  # avoid log(0) issues        
        if vmax <= 0:
            raise ValueError("LogNorm requires vmin and vmax to be positive.")
        norm = LogNorm(vmin=vmin, vmax=vmax)

    # Plot speed histogram
    h = ax.pcolormesh(xedges, yedges, masked_counts.T, norm=norm, cmap=cmap_obj, shading='auto')

    # # Plot streamlines
    # strm = ax.streamplot(X, Y, U, V, color=np.sqrt(U**2 + V**2), cmap='viridis', density=1.5)
    strm = ax.streamplot(X, Y, U, V, color='k', density=linedensity, linewidth=linewidth)

    # colorbar
    cbar = plt.colorbar(h, ax=ax)
    cbar.ax.set_ylabel(r'|v| [km/s]')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if disk_radius is not None:
        ax.set_xlim(-disk_radius.to_value(u.kpc), disk_radius.to_value(u.kpc))
        ax.set_ylim(-disk_radius.to_value(u.kpc), disk_radius.to_value(u.kpc))
    ax.set_title(f'Velocity Streamlines {axis.upper()}-Projection of {snap.name} ({snap.time.to_value(u.Myr):.2f} Myr)')
    
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_{axis}_velocity_streamlines.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_{axis}_velocity_streamlines.jpg")
    plt.show()



def magnetic_streamlines(snap, 
                         axis='z', 
                         bins=200,
                         interp_method='linear',
                         linedensity=1.5,
                         linewidth=1,
                         arrowstyle='->',
                         arrowsize=1,
                         disk_radius=None,
                         disk_half_height=None,
                         cmap='viridis',
                         vmin=None,
                         vmax=None,
                         norm=None,
                         save_path=None):
    """Plot velocity streamlines in the snapshot.

    Parameters
    ----------
    snap : Snapshot object
        The snapshot for which to plot velocity streamlines.
    axis : str, optional
        The axis along which to plot the streamlines ('x', 'y', or 'z'), by default 'z'.
    save_path : str, optional
        The path where the plot should be saved, by default None.
    """

    # get disk mask if applicable
    if ((disk_radius is not None) | (disk_half_height is not None)) & (axis != 'z'):
        raise ValueError('Disk mask can only be applied to z-axis projection.')
    if (disk_radius is not None) | (disk_half_height is not None):
        snap.set_box(xmax = disk_radius, ymax = disk_radius, zmax = disk_half_height)

    # Get coordinates and velocities
    x, y, z = snap.get_centered_coordinates().T * snap.arepo_length.to(u.kpc)
    bx, by, bz = snap.mag_field.T * 1/(snap.arepo_mag.cgs.scale * 10**9) # convert to nG
    bmag = np.sqrt(bx**2 + by**2 + bz**2)
    # speed = np.linalg.norm(snap.velocity, axis=1) * snap.arepo_velocity.to(u.km/u.s)

    if axis == 'z':
        if (disk_radius is not None) | (disk_half_height is not None): # apply disk mask if applicable
            X, Y, U = interpolate_to_grid(x[snap.box], y[snap.box], bx[snap.box], grid_size=(bins, bins), method=interp_method)
            X, Y, V = interpolate_to_grid(x[snap.box], y[snap.box], by[snap.box], grid_size=(bins, bins), method=interp_method)
            bmags, xedges, yedges = np.histogram2d(x[snap.box], y[snap.box], bins=bins, weights=bmag[snap.box])
            ncounts = np.histogram2d(x[snap.box], y[snap.box], bins=bins)[0]
        else: # no disk mask
            X, Y, U = interpolate_to_grid(x, y, bx, grid_size=(bins, bins), method=interp_method)
            X, Y, V = interpolate_to_grid(x, y, by, grid_size=(bins, bins), method=interp_method)
            bmags, xedges, yedges = np.histogram2d(x, y, bins=bins, weights=bmag)
            ncounts = np.histogram2d(x, y, bins=bins)[0]
        xlabel, ylabel = 'X [kpc]', 'Y [kpc]'
    elif axis == 'x':
        bmag = np.sqrt(by**2 + bz**2)
        X, Y, U = interpolate_to_grid(y, z, by, grid_size=(bins, bins), method=interp_method)
        X, Y, V = interpolate_to_grid(y, z, bz, grid_size=(bins, bins), method=interp_method)
        bmags, xedges, yedges = np.histogram2d(y, z, bins=bins, weights=bmag)
        ncounts = np.histogram2d(y, z, bins=bins)[0]
        xlabel, ylabel = 'Y [kpc]', 'Z [kpc]'
    elif axis == 'y':
        bmag = np.sqrt(bx**2 + bz**2)
        X, Y, U = interpolate_to_grid(x, z, bx, grid_size=(bins, bins), method=interp_method)
        X, Y, V = interpolate_to_grid(x, z, bz, grid_size=(bins, bins), method=interp_method)
        bmags, xedges, yedges = np.histogram2d(x, z, bins=bins, weights=bmag)
        ncounts = np.histogram2d(x, z, bins=bins)[0]
        xlabel, ylabel = 'X [kpc]', 'Z [kpc]'
    else:
        raise ValueError('Invalid axis. Choose from x, y, or z.')

    fig, ax = plt.subplots(figsize=(10, 8))

    # mean |B|
    counts = bmags / ncounts
    counts[ncounts == 0] = 0 # avoid divide-by-zero

    # Replace zeros with a very small number to show the lowest color
    # Or alternatively, mask them and set 'bad' color in cmap
    masked_counts = np.ma.masked_where(counts == 0, counts)

    # Setup colormap to show 'bad' (masked) values using lowest color
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad(cmap_obj(0.0))

    # set color limits
    if vmin is None:
        vmin = np.min(counts)
    if vmax is None:
        vmax = np.max(counts)

    if norm == 'log':
        if vmin <= 0:
            warnings.warn("LogNorm requires vmin to be positive. Setting vmin to 1e-10.", UserWarning)
            vmin = 1e-10  # avoid log(0) issues        
        if vmax <= 0:
            raise ValueError("LogNorm requires vmin and vmax to be positive.")
        norm = LogNorm(vmin=vmin, vmax=vmax)

    # Plot speed histogram
    h = ax.pcolormesh(xedges, yedges, masked_counts.T, norm=norm, cmap=cmap_obj, shading='auto')

    # # Plot streamlines
    # strm = ax.streamplot(X, Y, U, V, color=np.sqrt(U**2 + V**2), cmap='viridis', density=1.5)
    strm = ax.streamplot(X, Y, U, V, color='k', density=linedensity, linewidth=linewidth, arrowstyle=arrowstyle, arrowsize=arrowsize)

    # colorbar
    cbar = plt.colorbar(h, ax=ax)
    cbar.ax.set_ylabel(r'|B| [nG]')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if disk_radius is not None:
        ax.set_xlim(-disk_radius.to_value(u.kpc), disk_radius.to_value(u.kpc))
        ax.set_ylim(-disk_radius.to_value(u.kpc), disk_radius.to_value(u.kpc))
    ax.set_title(f'Magnetic Streamlines {axis.upper()}-Projection of {snap.name} ({snap.time:.2f})')
    
    if save_path is not None:
        plt.savefig(f'{save_path}{snap.name}_{axis}_magnetic_streamlines.jpg', dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}{snap.name}_{axis}_magnetic_streamlines.jpg")
    plt.show()


def region(snap, 
           threshold_quantity, 
           threshold_value, 
           threshold_unit, 
           axis = 'z', 
           cmap = 'viridis',
           alpha = 0.5,
           save_path = None):
    """Plot a region of the snapshot based on a threshold quantity.

    This function plots a region of the snapshot based on a threshold quantity, such as temperature, density, or number density.
    The region is defined by the threshold value and unit, and can be plotted in the x-y, y-z, or x-z plane.

    Parameters
    ----------
    snap : Snapshot object
        _description_
    threshold_quantity : str
        _description_
    threshold_value : float
        _description_
    threshold_unit : str | astropy unit
        _description_
    axis : str, optional
        _description_, by default 'z'
    save_path : str, optional
        _description_, by default None
    """    

    try:
        q = getattr(snap, threshold_quantity)
    except AttributeError:
        raise Exception(f"Snapshot does not have a '{threshold_quantity}' attribute. Available attributes: {', '.join(snap.available_attributes)}")
    
    # mask = 