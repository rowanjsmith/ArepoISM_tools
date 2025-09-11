# -*- coding: utf-8 -*-
"""
This is an example of how to use the plotting function with a Snapshot object.

Created on Thu Jun 12 14:48:01 2025

@author: zoefaes
"""

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from snapshot_tools import *
import plotter as plot
import os

# TUTORIAL
# TODO: add a few lines about how to use astropy units

# Specify your snapshot filepath
filepath = '/cosma8/data/dp058/dc-prab1/arepo2_len/snapshots/V1a_245.hdf5'

# Specify the path where you want to save your plots to
save_path = "~/figures/"
save_path = os.path.expanduser(save_path)

# Read in your snapshot
snap = Snapshot(filepath)

# Print snapshot time
print(f"Snapshot time: {snap.time.to_value(u.Myr)} Myr")
print(snap.time.to(snap.arepo_time))

# Specify the extent of the disk in your snapshot
# You can get an idea of the size of your disk by plotting the whole box initially. Try:
# plot.column_density_hist(snap, save_path=save_path)
# THESE QUANTITIES MUST HAVE UNITS OF LENGTH (use astropy quantities as shown here)
radius = 10 * u.kpc
half_height = 4 * u.kpc

# Examples of how to call the plotting functions
plot.mass_hist(snap, disk_radius=radius, disk_half_height=half_height, vmin=5e2, vmax=1e5, save_path=save_path)
plot.mass_hist(snap, save_path=save_path)
plot.mass_hist(snap, axis='x', vmin=1e-6, vmax=1e4, save_path=save_path)
plot.temperature_hist(snap, save_path=save_path)
plot.column_density_hist(snap, show_stars=True, disk_radius=radius, disk_half_height=half_height, vmin=1e20, vmax=1e25, cmap='binary', save_path=save_path)
plot.resolution(snap, save_path=save_path)
plot.mass_resolution(snap, save_path=save_path)
plot.ndensity_hist(snap, save_path=save_path)
plot.potential_hist(snap, save_path=save_path)
plot.velocity_profile(snap, save_path=save_path)
plot.ndensity_profile(snap, save_path=save_path)
plot.temperature_profile(snap, save_path=save_path)
plot.mass_profile(snap, save_path=save_path)
plot.phase_space(snap, save_path=save_path)
plot.jeans_mass_hist(snap, save_path=save_path)
plot.sfr_hist(snap, disk_radius=5*u.kpc, vmin=-1e-5, save_path=save_path)
plot.longitude_velocity(snap, save_path=save_path)

# Most functions will accept 'disk_radius' and 'disk_half_height' parameters to restrict your domain to the galaxy disk
# Set limits on the value ranges with 'vmin' & 'vmax'
# Change the colormap by specifying a matplotlib colormap with 'cmap'
# As the snapshots typically have a large number of gas particles, 2D histograms are used to map quantities instead of 
# e.g. a scatter plot. You can vary the resolution of your plot by changing the number of bins with 'bins'