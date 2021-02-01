#!/home/labmember/anaconda2/envs/py36/bin/python

# Author: Dhananjay Bhaskar <dhananjay_bhaskar@brown.edu>
# Sample local density and velocity uniformly on sphere and visualize

import os
import math
import multiprocessing

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

from matplotlib import cm, colors
from mpl_toolkits.mplot3d import axes3d, Axes3D
from matplotlib.colors import LightSource

from sklearn.neighbors import KDTree
from sklearn.neighbors import KernelDensity

from joblib import Parallel, delayed

### load data

folder = "SCM_sim_test_4"
in_file = "SCM_sim_data_4.npy"

savedict = np.load(folder + os.sep + in_file, allow_pickle=True)
savedata = savedict.item()

sph_rad = savedata['r']
pos_dat = savedata['pos_mat']
dir_dat = savedata['dir_mat']
vel_dat = savedata['vel_mat']

(num_agents, n_iter) = np.shape(vel_dat)

v_max = np.max(vel_dat[:])
v_min = np.min(vel_dat[:])
num_cores = multiprocessing.cpu_count()

### helper functions
# - https://stackoverflow.com/questions/4116658/faster-numpy-cartesian-to-spherical-coordinate-conversion
# - https://stackoverflow.com/questions/41699494/how-to-obscure-a-line-behind-a-surface-plot-in-matplotlib

def appendSpherical_np(xyz):
    
    ptsnew = np.hstack((xyz, np.zeros(xyz.shape)))
    xy = xyz[:,0]**2 + xyz[:,1]**2
    
    # r
    ptsnew[:,3] = np.sqrt(xy + xyz[:,2]**2)
    
    # theta: elevation angle measured from XY plane (+ve upwards)
    ptsnew[:,4] = np.arctan2(xyz[:,2], np.sqrt(xy)) 
    
    # phi
    ptsnew[:,5] = np.arctan2(xyz[:,1], xyz[:,0])
    
    return ptsnew

def plot_visible(azimuth, elev, points, X_coord, Y_coord, Z_coord, quivers, U_coord, V_coord, W_coord, spd, plot_quiver_flag):
    
    # transform viewing angle to normal vector in data coordinates
    a = azimuth*np.pi/180. -np.pi
    e = elev*np.pi/180. - np.pi/2.
    
    X = [ np.sin(e) * np.cos(a),np.sin(e) * np.sin(a),np.cos(e)] 
    
    # concatenate coordinates
    Z = np.c_[X_coord, Y_coord, Z_coord]
    
    # calculate dot product 
    # the points where this is positive are to be shown
    cond = (np.dot(Z,X) >= 0)
    
    # filter points by the above condition
    x_c = X_coord[cond]
    y_c = Y_coord[cond]
    z_c = Z_coord[cond]
    s_c = spd[cond]
    u_c = np.multiply(s_c, U_coord[cond])
    v_c = np.multiply(s_c, V_coord[cond])
    w_c = np.multiply(s_c, W_coord[cond])
    
    # set the new data points
    points.set_data(x_c, y_c)
    points.set_3d_properties(z_c, zdir="z")
    
    if plot_quiver_flag:
        segs = zip(*[x_c, y_c, z_c, u_c, v_c, w_c])
        reshaped_segs = [[[x,y,z],[u,v,w]] for x,y,z,u,v,w in segs]
        quivers.set_segments(reshaped_segs)
        quivers.set_sort_zpos(z_c)
    
    return cond

def norm_vec_list(vlist):
    
    vlist_norm = np.sqrt(np.power(vlist[:,0], 2) + np.power(vlist[:,1], 2) + np.power(vlist[:,2], 2))
    
    vlist_norm = np.tile(vlist_norm, (3,1))
    vlist_norm = np.transpose(vlist_norm)
    vlist = np.divide(vlist, vlist_norm)
    
    return vlist

### output folders

map_folder = folder + os.sep + "maps"
if not os.path.isdir(map_folder):
    os.mkdir(map_folder)
    
overlay_folder = folder + os.sep + "overlay"
if not os.path.isdir(overlay_folder):
    os.mkdir(overlay_folder)
    
abm_folder = folder + os.sep + "ABM"
if not os.path.isdir(abm_folder):
    os.mkdir(abm_folder)
    
npy_folder = folder + os.sep + "output_npy"
if not os.path.isdir(npy_folder):
    os.mkdir(npy_folder)
    
mat_folder = folder + os.sep + "output_mat"
if not os.path.isdir(mat_folder):
    os.mkdir(mat_folder)

### sample points on sphere and render images for local density and velocity

def sample_and_visualize(tp):
    
    pos = pos_dat[:,:,tp]
    drn = dir_dat[:,:,tp]
    vel = vel_dat[:,tp]

    # convert to spherical coords
    sph_coords = appendSpherical_np(pos)
    theta = sph_coords[:,4]               # latitude in [-pi/2, pi/2]
    phi = sph_coords[:,5]                 # longitude in [-pi, pi]

    # compute KD-Tree for fast neighbor search 
    # note: KD-Tree does not support haversine metric
    kd = KDTree(pos)

    theta_idx = np.argsort(theta)
    theta_sorted = theta[theta_idx]
    phi_idx = np.argsort(phi)
    phi_sorted = phi[phi_idx]

    # compute velocity vector components
    # reference: https://stackoverflow.com/questions/707985/
    # calculate-3d-vector-perpendicular-to-plane-described-by-a-point-and-true-north-h
    origin = [0, 0, 0]
    north_vec = [0, 0, sph_rad]

    east_dir = np.cross(north_vec - pos, pos - origin)  # points due east, tangent to the sphere
    north_dir = np.cross(east_dir, pos - origin)        # points due north, tangent to the sphere

    east_dir = norm_vec_list(east_dir)
    north_dir = norm_vec_list(north_dir)

    vec_east = np.multiply(east_dir[:,0], drn[:,0]) + np.multiply(east_dir[:,1], drn[:,1]) + \
                np.multiply(east_dir[:,2], drn[:,2])
    vec_east = np.multiply(vel, vec_east)    

    vec_north = np.multiply(north_dir[:,0], drn[:,0]) + np.multiply(north_dir[:,1], drn[:,1]) + \
                np.multiply(north_dir[:,2], drn[:,2])
    vec_north = np.multiply(vel, vec_north)

    # sampling mesh
    sample_theta = np.linspace(-np.pi/2, np.pi/2, 180)
    sample_phi = np.linspace(-np.pi, np.pi, 360)
    lats, lons = np.meshgrid(sample_theta, sample_phi)
    lats_flat = lats.ravel()
    lons_flat = lons.ravel()

    sample_X = sph_rad * np.multiply(np.sin(lats_flat + np.pi/2), np.cos(lons_flat + np.pi))
    sample_Y = sph_rad * np.multiply(np.sin(lats_flat + np.pi/2), np.sin(lons_flat + np.pi))
    sample_Z = sph_rad * np.cos(lats_flat + np.pi/2)

    # KDE
    kde = KernelDensity(bandwidth=0.02, metric='haversine')
    kde.fit(np.vstack([theta, phi]).T, sample_weight=np.divide(1.0, num_agents))

    latlon = np.vstack([lats_flat, lons_flat]).T
    density_est = np.exp(kde.score_samples(latlon))
    density_est = density_est.reshape((360,180))

    # compute velocity at sample points
    vel_east = np.zeros(np.size(sample_X))
    vel_north = np.zeros(np.size(sample_X))

    nearest_neighbor_dist = 0.005
    ind, dist = kd.query_radius(np.vstack([sample_X, sample_Y, sample_Z]).T, nearest_neighbor_dist, 
                                count_only=False, return_distance=True)

    for i in range(np.size(sample_X)):

        num_neighbors = np.size(ind[i])
        vec_east_sum = np.sum(vec_east[ind[i]])
        vec_north_sum = np.sum(vec_north[ind[i]])

        if num_neighbors > 0:
            vel_east[i] = np.divide(vec_east_sum, num_neighbors)
            vel_north[i] = np.divide(vec_north_sum, num_neighbors)

    vel_east = vel_east.reshape((360,180))
    vel_north = vel_north.reshape((360,180))

    # shift indices by pi/2 and pi
    vel_east = np.roll(vel_east, (180, 0), axis=(0, 1))
    vel_east = np.fliplr(vel_east)
    vel_north = np.roll(vel_north, (180, 0), axis=(0, 1))
    vel_north = np.fliplr(vel_north)
    
    ###################################
    ## plot maps
    ###################################
    
    ytic_vals = np.array(list(np.linspace(0, 180, 6) - 90))
    xtic_vals = np.array(list(np.linspace(-180, 180, 10)))

    plt.figure(figsize=(8,3*3), dpi=300)

    plt.subplot(311)
    plt.imshow(density_est.T)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=6)
    plt.clim(0, 1)

    plt.yticks(np.linspace(0, 180, 6), ytic_vals.astype(int), fontsize=7)
    plt.xticks(np.linspace(0, 360, 10), xtic_vals.astype(int), fontsize=7)

    plt.xlabel("Longitude", fontsize=8)
    plt.ylabel("Latitude", fontsize=8)
    plt.title("Local Density", fontsize=8)

    plt.subplot(312)
    plt.imshow(vel_north.T)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=6) 
    plt.clim(-v_max, v_max)

    plt.yticks(np.linspace(0, 180, 6), ytic_vals.astype(int), fontsize=7)
    plt.xticks(np.linspace(0, 360, 10), xtic_vals.astype(int), fontsize=7)

    plt.xlabel("Longitude", fontsize=8)
    plt.ylabel("Latitude", fontsize=8)
    plt.title("Velocity (North)", fontsize=8)

    plt.subplot(313)
    plt.imshow(vel_east.T)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=6) 
    plt.clim(-v_max, v_max)

    plt.yticks(np.linspace(0, 180, 6), ytic_vals.astype(int), fontsize=7)
    plt.xticks(np.linspace(0, 360, 10), xtic_vals.astype(int), fontsize=7)

    plt.xlabel("Longitude", fontsize=8)
    plt.ylabel("Latitude", fontsize=8)
    plt.title("Velocity (East)", fontsize=8)

    plt.subplots_adjust(wspace=0.2, hspace=0.4)

    plt.savefig(map_folder + os.sep + repr(tp).zfill(3) + ".png")
    plt.close()
    
    ###################################
    ## plot agents
    ###################################
    
    plot_quiver_flag = False
    
    R = sph_rad
    theta_val = theta_sorted[theta_idx.argsort()] + np.pi/2     # latitude in [0, pi]
    phi_val = phi_sorted[phi_idx.argsort()] + np.pi             # longitude in [0, 2*pi]

    X_coord = R * np.sin(theta_val) * np.cos(phi_val)
    Y_coord = R * np.sin(theta_val) * np.sin(phi_val)
    Z_coord = R * np.cos(theta_val)

    # sphere parameterization
    u, v = np.mgrid[0:2*np.pi:36j, 0:np.pi:18j]
    xs = R*np.cos(u)*np.sin(v)
    ys = R*np.sin(u)*np.sin(v)
    zs = R*np.cos(v)

    fig = plt.figure(figsize=(3,3), dpi=300)
    ax = fig.add_subplot(1,1,1, projection='3d')

    # plot empty plot, with points (without a line)
    points, = ax.plot([], [], [], 'ro', markersize=0.3, alpha=1.0, 
                      fillstyle="full", markerfacecolor="red", markeredgecolor='red', zorder=10)

    quivers = ax.quiver([], [], [], [], [], [], color='lightblue', linewidth=0.8, normalize=False, zorder=5)

    # set initial viewing angles
    azimuth, elev = 75, 21
    ax.set_xlim([-R, R])
    ax.set_ylim([-R, R])
    ax.set_zlim([-R, R])
    ax.view_init(elev, azimuth)

    plot_idx = plot_visible(azimuth, elev, points, X_coord, Y_coord, Z_coord, quivers,
                           drn[:,0], drn[:,1], drn[:,2], vel, plot_quiver_flag)

    if plot_quiver_flag:
        rndr = plt.gcf().canvas.get_renderer()
        quivers.draw(rndr)
        
    fig.canvas.draw_idle()

    ax.plot_surface(xs, ys, zs, linewidth=0.1, zorder=0, edgecolor='gray', color='white', shade=False)

    plt.axis("off")
    plt.savefig(abm_folder + os.sep + repr(tp).zfill(3) + ".png")
    plt.close()
    
    ###################################
    ## plot overlays on sphere
    ###################################
    
    lats_vals = lats + np.pi/2
    lons_vals = lons + np.pi

    x = R * np.sin(lats_vals) * np.cos(lons_vals)
    y = R * np.sin(lats_vals) * np.sin(lons_vals)
    z = R * np.cos(lats_vals)

    fig = plt.figure(figsize=(4*3,3), dpi=300)

    # plot density

    ax1 = fig.add_subplot(1,3,1, projection='3d')

    points, = ax1.plot([], [], [], 'ro', markersize=1.5, alpha=0.5, 
                      fillstyle="full", markerfacecolor="red", markeredgecolor='none', zorder=10)

    azimuth, elev = 75, 21
    ax1.set_xlim([-R, R])
    ax1.set_ylim([-R, R])
    ax1.set_zlim([-R, R])
    ax1.view_init(elev, azimuth)

    ls = LightSource(75, 0)
    rho_colors = ls.shade(density_est, cmap=cm.viridis, blend_mode='soft', vert_exag=1)
    rho_plt = ax1.plot_surface(x, y, z, rstride=1, cstride=1, linewidth=0, edgecolor='white',
                          facecolors=rho_colors, antialiased=False, shade=True)

    cbar1 = fig.colorbar(rho_plt, ax=ax1, shrink=0.75)
    cbar1.ax.tick_params(labelsize=7) 
    cbar1.mappable.set_clim(0, 1)

    plt.title("Local Density", fontsize=8)
    plt.axis("off")

    # plot velocity (north)

    ax2 = fig.add_subplot(1,3,2, projection='3d')

    points, = ax2.plot([], [], [], 'ro', markersize=1.5, alpha=0.5, 
                      fillstyle="full", markerfacecolor="red", markeredgecolor='none', zorder=10)

    azimuth, elev = 75, 21
    ax2.set_xlim([-R, R])
    ax2.set_ylim([-R, R])
    ax2.set_zlim([-R, R])
    ax2.view_init(elev, azimuth)

    ls = LightSource(75, 0)
    vel_north_colors = ls.shade(vel_north, cmap=cm.viridis, blend_mode='soft', vert_exag=1)
    vel_north_plt = ax2.plot_surface(x, y, z, rstride=1, cstride=1, linewidth=0, edgecolor='white',
                          facecolors=vel_north_colors, antialiased=False, shade=True)

    cbar2 = fig.colorbar(vel_north_plt, ax=ax2, shrink=0.75)
    cbar2.ax.tick_params(labelsize=7) 
    cbar2.mappable.set_clim(-v_max, v_max)

    plt.title("Velocity (North)", fontsize=8)
    plt.axis("off")

    # plot velocity (east)

    ax3 = fig.add_subplot(1,3,3, projection='3d')

    points, = ax3.plot([], [], [], 'ro', markersize=1.5, alpha=0.5, 
                      fillstyle="full", markerfacecolor="red", markeredgecolor='none', zorder=10)

    azimuth, elev = 75, 21
    ax3.set_xlim([-R, R])
    ax3.set_ylim([-R, R])
    ax3.set_zlim([-R, R])
    ax3.view_init(elev, azimuth)

    ls = LightSource(75, 0)
    vel_east_colors = ls.shade(vel_east, cmap=cm.viridis, blend_mode='soft', vert_exag=1)
    vel_east_plt = ax3.plot_surface(x, y, z, rstride=1, cstride=1, linewidth=0, edgecolor='white',
                          facecolors=vel_east_colors, antialiased=False, shade=True)

    cbar3 = fig.colorbar(vel_east_plt, ax=ax3, shrink=0.75)
    cbar3.ax.tick_params(labelsize=7) 
    cbar3.mappable.set_clim(-v_max, v_max)

    plt.title("Velocity (East)", fontsize=8)
    plt.axis("off")

    # save figure

    plt.savefig(overlay_folder + os.sep + repr(tp).zfill(3) + ".png")
    plt.close()
    
    ###################################
    ## save results
    ###################################

    feat_vec = np.vstack([lats.ravel(), lons.ravel(), density_est.ravel(), vel_north.ravel(), vel_east.ravel()]).T

    output_fname = npy_folder + os.sep + repr(tp).zfill(3) + ".npy"
    savedict = {'feat_vec':feat_vec, 'theta':sample_theta, 'phi':sample_phi}
    np.save(output_fname, savedict)

    output_fname = mat_folder + os.sep + repr(tp).zfill(3) + ".mat"
    sio.savemat(output_fname, savedict)
    
## parallel computation
Parallel(n_jobs=num_cores-1)(delayed(sample_and_visualize)(tp) for tp in range(5));
