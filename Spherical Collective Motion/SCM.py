#!/home/labmember/anaconda2/envs/py36/bin/python

# Spherical Collective Motion (collective motion of swarming agents on a sphere manifold)
# Author: Dhananjay Bhaskar <dhananjay_bhaskar@brown.edu>
# Last Modified: Jan 27, 2021

# Python implementation of Wei Li's framework (DOI: 10.1038/srep13603)

# To generate time-lapse movies:
# ffmpeg -r 10 -f image2 -pattern_type glob -i "*?png" -vcodec libx264 -crf 20 -pix_fmt yuv420p output.mp4
# ffmpeg -i output.mp4 -filter:v "crop=1200:1200:600:600" zoomed_in.mp4

import os
import math
import numpy as np

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d

np.random.seed(1337)

def f(u, v):

    v_norm_sqr = np.linalg.norm(v, 2)**2
    udotv = np.sum(np.multiply(u, v))

    proj = u - np.divide(np.multiply(udotv, v), v_norm_sqr)
    magnitude = np.linalg.norm(proj, 2)

    return np.divide(proj, magnitude)

def random_init(r, n):

    vec = np.random.randn(3, n)
    vec /= np.linalg.norm(vec, axis=0)

    pos = r * np.transpose(vec)

    normal_vec = np.vstack([2*pos[:,0], 2*pos[:,1], 2*pos[:,2]])
    normal_vec = np.transpose(normal_vec)

    const = np.multiply(normal_vec[:,0], pos[:,0]) + np.multiply(normal_vec[:,1], pos[:,1]) + \
                        np.multiply(normal_vec[:,2], pos[:,2])

    x = np.random.rand(n)
    y = np.random.rand(n)

    z = np.divide(const - np.multiply(normal_vec[:,0], x) - np.multiply(normal_vec[:,1], y), normal_vec[:,2])

    tangent_x = x - pos[:,0]
    tangent_y = y - pos[:,1]
    tangent_z = z - pos[:,2]

    tangent_vec = np.vstack([tangent_x, tangent_y, tangent_z])
    tangent_vec = np.transpose(tangent_vec)
    tangent_norm = np.sqrt(np.power(tangent_vec[:,0], 2) + np.power(tangent_vec[:,1], 2) + \
                           np.power(tangent_vec[:,2], 2))

    tangent_norm = np.tile(tangent_norm, (3,1))
    tangent_norm = np.transpose(tangent_norm)
    tangent_vec = np.divide(tangent_vec, tangent_norm)

    return (pos, tangent_vec)

def mapped_Euclidean_step_size(r, v_i):

    return np.multiply(r, np.tan(np.divide(v_i, r)))

def calc_zeta(d, p, r_0, GCD, r, eta, noise_flag):

    n = np.shape(p)[0]

    zeta = np.zeros((n, 3))

    for i in range(n):

        neighbor_set_idx = []

        for j in range(n):
            if GCD[i,j] < r_0:
                neighbor_set_idx.append(j)

        set_size = len(neighbor_set_idx)

        GCR = [0.0, 0.0, 0.0]

        for idx in neighbor_set_idx:

            if noise_flag == False:
                GCR += f(d[idx,:], p[j,:])
            else:
                theta = np.random.uniform(-eta, eta)
                x_j = f(d[idx,:], p[j,:])
                y_j = np.cross(p[j,:], x_j)
                y_j = np.divide(y_j, r)
                GCR += np.multiply(np.cos(theta), x_j) + np.multiply(np.sin(theta), y_j)

        zeta[i,:] = np.divide(GCR, set_size)

    return zeta

def plot_particles(r, pos, drn, vel, fname):

    phi = np.linspace(0, np.pi, 20)
    theta = np.linspace(0, 2 * np.pi, 40)
    xs = r * np.outer(np.sin(theta), np.cos(phi))
    ys = r * np.outer(np.sin(theta), np.sin(phi))
    zs = r * np.outer(np.cos(theta), np.ones_like(phi))

    fig, ax = plt.subplots(1, 1, figsize=(8,8), dpi=300, subplot_kw={'projection':'3d'})
    ax.plot_wireframe(xs, ys, zs, color='gray', rstride=1, cstride=1, alpha=0.6, linewidth=1.0)
    # ax.plot_surface(xs, ys, zs, rstride=1, cstride=1, color='c')

    ax.scatter(pos[:,0], pos[:,1], pos[:,2], s=8, c='r', zorder=10)
    ax.quiver(pos[:,0], pos[:,1], pos[:,2], \
              np.multiply(vel, drn[:,0]), np.multiply(vel, drn[:,1]), np.multiply(vel, drn[:,2]),
              normalize=False, linewidth=0.8)
    plt.axis("off")
    plt.savefig(fname)
    plt.close()

r = 0.25         # Default: 1.0
v_0 = 0.02
n = 800
alpha = 1
eta = math.pi/3  # Default: 0.1
r_0 = 0.05       # Default: 0.4
noise_flag = True

assert(r_0 < r * (math.pi/2))

(pos, drn) = random_init(r, n)
vel = v_0*np.ones((n))

outFolder = "SCM_sim_test_5"
if not os.path.isdir(outFolder):
    os.mkdir(outFolder)

n_t = 400
pos_mat = np.zeros((n, 3, n_t))
dir_mat = np.zeros((n, 3, n_t))
vel_mat = np.zeros((n, n_t))

for k in range(1, n_t):

    step_size = mapped_Euclidean_step_size(r, vel)
    step_size = np.tile(step_size, (3,1))
    step_size = np.transpose(step_size)

    newpos = pos + np.multiply(step_size, drn)

    newpos_norm = np.sum(np.abs(newpos)**2, axis=-1)**(1./2)
    newpos_norm = np.tile(newpos_norm, (3,1))
    newpos_norm = np.transpose(newpos_norm)

    newpos = np.multiply(r, np.divide(newpos, newpos_norm))

    GCD = np.zeros((n,n))
    for i in range(n):
        for j in range(i):
            dot_prod = np.sum(np.multiply(pos[i,:], pos[j,:]))
            GCD[i,j] = r * np.arccos(np.divide(dot_prod, r*r))
            GCD[j,i] = GCD[i,j]

    zeta = calc_zeta(drn, newpos, r_0, GCD, r, eta, noise_flag)

    newdrn = np.zeros((n,3))

    for i in range(n):
        newdrn[i,:] = f(zeta[i,:], newpos[i,:])
        vel[i] = v_0 * (np.linalg.norm(zeta[i,:], 2)**alpha)

    out_fname = outFolder + os.sep + repr(k).zfill(3) + ".png"
    plot_particles(r, pos, drn, vel, out_fname)

    pos = newpos
    drn = newdrn

    pos_mat[:,:,k-1] = newpos
    dir_mat[:,:,k-1] = newdrn
    vel_mat[:,k-1] = vel

    print(out_fname)

savedict = {'r':r, 'v_0':v_0, 'n':n, 'alpha':alpha, 'eta':eta, 'r_0':r_0, 'noise_flag':noise_flag,
            'pos_mat':pos_mat, 'dir_mat':dir_mat, 'vel_mat':vel_mat}
np.save("SCM_sim_data_5.npy", savedict)
