"""
comp.py

A collection of functions used to manipulate and IGS data primarily with respect
to comparisons with other expternal datasets.

"""

import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle 



def get_T(clusters, B, resolution, threshold = 0.5):
    """
    Get ensemble lamin contact probability (i.e. distances below threshold) 
    for a number of bins spanning a chromosome.

    Params:
    -------
        clusters: single chromosome copies, dataframe
        B: bin vector at target resolution
        resolution: target resolution in basepairs
        threshold: threshold for contact in microns
    
    Returns:
    --------
        T: lamin contact probability for each bin
    """
    lamin_tracks = []
    for i in range(len(B)): lamin_tracks.append([])

    #Get the lamin distances for each genomic bin
    for k in range(len(clusters)):
        cluster = clusters[k]
        P = [row["pos"] for index, row in cluster.iterrows()]
        L = [row["dist_to_lamin"] for index, row in cluster.iterrows()]

        P_inds = np.digitize(P, B*resolution) - 1

        for i in range(len(P_inds)):
            lamin_tracks[P_inds[i]].append(L[i])

    lamin_tracks = np.array(lamin_tracks)

    #Get the lamin contact probabilites for each bin
    T = np.zeros(len(B))

    for i in range(len(lamin_tracks)):
        l = np.array(lamin_tracks[i])
        if len(l) == 0: T[i] = np.nan
        else:
            at = len(np.nonzero(l>threshold)[0]) #above threshold
            bt = len(np.nonzero(l<=threshold)[0]) #below threshold            
            T[i] = bt/(bt+at)   
            
    #Mean-center the lamin contact probabilities to get the lamin proximity
    #score.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        T = T - np.nanmean(T)
    
    return T

def draw_comp_plot(E, s_inds, e_inds, T, B, resolution, chr_num):
    """
    Draw the plot comparing Hi-C Eig, DamID LAD, and IGS proximity.
    
    Params:
    -------
        E: vector of Hi-C eigenvalues as a fn of genomic position
        s_inds: start indices of LADs
        e_inds: end incides of LADs
        T: observed lamin proximity probabilities as fn of genomic pos
        B: bins corresponding to a given genomic size and resolution
        resolution: bin resolution in basepairs
        chr_num:: chromosome of interest
    
    Returns:
    --------
        fig, axs: the figure with subplots for each modality
    """

    fig, axs = plt.subplots(3)
    
    hic_ax, lad_ax, spatial_ax = 0, 1, 2

    axs[hic_ax].set_xlim(0,len(E)-1)
    axs[lad_ax].set_xlim(0,len(E)-1)
    axs[spatial_ax].set_xlim(0,len(E)-1)

    axs[hic_ax].set_ylim(-0.15,0.15)
    axs[lad_ax].set_ylim(0,1)
    axs[spatial_ax].set_ylim(-0.5,0.5)
    
    axs[hic_ax].set_ylabel('Hi-C Eig')
    axs[lad_ax].set_ylabel('LAD')
    axs[spatial_ax].set_ylabel('Lamin spatial proximity')    
    axs[spatial_ax].set_xlabel("Chr " + str(chr_num) + \
                               " Genomic Coordinate [Mb]")
    
    axs[hic_ax].spines['top'].set_visible(False)
    axs[hic_ax].spines['bottom'].set_visible(False)
    axs[hic_ax].spines['right'].set_visible(False)
    axs[lad_ax].spines['top'].set_visible(False)
    axs[lad_ax].spines['bottom'].set_visible(False)
    axs[lad_ax].spines['left'].set_visible(False)
    axs[lad_ax].spines['right'].set_visible(False)    
    axs[spatial_ax].spines['top'].set_visible(False)
    axs[spatial_ax].spines['right'].set_visible(False)
    
    axs[hic_ax].tick_params(axis='x', which='both',bottom=False,
                            labelbottom=False)
    axs[lad_ax].tick_params(axis='x', which='both',bottom=False,
                            labelbottom=False)
    axs[lad_ax].tick_params(axis='y', which='both',left=False,
                            labelleft=False)

    Ex = np.arange(0, len(E))
    z = np.zeros(len(E))
    
    axs[hic_ax].plot(Ex,E, color='black',linewidth=1,alpha=0.5)
    axs[hic_ax].plot(Ex,z, color='black',linewidth=1)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        axs[hic_ax].fill_between(Ex,z, E, where= E >=z,interpolate=True,
                                 label='B compartment')
        axs[hic_ax].fill_between(Ex,z, E, where= E < z,interpolate=True,
                                 label='A compartment')
    
    axs[spatial_ax].plot(B,T,color='black', linewidth=1,alpha=0.5)
    axs[spatial_ax].plot(B,z, color='black',linewidth=1)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        axs[spatial_ax].fill_between(B,z, T, where= T > z,interpolate=True,
                                     label='lamin proximal')
        axs[spatial_ax].fill_between(B,z, T, where= T < z,interpolate=True,
                                     label = 'lamin distal')
    
    #Draw lads as rectangles
    for i in range(0,len(e_inds)):
        if i == 0:
            rect = Rectangle(
                xy=(s_inds[i], 0.5),
                width = e_inds[i] - s_inds[i],
                height = 0.2,
                edgecolor = 'black',
                linewidth = 0.5,
                label = 'LAD'
                )
        else:
            rect = Rectangle(
                xy=(s_inds[i], 0.5),
                width = e_inds[i] - s_inds[i],
                height = 0.2,
                edgecolor = 'black',
                linewidth = 0.5,
                )

        axs[lad_ax].add_patch(rect) 
    
    return fig, axs
