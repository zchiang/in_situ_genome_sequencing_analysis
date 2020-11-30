"""
scaling.py

A collection of functions used to manipulate and analyze IGS data primarily
with respect to their scaling properties (e.g. genomic vs. spatial distance).

"""

import warnings
import numpy as np
import source.const as const
import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy.stats import ks_2samp


def get_cell_clusters(cell, chr_nums,genome="hg38"):

    """
    Return clusters (i.e. single chromosome copies) from a single cell given
    a list of chromosome numbers.

    Params:
    -------
        cell: cell of interest, dataframe
        chr_nums: chromosomes of interest, list of ints
        genome: target genome, to retreive constants, string

    Returns:
    --------
        cell_clusters: list of dataframes corresponding to single chr copies
    """
    
    KEYS = const.get_genome_keys(genome)
    
    cell_clusters = []
    
    for chr_num in chr_nums:
        chro = cell.loc[cell[KEYS["chr"]] == chr_num]
        if len(chro) == 0: continue #e.g. x chromsome not present
        cluster_nums = chro[KEYS["cluster"]].unique()

        if genome == "mm10":
            for cluster_num in cluster_nums:
                cell_clusters.append(chro.loc[chro[KEYS["cluster"]] == \
                                                   cluster_num])
                                              
        #Annoying but necessary logic due to cluster labeling in fibroblast data
        elif genome == "hg38":
            clusters_temp = [] 
            for cluster_num in cluster_nums:
                clusters_temp.append(chro.loc[chro[KEYS["cluster"]] == \
                                              cluster_num])

            clusters = sorted(clusters_temp,key=len,reverse=True)

            #If there are three or more clusters, discard all but the largest
            #two, corresponding to the putative chromosome territories. The
            #smaller clusters are the outliers.

            for i in range(len(clusters)):
                if len(clusters) > 1 and i < 2:
                    cell_clusters.append(clusters[i])

        else:
            raise ValueError("Genome not found.")
        
    return cell_clusters


def cluster_helper(cells,chr_nums,genome="hg38"):
    """
    Helper function to get a list of all of the clusters of interest from
    a set of cells, as well as separate lists of clusters of interest
    corresponing to each individual cell.
    
    Params:
    -------
        cells: list of cells, dataframes
        chr_nums: list of chromosomes of interest, ints
        genome: target genome, to retreive constants, string
    
    Returns:
    --------
        clusters: list of all clusters of interest, dataframes
        cells_clusters: list of clusters for each cell, dataframes
    """
    
    clusters, cells_clusters = [], []
    
    for cell in cells:
        cell_clusters = get_cell_clusters(cell,chr_nums,genome)
        for c in cell_clusters:
            clusters.append(c)
        cells_clusters.append(cell_clusters)
        
        
    return clusters, cells_clusters


def cell_cluster_helper(cell, genome="mm10"):
    """
    Helper function to retrieve lists of pairwise spatial and genomic dists.
    
    Params:
    -------
        cell: single cell, dataframe
    
    Returns:
    --------
        R_cell, P_cell: lists of spatial and genomic distances for the cell
    """
    
    R_cell, P_cell = [], []
    for cluster in cell:
        R_pdist, P_pdist = get_pdists(cluster,genome)
        R_cell.extend(R_pdist)
        P_cell.extend(P_pdist)
    
    return R_cell, P_cell 


def get_pdists(cluster,genome="hg38"):
    """
    Get all pairwise euclidean distances within a cluster (e.g. single 
    chromsome copy, chromosome arm), as well as their genomic distances.
    Params:
    -------
        cluster: reads of interest, dataframe
        genome: target genome, to retreive constants, string

    Returns:
    --------
        R_pdist: list of all pairwise euclidean distances
        P_pdist: list of all pairwise genomic distances
     """
    
    KEYS = const.get_genome_keys(genome)

    #Get spatial position vector and genomic position vector
    R = np.array([cluster[KEYS["x"]].values,
                  cluster[KEYS["y"]].values,
                  cluster[KEYS["z"]].values]).T
    P = np.array(cluster[KEYS["pos"]].values)
    
    
    R_pdist = distance.pdist(R)
    P_pdist = distance.pdist(np.array([P]).T)
    
    return R_pdist, P_pdist


def get_cdist(c1, c2, genome = "hg38"):
    """
    Helper function to compute the pairwise spatial distance between two
    clusters (single chromosome copy).
    
    Params:
    -------
        c1: first cluster, dataframe    
        c2: second cluster, dataframe
        genome: genome type used to return dataframe keys, string

    Returns:
    --------
        R_cdist: pairwise distance between the two clusters
        P1: list of genomic positions for the first cluster
        P2: list of genomic positions for the second cluster
    """
    #Init genome-specific dataframe keys
    KEYS = const.get_genome_keys(genome)

    #Get spatial position vectors for the clusters
    R1 = np.array([c1[KEYS["x"]].values,
                   c1[KEYS["y"]].values,
                   c1[KEYS["z"]].values]).T
    
    R2 = np.array([c2[KEYS["x"]].values,
                   c2[KEYS["y"]].values,
                   c2[KEYS["z"]].values]).T

    #Get genomic position vectors for the clusters
    P1 = np.array([c1[KEYS["pos"]].values]).T
    P2 = np.array([c2[KEYS["pos"]].values]).T

    #Compute the pairwise distances
    R_cdist = distance.cdist(R1,R2)
    P_cdist = distance.cdist(P1,P2)

    return (R_cdist, P_cdist)


def get_arm_distances(clusters, chr_num):
    """
    Find the pairwise genomic and spatial distances for each arm of a
    chromosome, aggregated over a list of chromosomes.

    Params:
    -------
        clusters: list of chromosomes of interest, dataframes
        chr_num: chromosome number needed to retrieve centromere consts, int
    
    Returns:
    --------
        R_intra, R_inter: the pairwise spatial distances within and between
                          chromsome arms
        P_intra, P_inter: the pairwise genomic distances within and between
                          chromosomes arms
    """

    #left and right centromere bounds
    centromere_pos = const.CENTROMERES_HG38[chr_num]

    #spatial and genomic distance within and between arms
    R_inter, P_inter, R_intra, P_intra = [],[],[],[]

    for cluster in clusters:
        #break cluster up into arms
        qarm = cluster.loc[cluster["hg38_pos"] < centromere_pos[0]]
        parm = cluster.loc[cluster["hg38_pos"] > centromere_pos[1]]

        #intra-arm distances
        R_pdist, P_pdist = get_pdists(parm)
        R_intra.append(R_pdist)
        P_intra.append(P_pdist)

        R_pdist, P_pdist = get_pdists(qarm)
        R_intra.append(R_pdist)
        P_intra.append(P_pdist)

        #inter-arm distances
        if len(parm) > 0 and len(qarm) > 0:
            R_cdist, P_cdist = get_cdist(parm, qarm)
            if len(R_cdist.ravel())>0: R_inter.append(R_cdist.ravel())
            if len(P_cdist.ravel())>0: P_inter.append(P_cdist.ravel())

    #combine the cluster distances
    R_intra, P_intra = np.concatenate(R_intra), np.concatenate(P_intra)
    R_inter, P_inter = np.concatenate(R_inter), np.concatenate(P_inter)

    return R_inter, P_inter, R_intra, P_intra


def draw_arm_curves(x, y, yerr):
    """
    Draw the genomic vs spatial distance curves for inter- and intra-arm dists.
    Params:
    -------
        x: tuple of arrays, inter and intra_arm genomic distances
        y: tuple of arrays, inter and intra_arm spatial distances
        yerr: tuple of arrays, spatial distance standard deviations

    Returns:
    --------
        fig, ax: the figure and axes for the plot
    """

    fig, ax = plt.subplots()
    
    ax.errorbar(x[0],y[0], yerr=yerr[0], linestyle="", marker='o', markersize=4,
                alpha=0.25, capthick=1, capsize=2, label="inter-arm",c='#1f77b4')

    ax.plot(x[0],y[0], linestyle="", marker='o', markersize=4, alpha=1, 
            c='#1f77b4',markeredgewidth=0)

    ax.errorbar(x[1],y[1], yerr=yerr[1], linestyle="", marker='o', markersize=4,
                alpha=0.25, capthick=1, capsize=2, label="intra-arm", c='r')

    ax.plot(x[1],y[1], linestyle="", marker='o', markersize=4, alpha=1, 
            c='r',markeredgewidth=0)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xlabel("Genomic Distance [Mb]")
    ax.set_ylabel("Mean Spatial Distance [μm]")
    
    return fig, ax


def draw_arm_violins(R_inter, R_intra):
    """
    Draw the distribution of inter and intra arm spatial distances.
    
    Params:
    -------
        R_inter, R_intra: the inter and intra arm distances
    
    Returns:
    --------
        fig, ax: the fig and axes for the plot
    """
    fig, ax = plt.subplots()

    for i in range(len(R_intra)):
        y = R_intra[i]
        x = np.random.normal(1, 0.02)
        ax.plot(x, y,  marker='.', color='r', markersize='2', alpha=0.05)
        ax.plot(x, y, marker='.', color='r', markersize='2', alpha=0.2)

    for i in range(len(R_inter)):
        y = R_inter[i]
        x = np.random.normal(2, 0.02)
        ax.plot(x, y, marker='.', color='#1f77b4', markersize='2', alpha=0.05)
        ax.plot(x, y, marker='.', color='#1f77b4', markersize='2', alpha=0.2)

    violins = ax.violinplot([R_intra, R_inter], vert=True,
                                 showmedians=True,showextrema=True)

    body = violins['bodies'][0]
    body.set_facecolor('r')
    body.set_edgecolor('r')
    
    violins['cmedians'].set_color(['r','#1f77b4'])
    violins['cmins'].set_color(['r','#1f77b4'])
    violins['cmaxes'].set_color(['r','#1f77b4'])
    violins['cbars'].set_color(['r','#1f77b4'])
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_ylabel("Pairwise Spatial distance (μm)")
    ax.set_xlabel("Intra-arm                  Inter-arm")

    #test significance
    p = ks_2samp(R_intra, R_inter)[1]
    
    title = "p = {:.2e} (K-S test)".format(p)
    ax.set_title(title)

    return fig, ax
