"""

matrixtools.py

A collection of functions used mainly for manipulating and analyzing distance
matrices constructed from IGS data.

"""

import warnings
import numpy as np
import pandas as pd
import source.const as const
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy.interpolate import griddata
from matplotlib.patches import Rectangle 

def init_empty_matrix(bins):
    """ Helper function to generate a matrix of empty lists
    
    Params:
    ----------
        bins: integer bin resolution of matrix

    Returns:
    ----------
        A: matrix of empty lists of resolution equal to bins """
    
    A = []
    for i in range(bins):
        A.append([])
        for j in range(bins):
            A[i].append([np.nan])

    return(A)

def flatten_matrix(A, func = np.nanmean):
    """ Helper function to flatten a matrix of lists

    Params:
    -------
        A: matrix of lists
        func: function to be applied to each matrix element

    Returns:
    --------
        A_flat: 2D matrix flattened according to func """
    
    L = len(A)
    A_flat = np.zeros((L,L))

    for i in range(L):
        for j in range(L):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                A_flat[i][j] = func(A[i][j])
            
    return A_flat

def get_bin_sizes(resolution, genome="hg38"):
    """
    Helper function to generate a list of bin sizes corresponding to the 
    chromosomes in a target genome. Each bin size reflects the number of bins
    that a given chromosome is discretized to at a given bin resolution.
    
    Params:
    -------
        resolution: the resolution at which the chromosomes will be binned,
                    in basepairs
        genome: genome type used to look up chromosome sizes
    
    Returns:
    --------
        bin_sizes_chr: list of bin sizes corresponding to each chromosome in the
                       target genome.
    """
    #Init constants
    SIZES = const.get_genome_sizes(genome)

    #Bin sizes for matrix per chromosome
    chr_count = len(SIZES)
    bin_sizes_chr = np.array([0])
    for i in range(1, chr_count):
        chr_size = SIZES[i]
        num_bins = int(np.ceil(SIZES[i]/resolution))
        bin_sizes_chr = np.append(bin_sizes_chr, num_bins)

    return bin_sizes_chr

def get_B_vector(cn,resolution,
                 sizes=const.SIZES_HG38):
    """
    Helper function to generate a list of bin edges for a chromosome of a
    given genomic length, discretized at an input resolution.
    
    Params:
    -------
        cn: chromosome number
        resolution: resolution used to calculate bin edges
        sizes: the sizes of the chromosomes in the target genome,
               indexed by chromosome number, in basepairs
    Returns:
    --------
        B: vector of bin edges
    """
    
    chr_size = sizes[cn]
    num_bins = int(np.ceil(chr_size/resolution))
    B = np.arange(0,num_bins)
    
    return B

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
    P1 = np.array(c1[KEYS["pos"]].values)
    P2 = np.array(c2[KEYS["pos"]].values)


    #Compute the pairwise distances
    R_cdist = distance.cdist(R1,R2)

    return (R_cdist, P1, P2)

def populate_tile(GM, cni, cnj, ci, cj, resolution):
    """
    Helper function to add pairwise distances from a single cell into a 
    population genome-wide distance matrix by indexing into a "tile" 
    corresponding to a submatrix of two chromosome pairwise distances.

    Params:
    -------
        GM: genome-wide distance matrix, matrix of lists
        cni, cnj: chromosome numbers of the first and second chromosomes
        ci, cj: dataframes for the first and second chromosomes
        resolution: the matrix resolution in basepairs

    Returns:
    --------
        GM: input matrix with new entries appended
        
    """

    #Get bin sizes for all chromosomes in the target genome
    bin_sizes_chr = get_bin_sizes(resolution)

    #Get bin edges for the individual chromosomes 
    B1 = get_B_vector(cni,resolution)
    B2 = get_B_vector(cnj,resolution)

    #Get pairwise distances and corresponding genomic positions
    (R_cdist, P1, P2) = get_cdist(ci, cj)
    
    #Get the bin indices of the observed genomic positions
    P1_inds = np.digitize(P1, B1*resolution)-1
    P2_inds = np.digitize(P2, B2*resolution)-1

    #Offset the indicies by the cumsum of bins preceding current chromosome
    offset1 = np.sum(bin_sizes_chr[:cni])
    offset2 = np.sum(bin_sizes_chr[:cnj])

    #Populate the distance matrix
    for i in range(R_cdist.shape[0]):
        for j in range(R_cdist.shape[1]):
            ii = P1_inds[i] + offset1 
            jj = P2_inds[j] + offset2
            if R_cdist[i][j] == 0: continue
            GM[ii][jj].append((R_cdist[i][j]))
            GM[jj][ii].append((R_cdist[i][j]))

    return GM

def make_genome_wide_matrix(cells,
                            resolution = 10*10**6,
                            genome="hg38"):
    """
    Make a population ensemble genome wide distance matrix given a list of 
    single cells. For each cell, iterate over chromosomes (distinguishing 
    between homologs), compute pairwise distances between their corrsponding
    reads, and index into population matrix to append distances.

    Params:
    -------
        cells: list of cells, dataframe
        resolution: matrix resolution, in basepairs
        genome: target genome, to retreive constants
    
    Returns:
    --------
        GM: the genome wide distance matrix; each pixel is a list of distances
            observed at that corresponding pair of genomic positions, binned
            at the input resolution.
    """

    #Init constants
    SIZES = const.get_genome_sizes(genome) #chromosome sizes
    KEYS = const.get_genome_keys(genome) #dataframe keys
    chr_count = len(SIZES)
    
    #Bin sizes for matrix per chromosome
    bin_sizes_chr = get_bin_sizes(resolution,genome)
    
    #Total and cumulative bins
    total_bins = np.sum(bin_sizes_chr)

    #Init genome wide distance matrix
    GM = init_empty_matrix(total_bins)
    
    for cell in cells:
        for i in range(1,chr_count):
            for j in range(i,chr_count):
                #Handle multiple clusters (e.g. homologs) for a chromosome 
                if i == j:
                    chro = cell.loc[cell[KEYS["chr"]] == i]
                    cluster_idxs = chro[KEYS["cluster"]].unique()
                    
                    for idx in cluster_idxs: #get intra cluster distances

                        ci = chro.loc[chro[KEYS["cluster"]]==idx]
                        if len(ci) < 1: continue #need more than 1 read
                        GM = populate_tile(GM, i, j, ci, ci,resolution)

                #Nonhomologous clusters
                else:
                    ci = cell.loc[cell[KEYS["chr"]] == i]
                    cj = cell.loc[cell[KEYS["chr"]] == j]
                    if len(ci) < 1 or len(cj) < 1: continue
                    GM = populate_tile(GM, i, j, ci, cj,resolution)

    return GM

def draw_genome_wide_matrix(A,
                            xlabel = "\nGenomic Coordinate [Mb]",
                            clabel = '\nSpatial Distance [um]',
                            resolution=10*10**6,
                            q = 0.01,
                            genome = 'hg38'):
    """
    Draw a genome wide distance matrix.
    
    Params:
    -------
        A: distance matrix
        xlabel, clabel: x-axis and colorbar labels
        resolution: matrix resolution in basepairs
        q: percentile cutoff for clim
        genome: target genome, to retreive constants

    Returns:
    --------
        fig, ax: the figure and axes where the matrix is drawn
    """

    fig, ax = plt.subplots()

    SIZES = const.get_genome_sizes(genome) #chromosome sizes, bp
    chr_count = len(SIZES)
    
    #Bin sizes for matrix per chromosome
    bin_sizes_chr = get_bin_sizes(resolution)
    total_bins = np.sum(bin_sizes_chr)
    sum_sizes = np.cumsum(bin_sizes_chr)

    clim = get_clims([A],q) #clim = qth and 1-qth percentile values
    clim = (0,clim[1])#startl linear scale at 0
    
    #Draw outline around chromosome territories
    for i in range(2, chr_count+1):
        offset0 = np.sum(bin_sizes_chr[:i-1])
        offset1 = np.sum(bin_sizes_chr[:i])
        
        ax.hlines(offset0-1,offset0,offset1, lw=1, color='black')
        ax.hlines(offset1-1,offset0,offset1, lw=1, color='black')
        ax.vlines(offset0,offset0-1,offset1-1, lw=1, color='black')
        ax.vlines(offset1,offset0-1,offset1-1, lw=1, color='black')

    cmap = plt.get_cmap('seismic_r')
    cmap.set_bad(color='lightgrey') #handle unmapped regions

    #Draw the matrix and interpolate adjacent unmapped regions
    cax = ax.imshow(A,cmap=cmap,interpolation='nearest')
    cax.set_clim(clim)
    cbar = fig.colorbar(cax, label=clabel)

    #Handle tick labels
    x_tick_labels = ["Chr 1       "]
    y_tick_labels = ["1"]
    for i in range(2,chr_count):
        if i == 23:
            x_tick_labels.append("X")
            y_tick_labels.append("X")
        elif i == 24:
            x_tick_labels.append("Y")
            y_tick_labels.append("Y")
        elif i > 1 and i < 9:
            x_tick_labels.append(str(i))
            y_tick_labels.append(str(i))
        elif i == 10:
            x_tick_labels.append(str(" ⋯ "))
            y_tick_labels.append(str(""))
        else:
            x_tick_labels.append(str(""))
            y_tick_labels.append(str(""))

    plt.xticks(sum_sizes)
    plt.yticks(sum_sizes)
    ax.set_xticklabels(x_tick_labels)
    ax.set_yticklabels(y_tick_labels)
    
    ax.set_xlabel(xlabel)
    ax.set_xlim(0,total_bins)
    ax.set_ylim(total_bins,0)
    plt.tight_layout()
    plt.show()   
    
    return fig, ax

    

def make_distance_matrix(cluster,
                         resolution = 2.5*10**6,
                         statistic = np.nanmean,
                         flatten = True,
                         genome = "mm10"):
    """ Function to generate a distance matrix from a chromosome cluster
    
    Params:
    -------
        cluster: pandas dataframe with chromosome cluster information
        resolution: matrix resolution in base pairs
        statistic: function implementing desired distance metrix
        flatten: If true, apply the distance function, 
                 if false, return matrix with variable length 
        genome: human or mouse genome, string

    Returns:
    --------

    A: matrix of empty lists of resolution equal to bins """

    chr_num = cluster["chr"].unique()[0]

    if type(chr_num) != np.int64:
        raise ValueError("Cluster includes multiple chromosomes.")

    SIZES = const.get_genome_sizes(genome)
    KEYS = const.get_genome_keys(genome)
        
    chr_size = SIZES[chr_num]
    num_bins = int(np.ceil(chr_size/resolution)) 
    
    A = init_empty_matrix(num_bins) 
    B = np.arange(0, num_bins+1) #bin vector
    
    #Get spatial position vector and genomic position vector
    R = np.array([cluster[KEYS["x"]].values,
                  cluster[KEYS["y"]].values,
                  cluster[KEYS["z"]].values]).T
    P = np.array(cluster[KEYS["pos"]].values)
    
    #Bin position vector, then populate binned matrix from unbinned matrix of
    #pdists using binned indices
    P_inds = np.digitize(P, B*resolution) - 1    
    R_pdist = distance.pdist(R)
    R_sf = distance.squareform(R_pdist)

    #Do the matrix binning
    for i in range(len(P_inds)):
        for j in range(i+1, len(P_inds)):
            ii, jj = P_inds[i], P_inds[j]
            A[ii][jj].append((R_sf[i][j]))
            A[jj][ii].append((R_sf[i][j]))

    if flatten:
        #Now flatten lists into 2D matrix according to some statistic
        A = flatten_matrix(A, func=statistic)
    
    return A

def make_ensemble_matrix(data, chr_num,
                         resolution = 2.5*10**6,
                         statistic=np.nanmean,
                         genome="mm10"):
    
    """
    Build an ensemble distance matrix for a particular chromosome given a 
    dataframe.
    
    Params:
    -------
         data: the dataframe(must be prefiltered on stage, parent, and chr)
         chr_num: desired chromosome number for ensemble distance
         resolution: matrix resolution in base pairs
         statistic: function to get distance metric 
         genome: species corresponding to the data
         
    Returns:
    --------
        A_ensemble: the ensemble distance matrix.
    
    """
    #Get all chromosome copies
    matrices, clusters = [], []
    chr_data = data.loc[data["chr"] == chr_num]
    cell_indexes = data["cell_index"].unique()

    for ci in cell_indexes:
        cell = chr_data.loc[chr_data["cell_index"]==ci]
        cluster = cell.loc[cell["chr"] == chr_num]
        if len(cluster) > 0: clusters.append(cluster)

    #Get matrix parameters
    SIZES = const.get_genome_sizes(genome)
    chr_size = SIZES[chr_num]
    num_bins = int(np.ceil(chr_size/resolution))

    #Make the marix
    A_ensemble = init_empty_matrix(num_bins)
    for cluster in clusters:
        A = make_distance_matrix(cluster, resolution=resolution, flatten=False)
        for i in range(len(A)):
            for j in range(len(A)):
                A_ensemble[i][j] = np.concatenate((A_ensemble[i][j],A[i][j]))

    #Flatten the matrix using distance metric function
    A_ensemble = flatten_matrix(A_ensemble, func=statistic)
    
    return A_ensemble

def check_coverage(A, chr_num, expected = 3):
    """
    Helper function to check the coverage (missing rows) of a distance matrix
    
    Params:
    -------
        A: the 2D matrix
        expected: number of rows to ignore in the coverage calculation
                (default of 2, to exclude centromeric/telomeric gaps)
    Returns:
    --------

        coverage: fraction of valid rows
    """
    rows_missing = 0
    for row in A:
        if np.isnan(row).all(): rows_missing += 1

    coverage = 1 - (rows_missing - expected) / (len(A) - expected)
    
    return coverage

def get_diag_lims(A):
    """
    Helper function to find the first non-nan values on the diagonal of a matrix
    
    Params:
    -------
        diag: a matrix diagonal

    Returns:
    --------
        li, ri: tuple corresponding to the (first, last) non-nan diag indices
    """
    li, ri = 0, len(A)
    for i in range(len(A)):
        if not np.isnan(A[i]).all():
            li = i
            break
    for i in range(len(A)-1,-1,-1):
        if not np.isnan(A[i]).all():
            ri = i
            break
    return li, ri

def interpolate_missing(A):
    """
    Helper function to interpolate missing rows in a distance matrix while 
    respecting boundary conditions (e.g. centromere).

    Params:
    -------
        A: the distance matrix
    
    Returns:
    --------
        A: the interpolated distance matrix

    """

    li, ri = get_diag_lims(A)
    diag = np.diagonal(A).copy()

    #interpolate diagonal separately (to handle rows with single reads)
    nan_inds = np.nonzero(np.isnan(diag))[0]
    valid_inds = np.nonzero(~np.isnan(diag))[0]
    diag[nan_inds] = np.interp(nan_inds, valid_inds, diag[valid_inds])

    #interpolate diagonal nans except for those in missing rows
    for i in range(len(A)):
        if i >= li and i <= ri: A[i,i] = diag[i]

    #attribution: Gilly
    #https://stackoverflow.com/questions/6518811/interpolate-nan-values-in-a-numpy-array
    x, y = np.indices(A.shape)
    points = (x[~np.isnan(A)], y[~np.isnan(A)])
    values = A[~np.isnan(A)]
    xi = (x[np.isnan(A)], y[np.isnan(A)])
    A[np.isnan(A)] = griddata(points, values, xi, method='linear')
    
    return A

def get_clims(matrices, q, round_to_int=True):
    """
    Get clims for visualizing a distance matrix, or consensus clims for
    comparisons across multiple matrices.

    Params:
    -------
        matrices: list of distance matrices
        q: percentile cutoff for dynamic range, i.e. (q, 100-q)
        round: optionally round to nearest integer
    Returns:
    --------
        (lower, upper): tuple corresponding to clims
    """
    lowers, uppers = [], []
    for A in matrices:
        lowers.append(np.nanpercentile(A, q = q))
        uppers.append(np.nanpercentile(A, q = 100 - q))

    lower, upper = np.mean(lowers), np.mean(uppers)
    if round_to_int: lower, upper = np.round(lower,0), np.round(upper,0)
    
    return (lower, upper)
    
def draw_scd_bars(axs, peaks):
    
    """
    Draw visually distinct rectangles corresponding to SCD peaks underneath 
    distance matrix.
    
    Params:
    -------
        axs: axes to draw on, first is the distance matrix and second is the scd bars
        peaks: locations of the scd boundary calls
        
    Returns:
    --------
        axs: axes with scd bars drawn underneath the matrix (pre-alignment)
    """
    
    for i in range(1, len(peaks)):
        
        left, right = peaks[i-1], peaks[i]
    
        rect = Rectangle(
            xy=(left, 0.8),
            width = right - left - 1,
            height = 0.05,
            facecolor = const.DISTINCT_COLORS[i],
            edgecolor = 'black',
            linewidth = 0.5
        )
        axs[1].add_patch(rect) 
        
    return axs

def align_scd_bars(fig, axs ,cax, A, res):
  
    """
    Align SCD bars to a distance matrix.
    
    Params:
    -------
        fig: fig containing distance matrix and scd bars
        axes: axes to draw on, 1st is the distance matrix and 2nd is the scd bars
        cax: distance matrix image
        A: raw distance matrix
        res: resolution of the matrix, in megabases
        
    Returns:
    --------
        fig, axs: fig and ax with colorbars aligned to matrix
    """
    
    axs[1].set_xlim(0,len(A)*res)
    
    for loc in ['top','bottom','left','right']:
        axs[1].spines[loc].set_visible(False)
        
    axs[1].tick_params(axis='x', which='both',bottom=False,labelbottom=False)
    axs[1].tick_params(axis='y', which='both',left=False,labelleft=False)

    #ersatz colorbar needed to align with real colorbar
    cbar2 = fig.colorbar(cax, ax=axs[1],label='\nSpatial Distance [um]')
    #align on aspect ratio
    asp = np.diff(axs[1].get_xlim())[0] / np.diff(axs[1].get_ylim())[0]
    axs[1].set_aspect(asp)
    
    plt.tight_layout()
    #ersatz colorbar no longer needed
    cbar2.remove()
    
    return fig, axs

def draw_distance_mat(fig, ax, A, res, clim,
                     chr_num = 11,
                     cmap_str = 'seismic_r',
                     nancolor = 'grey'):
    
    """
        Draw a distance matrix formatted for display in a figure.
    
    Params:
    -------
        fig: fig for distance matrix and scd bars
        axes: axes to draw on, 1st is the distance matrix and 2nd is the scd bars
        A: raw distance matrix
        res: resolution of the matrix, in megabases
        clim: colormap limits
        
    Returns:
    --------
        fig, axs: fig and ax with colorbars aligned to matrix
    """
    
    cmap = plt.get_cmap(cmap_str)
    cmap.set_bad(color=nancolor) #nan pixels, i.e. unmapped regions
    
    extent = [0, res*len(A), res*len(A),0]


    #draw the distance matrix
    cax = ax.imshow(A, cmap=cmap, extent=extent, norm=mpl.colors.LogNorm())
    cax.set_clim(clim)

    #draw the colorbar
    cbar = fig.colorbar(cax, ax=ax)
    
    ax.set_xticks([0,np.ceil(res*len(A))])
    ax.set_yticks([np.ceil(res*len(A)),0])
    
    label = "\n Chr " + str(chr_num) + " Genomic Coordinate [Mb]"
    ax.set_xlabel(label)
    ax.set_ylabel(label)
    
    return fig, ax, cax


