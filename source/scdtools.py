"""
scdtools.py

Functions mainly concerned with manipulating and analyzing single-cell domains
in IGS data.

"""

import source.matrixtools as mt
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from source import const
from scipy.signal import find_peaks
from scipy.signal import resample
from matplotlib.patches import Rectangle


def bootstrap_2samp(data1,data2,n=5000,func=np.median):
    """
    Generate n pairs of bootstrap samples from the pair of datasets, evaluating 
    func at each resampling. Returns a funnction that can be called to obtain 
    confidence intervals.

    Params:
    -------
        data1, data2: the data, 1D lists
        n: number of simulations
        func: function to be evaluated in simulation
    
    Returns:
    --------
        ci: function which takes a % CI on [0,1] and returns the two-sided CI
    """
    
    simulations = list()
    sample_size1 = len(data1)
    sample_size2 = len(data2)
    
    for i in range(n):
        resample1 = np.random.choice(data1,size=sample_size1,replace=True)
        resample2 = np.random.choice(data2,size=sample_size2,replace=True)
        simulations.append(func(resample1,resample2))
    simulations.sort()
    
    def ci(p):
        #Return 2-sided symmetric confidence interval specified by p
        
        u_pval = (1+p)/2.
        l_pval = (1-u_pval)
        l_indx = int(np.floor(n*l_pval))
        u_indx = int(np.floor(n*u_pval))
        
        return(simulations[l_indx],simulations[u_indx])
    return(ci)


def cohens_d(d1,d2):
    """
    Calculate effect size by Cohen's d
    
    Params:
    -------
        d1, d2: the two distributions
    
    Returns:
    --------
        d: cohen's d for the two distributions
    """
        
    mean1 = np.mean(d1)
    mean2 = np.mean(d2)
    var1 = np.var(d1,ddof=1)
    var2 = np.var(d2,ddof=1)
    n1 = len(d1)
    n2 = len(d2)
    
    s_pooled = np.sqrt((var1*(n1-1)+var2*(n2-1))/(n1+n2-2))
    
    d = (mean1-mean2)/s_pooled
    
    return d


def sliding_window_dist(A, ws):
    """
    Return an insulation score for each entry along the diagonal of a distance 
    matrix, using a sliding window. Adapted from Su et al. 2020.
           
    Params:
    -------
    
        A: The distance matrix.
        ws: The size of the sliding window, in pixels.

    Returns:
    --------

        dists: The vector of insulation scores along the diagonal of the distance 
               matrix.

    """
    dists = np.zeros(len(A))
    for i in range(len(A)):
        #must consider at least half of window
        if (i - ws//2) < 0 or (i + ws//2) > len(A): continue

        ls = slice(max(0, i - ws), i) #left slice upstream of current entry
        rs = slice(i, min(i+ws,len(A)))

        #get non-nan intra-slice distances
        a = A[ls, ls]
        a = np.triu(a)
        a = a[np.isnan(a) == False]

        b = A[rs,rs]
        b = np.triu(b)
        b = b[np.isnan(b) == False]

        #combined upstream and downstream intra-slice distances
        intra_dist = np.concatenate([a[a>0],b[b>0]])

        #inter-slice distances
        inter_dist = A[ls, rs]
        inter_dist = inter_dist[np.isnan(inter_dist) == False]

        if len(intra_dist) == 0 or len(inter_dist) == 0: continue

        #calculate insulation score, here we use a score analagous to a
        #reflection coefficient, as in Su et al. 2020
        normed_insulation = (np.nanmean(intra_dist) - np.nanmean(inter_dist)) / (np.nanmean(intra_dist) + np.nanmean(inter_dist))

        dists[i] = normed_insulation
        
    return dists


def handle_edges(peaks, A, threshold=3):

    """
    Add matrix boundaries to list of peaks. To compensate for edge noise, also 
    check for peaks near matrix boundaries and remove peak if within a threshold 
    of matrix boundary distance. 
    
    Params:
    -------
        peaks: list of peaks
        A: distance matrix
        threshold: max distance from matrix boundary to remove a peak, in bins
        
    Returns:
    --------
        peaks: boundary-adjusted list of peaks 
    """
    
    li, ri = mt.get_diag_lims(A)

    left_edge_width = peaks[0] - li
    right_edge_width = ri - peaks[-1]
        
    if left_edge_width > threshold:
        peaks = np.insert(peaks,0, li)
    else:
        peaks[0] = li
        
    if right_edge_width > threshold:
        peaks = np.append(peaks,ri)
    else:
        peaks[-1] = ri
    
    return peaks


def get_scd_peaks(A, ws=5, distance=3, prominence=0.05, make_plots=False):
    """
    Find insulation maxima in a distance matrix using a sliding window. Adapted 
    from Su et al 2020.
    
    Params:
    -------
        A: distance matrix
        ws: sliding window size 
        distance: min distance between peaks
        prominence: prominence required for peaks

    Returns:
    --------
        peaks: a list of peaks corresponding to bins in the input matrix

    """
                  
    #get insulation score at each coordinate
    dists = sliding_window_dist(A, ws)
    
    #find maxima of insulation score
    peaks = find_peaks(-dists, distance=distance, prominence=prominence)[0]

    #correct at edges
    peaks = handle_edges(peaks,A)
    
    if make_plots==True:
        fig = plt.figure()
        ax = fig.add_subplot(111)

        x = np.arange(0, len(M))
        ax.plot(x, dists)
        ax.scatter(peaks,dists[peaks], c='red')

        fig = plt.figure()
        ax = fig.add_subplot(111)
        cax = ax.matshow(M,interpolation='nearest',cmap='seismic_r') 
        ax.plot(peaks, peaks, linestyle="", marker = 'o', color='black')   
            
    return peaks, dists


def get_scd_sizes(matrices):
    """
    Helper function to get list of distances between SCD boundaries from a 
    list of chromosome distance matrices
    
    Params:
    -------
        matrices: the list of distance matrices
    
    Returns:
    --------
        sizes: array of distances between pairs of SCD boundaries

    """
    
    sizes = []
    #get the sizes
    for i in range(len(matrices)):
        peaks, dists = get_scd_peaks(matrices[i])
        for j in range(1, len(peaks)):
            size = peaks[j] - peaks[j-1]
            sizes.append(size)

    sizes = np.array(sizes)

    return sizes
    

def draw_sizes_plot(sizes,res,bins,title):
    """
    Helper function to draw the histogram of distances b/w SCD boundaries

    Params:
    ------
        sizes: the list of distance b/w scd boundaries
        res: the matrix resolution in Mb
        bins: number of bins to be plotted in histogram
        title: plot title which may contain additional information

    Returns:
    --------
        fig, ax: the figure and axes for histogram

    """

    fig, ax = plt.subplots()

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel("SCD size [Mb]")
    ax.set_ylabel("Frequency")
    ax.set_title(title)
    
    ax.hist(sizes*res,bins=bins,alpha=0.3,histtype='bar',ec='black')

    return fig, ax


def get_boundary_strengths(data, matrices, num_chrs):
    """
    Helper function to get boundary strengths for ensemble and single cell
    distance matrices.
    Params:
    -------
        data: the dataframe (assumed to be prefiltered on stage, parent, and chr)
        matrices: the list of high-coverage distance matrices
        num_chrs: the chromosome numbers to consider
    
    Returns:
    --------
        sc_strengths, ensemble_strengths: respectively, lists of the boundary 
                      strengths in single cells, aggregated over all high-
                      coverage matrices, and the boundary strengths in ensemble
                      matrices, aggregated over all chromosomes numbers.
    """
    #single cell and ensemble boundary strengths
    sc_strengths, ensemble_strengths = [], []
    
    #compute ensemble boundary strengths first
    for i in range(1,num_chrs+1):
        E = mt.make_ensemble_matrix(data, i) #ensemble distance matrix
        peaks, dists = get_scd_peaks(E) #boundary indices and ins. scores
        dists = dists*-1 #express more insulated boundaries as positive
        #get the boundary strengths; do not consider the edge peaks which are
        #included in the peaks vector for convenience only
        ensemble_strengths.extend(dists[peaks[1:-1]])
        
    #compute single-cell boundary strengths in each high-coverage matrix
    for A in matrices:
        peaks, dists = get_scd_peaks(A) #boundary indices and ins. scores
        dists = dists*-1 #express more insulated boundaries as positive
        #get the boundary strengths; do not consider the edge peaks which are
        #included in the peaks vector for convenience only
        sc_strengths.extend(dists[peaks[1:-1]])

    return sc_strengths, ensemble_strengths
    

def draw_strengths_plot(sc_strengths, ensemble_strengths, title):
    """
    Helper function to draw the violin plots for boundary strengths.

    Params:
    -------
        sc_strengths: list of single-cell boundary strengths
        ensemble_strengths: list of ensemble boundary strengths
        title: plot title which may contain additional information

    Returns:
        fig, ax: the figure and axes for violin plots
    --------
    
    """
    fig, ax = plt.subplots()
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    color = ["green","green"]    
    plot_data = [np.array(ensemble_strengths),np.array(sc_strengths)]
    
    ax = sns.violinplot(data=plot_data,orient='v', palette=color, linewidth=2.5,
                       showfliers=False, showextrema=False)

    ax.set_title(title)
    ax.set_xticklabels(["Ensemble", "Single Cells"])
    ax.set_ylabel("Boundary strength")

    return fig, ax

def draw_scd_ensemble_subplot(Ms, Cs, res):
    """
    Helper function to draw a SCD ensemble plot for a single chromosome number.

    Params:
    -------
        A_i: list of distance matrices corresponding to current chromosome number
        C_i: list of single-chromosome clusters (datframes) corresponding to
             current chromosome number
        res: the matrix resolution, in megabases
    
    Returns:
    --------
        fig, ax: the figure and axes for the SCD ensemble
    """

    #Init the plot
    fig, axs = plt.subplots(len(Ms),1)

    #For each distance matrix
    for j in range(len(Ms)):

        #Get the peaks and convert to megabases
        peaks, dists = get_scd_peaks(Ms[j])
        peaks = peaks*res

        for k in range(1,len(peaks)):
            #Draw a distict color rectangle corresponding to the genomic
            #distance between two SCD peaks.
            rect = Rectangle(
                xy=(peaks[k-1], 0.1),
                width = (peaks[k] - peaks[k-1]),
                height = 0.5,
                facecolor =  const.DISTINCT_COLORS[k+1],
                edgecolor ='black',
                linewidth = 1
            )    
            #Draw the rectangle
            axs[j].add_patch(rect)
            
        axs[j].set_xlim(-1,len(Ms[j])*res)
        axs[j].tick_params(axis='y', which='both',left=False,labelleft=False)

        #Annotate plot with current embryo
        embryo_num = str(Cs[j]["embryo_id"].unique()[0])
        if j == 0: 
            l = axs[j].set_ylabel("Embryo " + embryo_num+"                ")
            l.set_rotation(0)
        else:
            l = axs[j].set_ylabel("\n"+str(embryo_num)+"   ")
            l.set_rotation(0)

        #Only draw x axis label on final row
        if j != len(Ms)-1: axs[j].set_xticks([])
        else:  axs[j].set_xlabel("Genomic Position [Mb]")

    return fig, axs
    

def get_scaled_dists(peaks, P, L, resolution, num_scaled_bins):
    """
    Scale positions of reads b/w two SCD peaks to [0,1] across N bins,
    and return a list of read attributes (e.g. lamin distance) at those scaled
    positions. 

    Params:
    -------
        peaks: list of peak coordinates
        P: genomic position vector of reads
        L: vector of target attribute (e.g. lamin distance), same length as P
        resolution: matrix resolution in basepairs
        num_scaled_bins: number of scaled bins, fixed for all scds

    Returns:
    --------
    
    scaled_dists: vector of bins, each bin contains a read attributes 
                  (e.g. lamin distance) at those scaled positions. 

    """

    #Bins to scale positions of reads b/w two SCD peaks.
    scaled_dists = []
    for i in range(num_scaled_bins): 
        scaled_dists.append([])

    #Do the scaling and record read attribute for each pair of peaks
    for i in range(1, len(peaks)):

        bin_edges = np.linspace(peaks[i-1], peaks[i], num=num_scaled_bins+1)

        # For each genomic pos in the position vector that falls b/w two
        # scd peaks, bin it based on its relative position b/w the peaks, and
        # retrieve its corresponding attribute.
        for j in range(len(P)):

            pos = P[j] # get the genomic position

            #position must fall between the peaks
            if pos >= peaks[i-1] and pos <= peaks[i]:
                
                ind = np.digitize(pos, bin_edges) - 1 #scaled bin index
                
                scaled_dists[ind].append(L[j]) #bin the attribute
    
    return scaled_dists


def get_shifted_scaled_dists(A, peaks, P, L, resolution,
                             num_scaled_bins, size_coeff=10):

    """
    Sample alternate valid placements for SCDs in single-cell distance matrices, 
    and generate their scaled distance profiles to construct an expected 
    distance profile.

    Params:
    -------
        A: distance matrix
        peaks: list of peak coordinates
        P: genomic position vector of reads
        L: vector of target attribute (e.g. lamin distance), same length as P
        resolution: matrix resolution in basepairs
        num_scaled_bins: number of scaled bins, fixed for all scds
        size_coeff: determines number of samples per chromosome, which
                    scales with chromosome size

    Returns:
    --------
    
    scaled_shifted_dists: Expected vector of binned distances, each bin contains 
                          read attributes (e.g. lamin distance) generated from
                          scaling their genomic positions in alternate (shifted)
                          SCD boundary placements.
    """
    
    li, ri = mt.get_diag_lims(A) #left/right indices for ~nan matrix boundaries
    
    #Bin vector to scale positions of reads b/w two SCD peaks
    scaled_shifted_dists = []    
    for i in range(num_scaled_bins): 
        scaled_shifted_dists.append([])

    #For each pair of peaks, sample valid alternate placements in distance
    #matrix and generate control distance profiles.
    for i in range(1, len(peaks)):

        #Keep distance between peaks constant
        scd_size = peaks[i] - peaks[i-1]
        #Min/max valid alternate placements for the left peak
        left_bound, right_bound = li*resolution, ri*resolution - scd_size
        #All valid alternate placements for the left boundary
        vb = np.arange(left_bound, right_bound +resolution, step=resolution)

        #Test a random sample of those boundaries, with sample size proportional
        #to chromosome size.
        sample_size = size_coeff*int(ri-li)
        R = np.random.choice(vb, size = sample_size)

        for r in R:
            #Define the shifted peaks
            left_peak, right_peak = r, r + scd_size
            shifted_peaks = [left_peak, right_peak]
            
            #Get the shifted distance profile
            sd = get_scaled_dists(shifted_peaks, P, L,
                                  resolution, num_scaled_bins)
            #Record the shifted profile
            for j in range(len(scaled_shifted_dists)):
                scaled_shifted_dists[j].extend(sd[j])

    return scaled_shifted_dists


def mirror_distance_vector(scd_scaled_dists):
    """
    Function to take a binned list of scaled scd distances on [0,1] and 
    mirror the right half to result in a binned list on [0,0.5], which is ~twice
    as deep per bin and has half as many bins.
    
    Params:
    ------
        scd_scaled_dists: list of scaled scd distances on [0,1]
    
    Returns:
    --------
        mirrored: binned list of scaled scd distance on [0,0.5]

    """
    
    interval = len(scd_scaled_dists)//2 #assumes even number of bins

    left = scd_scaled_dists[:interval] #left half
    right = scd_scaled_dists[-1:interval-1:-1] #right half, reversed

    #Do the mirroring
    mirrored, mirrored_shifted = [], []
    for i in range(interval):
        mirrored.append([])
        mirrored[i].extend(left[i])
        mirrored[i].extend(right[i])

    return mirrored


def get_oe_distances(observed_dists, expected_dists,
                     func = np.median, CI_bound=0.95):

    """
    Function to take a list of observed and expected scaled scd distance vectors
    and return a single flat observed over expected vector according to some
    function.
    
    Params:
    -------
        observed_dists: binned list of observed distances on scaled interval
        expected_dists: binned list of expected distances on scaled inverval
        func: metric on which to calculate the observed over expected
    
    Returns:
    --------
        oe: flat observed-over-expected vector corresponding to the given
            metric.
        CIs: list of CI tuples corresponding to 95% CI for oe
    """       

    def oe_func(observed,expected,func=np.median):
        o = func(observed)
        e = func(expected)
        return o/e

    OEs, CIs, = [],[]
    for i in range(len(observed_dists)):

        #get observed over expected for function
        o = np.array(observed_dists[i])
        e = np.array(expected_dists[i])
        oe = oe_func(o,e)
        OEs.append(oe)
        
        #get o/e confidence interval
        boot = bootstrap_2samp(o,e,func=oe_func)
        ci = boot(CI_bound)
        CIs.append(ci)

    OEs, CIs = np.array(OEs), np.array(CIs)
    
    return OEs, CIs


def make_oe_lamin_plot(OEs, CIs,bins=22):
    """
    Draw the observed-over-expected lamin profile.
    
    Params:
    -------
        OEs: the observed-over-expected values
        CIs: the confidence intervales for the OEs
        bins: the number of bins used to construct the profile
    
    Returns:
    --------
        fig, ax: figure and axes for the lamin profile.
    """

    fig, ax = plt.subplots()
    x = np.arange(len(OEs))/(bins-2)
    ax.hlines(1,0,1, linestyle="--")
    ax.set_xlim(0,0.5)

    errors = []
    for i in range(len(OEs)):
        errors.append(np.abs(CIs[i]-OEs[i]))

    errors = np.array(errors)

    ax.plot(x, OEs,'-')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(0,0.5)

    ax.hlines(1,0,1, linestyle="--")

    ax.fill_between(x, OEs - errors.T[0], OEs+errors.T[1],
                    color='gray', alpha=0.2)

    ax.set_xlabel("Scaled dist. to domain boundary")
    ax.set_ylabel("Median dist to lamin [obs/exp]")

    return fig, ax
