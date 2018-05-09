import pandas as pd
from sklearn import manifold
from sklearn.metrics.pairwise import pairwise_distances
import numpy as np
import matplotlib.pyplot as plt
import os
import glob

def tSNEgen(MF, CDs, tol_fact, filter, metric='chebyshev', fetch_ann='online', p=30, ea=12):
    """Performs tSNE analysis on the molecular data collected using spaceM.
    The documentation page of the sklearn implementation of tSNE:
    http://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html

    Args:
        MF (str): path to the Main Folder.
        CDs (list): correlation distance tresholds used for filtering background annotation images, only used when
            filter is 'correlation'. Default value is 0.75.
        tol_fact (float): tolerance factor to use for the filter 'mean'.
        filter (str): filter strategy to select background and on-sample annotation images:
            'mean': compare the mean ion intensity from off and on-sample pixels. Consider annotation as coming from
                the sample if  mean on-sample intensity > tol_fact * mean off-sample intensity.
            'correlation': compute the correlation distance between the intensity thresholded annotation image and
                the cell distribution binary mask. The annotation is considered as coming from the sample if the
                correlation distance is inferior to CDs[i]. The cell distribution mask has pixel equal to 1 if its
                corresponding ablation mark is touching a cell and 0 if not. The treshold value to binarize the
                annotation image is found using an optimizer which minimzes the correlation distance with the cell
                distribution mask. This removes the negative effect that an uneven ion intensity distribution will
                have on the correlation distance with the cell distribution mask.
        metric (str): The metric to use when calculating distance between instances in a feature array. Metric value
            must be one of the options allowed by scipy.spatial.distance.pdist for its metric parameter, or a metric
            listed in pairwise.PAIRWISE_DISTANCE_FUNCTIONS.
        fetch_ann (str): method for fetching annotations:
            'online': (default) queries metaspace using the name of the .imzml data present in the MALDI input folder
                as dataset name,
            'offline': reads annotation images from a provided dataframe.
        p (float): perplexity value to use for the tSNE algorithm. The perplexity is related to the number of nearest
            neighbors that is used in other manifold learning algorithms. Larger datasets usually require a larger
            perplexity. Consider selecting a value between 5 and 50. The choice is not extremely critical since t-SNE
            is quite insensitive to this parameter.
        ea (float): early exaggeration value to use for the tSNE algorithm. Controls how tight natural clusters in the
            original space are in the embedded space and how much space will be between them. For larger values,
            the space between natural clusters will be larger in the embedded space. Again, the choice of this
            parameter is not very critical. If the cost function increases during initial optimization, the early
            exaggeration factor or the learning rate might be too high.

    """

    if fetch_ann == 'online' and filter == 'correlation':
        MOLcsv_p = MF + 'Analysis/scAnalysis/Molecular_features/CD={}/MOLonlyData.csv'.format(CDs[0])
        MOLallcsv_p = MF + 'Analysis/scAnalysis/Molecular_features/CD={}/MOLallData.csv'.format(CDs[0])
    elif fetch_ann == 'online' and filter == 'mean':
        MOLcsv_p = MF + 'Analysis/scAnalysis/Molecular_features/tol_fact={}/MOLonlyData.csv'.format(tol_fact)
        MOLallcsv_p = MF + 'Analysis/scAnalysis/Molecular_features/tol_fact={}/MOLallData.csv'.format(tol_fact)
    if fetch_ann == 'offline':
        MOLcsv_p = MF + 'Analysis/scAnalysis/Molecular_features/offline/MOLonlyData.csv'
        MOLallcsv_p = MF + 'Analysis/scAnalysis/Molecular_features/offline/MOLallData.csv'
    MOLdf = pd.read_csv(MOLcsv_p)
    fluo_data = pd.read_csv(MOLallcsv_p).fluoMarksMean_lu.as_matrix()
    tsne_input = np.nan_to_num(np.log10(MOLdf.iloc[:,2:].as_matrix()))
    # perp = [5,10,15,20,25,30,40,50,75,100]
    # for p in perp:
    # p = 30
    dist = pairwise_distances(tsne_input, metric=metric)
    tsne = manifold.TSNE(n_components=2, metric='precomputed', early_exaggeration=ea, perplexity=p)
    X_tsne = tsne.fit_transform(np.nan_to_num(dist))
    x_min, x_max = np.min(X_tsne, 0), np.max(X_tsne, 0)
    X2_tsne = (X_tsne - x_min) / (x_max - x_min)
    tsne_colors_i = fluo_data

    contrast_cut = 5
    pc_top = np.percentile(tsne_colors_i, 100-contrast_cut)
    pc_down = np.percentile(tsne_colors_i, contrast_cut)

    tsne_colors_f = []
    for i in tsne_colors_i:
        if i >= pc_top:
            val = pc_top
        elif i <= pc_down:
            val = pc_down
        else:
            val = i
        tsne_colors_f = np.append(tsne_colors_f, val)

    plt.figure()
    plt.scatter(X2_tsne[:, 0], X2_tsne[:, 1], 50, np.log10(tsne_colors_f), cmap='viridis', edgecolors='none')
    plt.xlabel('tSNE dim 1', fontsize=20)
    plt.ylabel('tSNE dim 2', fontsize=20)
    plt.axis('equal')
    coords = pd.DataFrame({'tSNE1': X2_tsne[:, 0],
                  'tSNE2': X2_tsne[:, 1],
                  'ObjectNumber': MOLdf['ObjectNumber_lu']})
    coords.to_csv(MF + 'Analysis/tSNE/metric={}_perp={}_KLD='.format(metric, p) +
                  str(tsne.kl_divergence_)[:5] + '_' + fetch_ann + '.csv', index=False)
    plt.savefig(MF + 'Analysis/tSNE/metric={}_perp={}_KLD='.format(metric, p) +
                str(tsne.kl_divergence_)[:5]  + '_' + fetch_ann + '.png', dpi=200)
    plt.close('all')

def genCYTinput(MF):
    """Generate csv input for the CYT MATLAB package. CYT is a graphical software package from Dana Pe’er lab at
    Columbia University lab that combines multiple analysis and visualization tools, such as viSNE, Wanderlust and DREMI
    Link: http://www.c2b2.columbia.edu/danapeerlab/html/cyt.html

    Args:
        MF (str): path to Main Folder.

    """

    MFA = MF + 'Analysis/'
    os.chdir(MFA + 'tSNE/')
    MORPHnMOL_df = pd.read_csv(MFA + 'scAnalysis/MORPHnMOL.csv')

    for str in glob.glob('*.csv'):
        tSNExy_df = pd.read_csv(str)
        fName = str.split('.csv')[0]
        if not os.path.exists(fName):
            os.makedirs(fName)
        os.rename(MFA + 'tSNE/' + fName + '.csv', MFA + 'tSNE/' + fName + '/' + fName + '.csv')
        os.rename(MFA + 'tSNE/' + fName + '.png', MFA + 'tSNE/' + fName + '/' + fName + '.png')
        CYT_in = pd.concat([MORPHnMOL_df.set_index('ObjectNumber'), tSNExy_df.set_index('ObjectNumber')], axis=1)
        CYT_in.reset_index().to_csv(fName + '/cyt_input.csv')
        CYT_in.corr(method='spearman').to_csv(fName + '/CorrMatrix_ClusterGrammer_input.txt', sep="\t")



