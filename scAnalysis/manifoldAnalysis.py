import pandas as pd
from sklearn import manifold
from sklearn.metrics.pairwise import pairwise_distances
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import re

def tSNEgen(MF, CDs, tol_fact, filter, metric='chebyshev', fetch_ann='online', p=30, ea=12):
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



