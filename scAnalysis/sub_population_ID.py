# ----------------------------------------------------------------
# Attempt to create t-sne plot.
# 1 -  get annotation image pixel indexes which overlap with one/multiple cell
# 2 -  fetch annotation names which a correltaion distance with cohobin inferior to CHOOSE: 0.8/1 and a fdr = CHOOSE: 0.2
# 3 -  Loop over the annotation names, get ion image and create vectors for annotation image.
# 4 -  Input the vectors in embedding algorithms


from sm_analytics_python.sm_annotation_utils import sm_annotation_utils as smau
import glob, os
from scipy.spatial.distance import correlation
import numpy as np
import matplotlib.pyplot as plt
import tifffile as tiff
from sklearn import manifold
from sklearn.cluster import DBSCAN, AgglomerativeClustering, MeanShift
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import pairwise_distances
import matplotlib.gridspec as gridspec
from scipy import spatial
from scipy.optimize import basinhopping
from scipy.stats import pearsonr, spearmanr
from skimage.morphology import erosion,dilation,disk
import tifffile as tiff
import pandas as pd

MFA = MF + 'Analysis/'

marksMask = np.load(MFA + 'SURF/transformedMarksMask.npy')
cellMask = tiff.imread(MFA + 'CellProfilerAnalysis/Labelled_cells.tif')
fluo = tiff.imread(MFA + 'CellProfilerAnalysis/img_t1_z1_c2.tif')
bf = tiff.imread(MFA + 'CellProfilerAnalysis/img_t1_z1_c1.tif')
window = 100
fdr = 0.5
corr_dit_thres = 0.5
#img  = plt.imread(MFA + 'CellProfilerAnalysis/Contour_cells.png')
coordX, coordY = np.load(MFA + 'SURF/transformedMarks.npy')

foldername = MFA + 'tSNE/'
if not os.path.exists(MFA + 'tSNE/'):
    os.makedirs(MFA + 'tSNE/')

sm = smau.SMInstance()
os.chdir(MF + 'Input/MALDI/')
ds_name = glob.glob('*.imzML')[0].replace('.imzML', '')
d = sm.dataset(ds_name)
results1 = sm.msm_scores([d], d.annotations(0.5, database='HMDB'), db_name='HMDB').T
results2 = sm.msm_scores([d], d.annotations(0.5, database='ChEBI'), db_name='ChEBI').T
results = pd.concat([results1, results2]).drop_duplicates()

# cohibin_img = d.isotope_images('C15H17N5O6S2', '-H')
# cohi_vec255 = np.rot90(cohibin_img[0],2).T
# plt.imshow(img)
#
# for key, value in norm_MM.items():
#         print(key)
#         plt.show()
#         plt.scatter(norm_MM[key]['y'], norm_MM[key]['x'], 20, cohi_vec255.ravel(), edgecolors='none')

# if not os.path.exists(MFA + 'tSNE/' + 'marks_flitered_fluo.npy'):
norm_MM = {} #Normalized markMask -->
for i, data in enumerate(marksMask):
    print(i)
    norm_MM[str(i)] = {}
    norm_MM[str(i)]['x'] = np.array(marksMask[i, 0] - np.min(coordX) + window).astype(np.int64)
    norm_MM[str(i)]['y'] = np.array(marksMask[i, 1] - np.min(coordY) + window).astype(np.int64)

plt.imshow(cellMask)
plt.imshow(fluo)

cell_marks = {}
cell_fluo = {}
marks_fluo = {}

for i, label in enumerate(np.unique(cellMask)):
    if label > 0:

        cell_marks[str(i)] = []
        cell_fluo[str(i)] = []
        cellMask_bw = cellMask == label
        #
        # print(np.mean(fluo[cellMask == label]))
        for key, value in norm_MM.items():
            if True in np.unique(cellMask_bw[norm_MM[key]['x'], norm_MM[key]['y']]):
                print(i, key)

                plt.show()
                marks_fluo[str(i)] = np.append(cell_fluo[str(i)], np.mean(fluo[norm_MM[key]['x'], norm_MM[key]['y']]))#TODO: remove np.fliplr
                cell_marks[str(i)] = np.append(cell_marks[str(i)], key)
                cell_fluo[str(i)] = np.append(cell_fluo[str(i)], np.mean(fluo[cellMask == label]))#TODO: remove np.fliplr
                #
                plt.scatter(norm_MM[key]['y'], norm_MM[key]['x'],np.mean(fluo[cellMask == label])/20)

#Loop over annotation image and create violin plots of ion int. vs binned fluo int.
# nbins = 15
# fluo_bins = np.histogram(fluo, nbins)[1]
# tree = spatial.KDTree(list(zip(np.zeros(np.shape(fluo_bins)).ravel(), fluo_bins.ravel())))
# for i, row in enumerate(results.reset_index().itertuples()):
#     images = d.isotope_images(row.sf, row.adduct)
#     an_vec255 = (np.fliplr(images[0]).ravel())
#     plot_data = {}
#     print(i)
#     for j in fluo_bins:
#         plot_data[str(j)] = []
#     for key, value in cell_marks.items():
#         if not np.shape(cell_fluo[key])[0] == 0:
#             dist, ind = tree.query(list(zip(np.zeros(np.shape(cell_fluo[key])).ravel(),cell_fluo[key])))
#             # plt.scatter(np.ones(np.shape(cell_fluo[key]))*fluo_bins[ind], an_vec255[[int(i) for i in cell_marks[key]]] )
#             plot_data[str(fluo_bins[ind][0])] = np.append(plot_data[str(fluo_bins[ind][0])], an_vec255[[int(i) for i in cell_marks[key]]])
#
#     keys =  plot_data.keys()
#     del_key = []
#     for key in keys:
#         if np.shape(plot_data[key])[0] == 0:
#             del_key = np.append(del_key, key)
#         elif np.mean(plot_data[key]) == 0:
#             del_key = np.append(del_key, key)
#     for key in del_key:
#         del plot_data[key]
#
#     pos = np.array([float(key) for key, value in plot_data.items()])
#     data = np.array([plot_data[key][plot_data[key] > 0] for key, value in plot_data.items()])
#     data_mean = np.array([np.mean(value) for key, value in plot_data.items()])
#     pear_score = pearsonr(pos[np.where(data_mean > 0)[0]], data_mean[np.where(data_mean > 0)[0]]) # Pearson correlation
#     plt.figure()
#     figManager = plt.get_current_fig_manager()
#     figManager.window.showMaximized()
#     vio_data = plt.violinplot(data, pos, widths=10,
#                           showmeans=True)
#     plt.ylabel('Ion intensity (A.U)', fontsize = 50)
#     plt.xlabel('Cell mean fluorescence intensity (A.U)', fontsize = 50)
#     # plt.yscale("log", nonposy='clip')
#     plt.title(row.sf + ' Pearson = ' + str(pear_score[0])[:6], fontsize = 50)
#     plt.savefig('C:/Users/Luca/Desktop/GM6/pearson=' + str(pear_score[0])[:6]
#                 + '_sf=' + row.sf + '.png', format='png', dpi=100)
#     plt.close('all')

# cohibin_img = d.isotope_images('C37H68O4', '+H')
# cohi_vec255 = (np.fliplr(cohibin_img[0]).ravel())
np.save(MFA + 'tSNE/' + 'marks_flitered_fluo.npy', [norm_MM, cell_marks, cell_fluo, marks_fluo])
# else:
#     norm_MM, cell_marks, cell_fluo, marks_fluo = np.load(MFA + 'tSNE/' + 'marks_flitered_fluo.npy')


# if not os.path.exists(MFA + 'tSNE/' + 'tsne_inputs.npy'): #TODO uncomment this

cellMask_bw_all = cellMask > 0
# plt.scatter(Y-np.min(Y),X-np.min(X))
pmi = []  # Positive Mark Index
mbi = []  # Mark brightfield intensity
overLaps = []
for i in range(np.shape(marksMask)[0]):
    status = 0
    cell_mark_OL = 0
    bi = []
    for j in range(np.shape(marksMask[i][0])[0]):
        bi = np.append(bi, bf[
            int(marksMask[i][0][j] - np.min(coordX) + window), int(marksMask[i][1][j] - np.min(coordY) + window)])
        if cellMask_bw_all[
            int(marksMask[i][0][j] - np.min(coordX) + window), int(marksMask[i][1][j] - np.min(coordY) + window)]:
            status = 1
            cell_mark_OL += 1
            # if status == 1:
    pmi = np.append(pmi, status)
    mbi = np.append(mbi, np.mean(bi))
    overLaps = np.append(overLaps, cell_mark_OL)
    print(i)

pmi = np.reshape(pmi, [int(np.sqrt(coordX.shape[0])), int(np.sqrt(coordX.shape[0]))]).ravel()
mbi = np.reshape(mbi, [int(np.sqrt(coordX.shape[0])), int(np.sqrt(coordX.shape[0]))]).ravel()


def err_func(params, pmi, an_vec255):
    threshold = params
    corr = correlation(pmi, an_vec255 > threshold)
    return corr

    # step = []
    # corr = []
    # for i in np.linspace(0, np.max(an_vec255)-1, 100):
    #     corr = np.append(corr, correlation(pmi, an_vec255 > i))
    #     step = np.append(step, i)
    #     print(i)
    # print(step[corr == np.min(corr)])
    #
    # minF = basinhopping(err_func, x0=[947], \
    #                     niter=50, T=10, stepsize=10, \
    #                     minimizer_kwargs={'args': ((pmi, an_vec255))}, \
    #                     take_step=None, accept_test=None, callback=None, interval=200, disp=True, \
    #                     niter_success=50)
    # print(minF.x[0])


def running_mean(l, N):
    sum = 0
    result = list( 0 for x in l)

    for i in range( 0, N ):
        sum = sum + l[i]
        result[i] = sum / (i+1)

    for i in range( N, len(l) ):
        sum = sum - l[i-N] + l[i]
        result[i] = sum / N

    return result

# for CD in np.linspace(0, 1, 50):
# pmi = np.flipud(d.isotope_images('C13H8O6', '+H')[0]).ravel() #TODO small hack, remove it after
tsne_data = {}
tsne_wholeCell_fluo = {}
tsne_sf = []
tsne_cellMarks_fluo = {}
for key, value in cell_marks.items():
    if not np.shape(cell_fluo[key])[0] == 0:
        tsne_data[key] = []
        tsne_wholeCell_fluo[key] = np.mean([int(i) for i in cell_fluo[key]])
        tsne_cellMarks_fluo[key] = np.mean([int(i) for i in marks_fluo[key]])
best_corrS = []
for i, row in enumerate(results.reset_index().itertuples()):
    images = d.isotope_images(row.sf, row.adduct)
    an_vec255 = np.rot90(np.flipud(images[0]),3).ravel() #TODO: np.fliplr
    step = []
    corr = []
    for j in np.linspace(0, np.max(an_vec255) - 1, 100):
        corr = np.append(corr, correlation(pmi, an_vec255 > j))
        step = np.append(step, j)
        # print(i)
    threshold = np.mean(step[corr == np.min(corr)])
    minF = basinhopping(err_func, x0=[threshold], \
                        niter=50, T=10, stepsize=10, \
                        minimizer_kwargs={'args': ((pmi, an_vec255))}, \
                        take_step=None, accept_test=None, callback=None, interval=200, disp=False, \
                        niter_success=50)
    # print()
    threshold = minF.x[0]
    best_corr = correlation(pmi, an_vec255 > threshold)
    best_corrS = np.append(best_corrS, best_corr)
    print(i, best_corr, correlation(pmi, an_vec255), minF.x[0])
np.save(MFA + 'tSNE/best_corrS.npy', [best_corrS])

CDs = [0.5, 0.6, 0.7, 0.8, 2]
for CD in CDs:
    CD_fname = MFA + 'tSNE/CD={}/'.format(CD)
    if not os.path.exists(CD_fname):
        os.makedirs(CD_fname)

    for i, best_corr in enumerate(best_corrS):
        if best_corr < CD:
            sf = results.reset_index().as_matrix()[i, 0]
            adduct = results.reset_index().as_matrix()[i, 1]
            an_vec255 = np.rot90(np.flipud(d.isotope_images(sf, adduct)[0]),3).ravel()
            for key, value in cell_marks.items():
                if not np.shape(cell_fluo[key])[0] == 0:
                    tsne_data[key]  = np.append(tsne_data[key], np.mean(an_vec255[[int(k) for k in cell_marks[key]]]))
            tsne_sf = np.append(tsne_sf, sf)
    key = np.array(list(tsne_data.keys()))[0]
    print(np.shape(tsne_data[key]))

    if np.shape(tsne_data[key])[0]>0:
        # foldername = MFA + 'tSNE/CD_' + str(CD) + '_Nannot = ' + str(np.shape(tsne_data[key])[0]) + '/'
        # if not os.path.exists(foldername):
        #     os.makedirs(foldername)
        tsne_input = np.array([val for key, val in tsne_data.items()])
        fluoCell_input = np.array([val for key, val in tsne_wholeCell_fluo.items()])
        fluoMarks_input = np.array([val for key, val in tsne_cellMarks_fluo.items()])
        keys_input = np.array([key for key, val in tsne_data.items()])

        np.save(CD_fname + 'tsne_inputs_CD=' + str(CD) + '_nAnnot = ' + str(np.shape(tsne_data[key])[0]) + ' .npy',
                [tsne_input, fluoCell_input, fluoMarks_input, keys_input, tsne_sf, tsne_data, best_corrS, CD, pmi, mbi])

            # tsne_input, fluoCell_input, fluoMarks_input, keys_input, tsne_sf, tsne_data, best_corrS, CD, pmi, mbi = np.load(
            #     MFA + 'tSNE/' + 'tsne_inputs.npy')
            # show_corr_plots = 0
            # if show_corr_plots == 1: #TODO, remove this (repeated later)
            #
            #     if not os.path.exists(MFA + 'CellProfilerAnalysis/corr_plots/'):
            #         os.makedirs(MFA + 'CellProfilerAnalysis/corr_plots/')
            #     if not os.path.exists(MFA + 'CellProfilerAnalysis/pdf_plots/'):
            #         os.makedirs(MFA + 'CellProfilerAnalysis/pdf_plots/')
            #     spearmanS = []
            #     pearsonS = []
            #     sfS = []
            #     for i in range(0, np.shape(tsne_data[[key for key in tsne_data][0]])[0]):
            #         ion_ind = i
            #         fluo_noZero = fluoMarks_input[np.where(tsne_input[:, ion_ind] > 0)[0]]
            #         int_noZero = tsne_input[np.where(tsne_input[:, ion_ind] > 0)[0]]
            #         # sf_noZero =  tsne_sf[np.where(tsne_input[:, ion_ind] > 0)[0]]
            #         sort_ind = np.argsort(fluo_noZero)
            #         fluo_rel_sorted = np.sort(fluo_noZero) / np.max(fluo_noZero)
            #         int_sorted = int_noZero[sort_ind, ion_ind]
            #         sf = tsne_sf[ion_ind]
            #         plt.figure()
            #         plt.hist(int_sorted, 15)
            #         plt.xlabel('Ion Intensity', fontsize=20)
            #         plt.ylabel('PDF', fontsize=20)
            #         plt.savefig(MFA + 'CellProfilerAnalysis/pdf_plots/' + sf + '.png', format='png', dpi=100)
            #         plt.close('all')
            #
            #         running_mean_l = 150
            #         plt.figure()
            #         aa = plt.scatter(np.log10(fluo_rel_sorted + 1), np.log10(int_sorted), 50, edgecolors='none')
            #         if np.shape(int_sorted)[0]> running_mean_l:
            #             ab = plt.scatter(np.log10(fluo_rel_sorted + 1), np.log10(running_mean(int_sorted, running_mean_l)), 50, 'r',
            #                              edgecolors='none')
            #         plt.legend((aa, ab), ('Raw data', 'Smoothed'))
            #         plt.title('Putative annotation: ' + sf, fontsize=25)
            #         plt.xlabel('Log10( Relative Cell Fluorescence Intensity )', fontsize=20)
            #         plt.ylabel('Log10( Ion intensity )', fontsize=20)
            #         pear = pearsonr(np.log10(fluo_rel_sorted + 1), np.log10(int_sorted))[0]
            #         spear = spearmanr(np.log10(fluo_rel_sorted + 1), np.log10(int_sorted))[0]
            #         plt.savefig(MFA + 'CellProfilerAnalysis/corr_plots/spearman=' + str(spear) + 'pearson=' + str(
            #             pear) + 'sf=' + sf + '.png', format='png', dpi=100)
            #         plt.close('all')
            #         sfS = np.append(sfS, sf)
            #         spearmanS = np.append(spearmanS, spear)
            #         pearsonS = np.append(pearsonS, pear)
            #     np.save(MFA + 'CellProfilerAnalysis/corr_plots/sfS_spearmanS_pearsonS.npy', [sfS, spearmanS, pearsonS])

    # else:
    #     tsne_input, fluoCell_input, fluoMarks_input, keys_input, tsne_sf, tsne_data, best_corrS, CD = np.load(MFA + 'tSNE/tsne_inputs_CD=' + str(CD) + '.npy')

    show_corr_plots = 1
    if show_corr_plots == 1:

        corrPlots_fname = CD_fname + 'Correlation_plots/'
        pdfPlots_fname = CD_fname + 'PDF_plots/'

        if not os.path.exists(corrPlots_fname):
            os.makedirs(corrPlots_fname)
        if not os.path.exists(pdfPlots_fname):
            os.makedirs(pdfPlots_fname)

        for i in range(0, np.shape(tsne_data[[key for key in tsne_data][0]])[0]):
            ion_ind = i
            fluo_noZero = fluoMarks_input[np.where(tsne_input[:, ion_ind] > 0)[0]]
            int_noZero = tsne_input[np.where(tsne_input[:, ion_ind] > 0)[0]]
            # sf_noZero =  tsne_sf[np.where(tsne_input[:, ion_ind] > 0)[0]]
            sort_ind = np.argsort(fluo_noZero)
            fluo_rel_sorted = np.sort(fluo_noZero) / np.max(fluo_noZero)
            int_sorted = int_noZero[sort_ind, ion_ind]
            sf = tsne_sf[ion_ind]
            plt.figure()
            plt.hist(int_sorted, 15)
            plt.xlabel('Ion Intensity', fontsize=20)
            plt.ylabel('PDF', fontsize=20)
            plt.savefig(pdfPlots_fname + sf + '.png', format='png', dpi=100)
            plt.close('all')

            plt.figure()
            aa = plt.scatter(np.log10(fluo_rel_sorted + 1), np.log10(int_sorted), 50, edgecolors='none')
            ab = plt.scatter(np.log10(fluo_rel_sorted + 1), np.log10(running_mean(int_sorted, 150)), 50, 'r',
                             edgecolors='none')
            plt.legend((aa, ab), ('Raw data', 'Smoothed'))
            plt.title('Putative annotation: ' + sf, fontsize=25)
            plt.xlabel('Log10( Relative Cell Fluorescence Intensity )', fontsize=20)
            plt.ylabel('Log10( Ion intensity )', fontsize=20)
            pear = pearsonr(fluo_rel_sorted, int_sorted)[0]
            spear = spearmanr(fluo_rel_sorted, int_sorted)[0]
            plt.savefig(corrPlots_fname + 'spearman=' + str(spear) + 'pearson=' + str(
                pear) + 'sf=' + sf + '.png', format='png', dpi=100)
            plt.close('all')

            'kulsinski'
    distances = ['canberra', 'chebyshev', 'hamming', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener',
                   'sokalsneath', 'sqeuclidean']
    distances = ['canberra', 'chebyshev', 'seuclidean', 'sqeuclidean']

    for h in distances:
        #h = 'chebyshev'
        print('Computing ' + h)
        dist = pairwise_distances(tsne_input, metric=h)
        # dist = np.nan_to_num(np.clip(dist, 0, 1))
        tsne = manifold.TSNE(n_components=2, metric='precomputed', perplexity = 50, method= 'barnes_hut')
        X_tsne = tsne.fit_transform(dist)
        x_min, x_max = np.min(X_tsne, 0), np.max(X_tsne, 0)
        X2_tsne = (X_tsne - x_min) / (x_max - x_min)
        # plt.figure()
        # plt.imshow(np.dstack([X2_tsne[:,ii].reshape(50,50) for ii in range(3)]), interpolation='nearest')
        # plt.show()
        # labels_true = pmi
        X = X2_tsne
        X = StandardScaler().fit_transform(X)
        # labels_true = pmi
        X = X2_tsne
        X = StandardScaler().fit_transform(X)
        db = DBSCAN(eps=0.13, min_samples=6).fit(X)
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        labels = db.labels_
        tsne_colors_i = fluoMarks_input
        # for sf in tsne_sf:
        # sf_ind = np.where(tsne_sf == sf)[0][0]
        # tsne_colors_i = tsne_input[:,sf_ind]
        pc = np.percentile(tsne_colors_i, 97)

        tsne_colors_f = []
        for i in tsne_colors_i:
            if i >= pc:
                val = pc
            else:
                val = i
            # print(val)
            tsne_colors_f = np.append(tsne_colors_f, val)
        if not os.path.exists(foldername + 'tSNE_CD=' + str(CD) + '/'):
            os.makedirs(foldername + 'tSNE_CD=' + str(CD) + '/')

        plt.figure()
        plt.scatter(X2_tsne[:, 0], X2_tsne[:, 1], 100, np.log10(tsne_colors_f) ,cmap = 'viridis', edgecolors='none')
        plt.xlabel('tSNE dim 1', fontsize=20)
        plt.ylabel('tSNE dim 2', fontsize=20)
        plt.axis('equal')
        plt.title(h)
        # plt.title(sf)
        plt.savefig(CD_fname + h + '_tSNE_fluo.png', dpi = 200)
        plt.close('all')

plt.figure()
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
gs = gridspec.GridSpec(1, 2)
ax0 = plt.subplot(gs[0, 1])
# bandwidth = estimate_bandwidth(X, quantile=0.15)
# ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
# ms.fit(X)
# labels = ms.labels_
# cluster_centers = ms.cluster_centers_

# labels_unique = np.unique(labels)
# n_clusters_ = len(labels_unique)

# print("number of estimated clusters : %d" % n_clusters_)
# number of estimated clusters : 3
# plt.figure(1)
# plt.clf()

# n_clusters = 3
# clustering = AgglomerativeClustering(linkage='ward', n_clusters=n_clusters)
# clustering.fit(X2_tsne)
# labels = clustering.labels_
ax0.scatter(X2_tsne[:, 0], X2_tsne[:, 1], 100, labels, edgecolors='none')
# for k, col in zip(range(n_clusters_), colors):
#     my_members = labels == k
#     cluster_center = cluster_centers[k]
#     ax0.plot(X[my_members, 0], X[my_members, 1], col + '.')
#     ax0.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
#              markeredgecolor='k', markersize=14)

# plt.scatter(X2_tsne[:, 0], X2_tsne[:, 1], 100, clustering.labels_, edgecolors='none')

# if len(set(labels)) > 1:
#     # Number of clusters in labels, ignoring noise if present.
#     n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
#
#     # print('Estimated number of clusters: %d' % n_clusters_)
#     # print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
#     # print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
#     # print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
#     # print("Adjusted Rand Index: %0.3f"
#     #       % metrics.adjusted_rand_score(labels_true, labels))
#     # print("Adjusted Mutual Information: %0.3f"
#     #       % metrics.adjusted_mutual_info_score(labels_true, labels))
#     # print("Silhouette Coefficient: %0.3f"
#     #       % metrics.silhouette_score(X, labels))
#
#
#     unique_labels = set(labels)
#     colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
#     #Plot DBSCAN
#     for k, col in zip(unique_labels, colors):
#         if k == -1:
#             # Black used for noise.
#             col = 'k'
#         class_member_mask = (labels == k)
#         xy = X[class_member_mask & core_samples_mask]
#         ax0.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
#                  markeredgecolor='k', markersize=14)
#         xy = X[class_member_mask & ~core_samples_mask]
#         ax0.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
#                  markeredgecolor='k', markersize=6) #  #
ax1 = plt.subplot(gs[0, 0])
ax1.scatter(X2_tsne[:, 0], X2_tsne[:, 1], 100, np.log10(fluoMarks_input), edgecolors='none')
ax1.set_title('t-SNE, metric = ' + h)
# plt.savefig(MFA + '/tSNE/tSNEt.png')

n_clusters = np.shape(np.unique(labels))[0]-1
for clustN in range(n_clusters):
    print('Cluster ' + str(clustN) + ' = ' + str(np.shape(np.where(labels == clustN))[1]))

# plt.savefig(foldername + h + '.png', format='png',
#             dpi=150)
# plt.close('all')
#find the cluster of highest intensity

from matplotlib import cm
import tifffile as tiff

def maptSNE2microscopy(contour_p, label_img_p, save_p,  keys_input, fluoCell_input, tsne_xy, tsne_labels, cmap = 'brg'):

    cellMask = tiff.imread(label_img_p)
    img_highlighted = plt.imread(contour_p)

    # contour_p = MFA + 'CellProfilerAnalysis/Contour_cells_adjusted.png'
    # cellMask = tiff.imread(MFA + 'CellProfilerAnalysis/Labelled_cells.tif')
    CD = 0.8
    tsne_input, fluoCell_input, fluoMarks_input, keys_input, tsne_sf, tsne_data, best_corrS, CD = np.load(
        MFA + 'tSNE/tsne_inputs_CD=' + str(CD) + '.npy')

    tsne_xy = [np.random.rand(keys_input.size)*10,np.random.rand(keys_input.size)*10]
    tsne_labels = [int(lbl) for lbl in np.random.rand(keys_input.size)*2]

    splot = plt.scatter(tsne_xy[0], tsne_xy[1], 20, tsne_labels, cmap = cmap)
    face_colors = splot.get_facecolors()

    label_256 = [int(lbl) for lbl in (tsne_labels / np.max(tsne_labels)) * 255]
    face_colors = cm.brg([np.unique(label_256)])[0]
    for j, tsne_lbl in enumerate(np.unique(tsne_labels)):
        # print(tsne_lbl)

        cellMask_ind = np.zeros(np.shape(cellMask))
        dim1, dim2 = np.shape(cellMask)
        cell_indexes = keys_input[np.where(tsne_labels == tsne_lbl)[0]]

        for i, cell_ind in enumerate(cell_indexes):
            cellMask_ind[np.where(cellMask == np.unique(cellMask)[int(cell_ind)])] = 1
            print('{}/{}'.format(i, cell_indexes.size - 1))

        cellMask_ind_contour = dilation(cellMask_ind, disk(3)) - erosion(cellMask_ind, disk(3))
        x, y = np.where(cellMask_ind_contour[100:dim1 - 100, 100:dim2 - 100] == 1)
        img_highlighted[:, :, 0][x, y] = face_colors[j, 0]
        img_highlighted[:, :, 1][x, y] = face_colors[
            j, 1]  # lbl_c --> color from tSNE scatter plot of the figure; need to find the command to find it back
        img_highlighted[:, :, 2][x, y] = face_colors[j, 2]
    # save_p = MFA + 'tSNE/tSNE_cluster_colored.tif'
    tiff.imsave(save_p, img_highlighted)
#
# clust_Fluo = []
# for clustN in range(n_clusters):
#     # print(clustN)
#     fluos = []
#     for key in np.array(keys_input[np.where(labels == clustN)]):
#         fluos = np.append(fluos, cell_fluo[key])
#     clust_Fluo = np.append(clust_Fluo, np.mean(fluos))
# highFluoClust_ind = np.where(clust_Fluo == np.max(clust_Fluo))[0][0]
# # Overlay cells from cluster of high fluo int on the FLuo microscopy
# # Green cells-> high fluo int. cluster
# # Red cells -> low fluo int. cluster
# cellInd_highFluoClust = np.array(keys_input[np.where(labels == highFluoClust_ind)])
# cellInd_lowerFluoClust = np.array(keys_input[np.where(labels != highFluoClust_ind)])
# cellMask_highFluoClust = np.zeros(np.shape(cellMask))
# cellMask_lowerFluoClust = np.zeros(np.shape(cellMask))
# for i in cellInd_highFluoClust:
#     label_value = np.unique(cellMask)[int(i)]
#     cellMask_highFluoClust = cellMask_highFluoClust + (cellMask == label_value)
#     print(i)
# for i in cellInd_lowerFluoClust:
#     label_value = np.unique(cellMask)[int(i)]
#     cellMask_lowerFluoClust = cellMask_lowerFluoClust + (cellMask == label_value)
#     print(i)
# dim1, dim2 = np.shape(cellMask)[:2]
# cellMask_highFluoClust_contour = dilation(cellMask_highFluoClust, disk(1)) - erosion(cellMask_highFluoClust,
#                                                                  disk(3))
# cellMask_lowerFluoClust_contour = dilation(cellMask_lowerFluoClust, disk(1)) - erosion(cellMask_lowerFluoClust,
#                                                                  disk(3))
#
#
#
#
# img_highlighted =  plt.imread(MFA + 'CellProfilerAnalysis/Contour_cells.png')
# x1, y1 = np.where(cellMask_highFluoClust_contour[100:dim1 - 100, 100:dim2 - 100] > 0)
# x2, y2 = np.where(cellMask_lowerFluoClust_contour[100:dim1 - 100, 100:dim2 - 100] > 0)
# img_highlighted[:, :, 0][x1, y1] = 0
# img_highlighted[:, :, 1][x1, y1] = 1
# img_highlighted[:, :, 2][x1, y1] = 0
# img_highlighted[:, :, 0][x2, y2] = 0
# img_highlighted[:, :, 1][x2, y2] = 0
# img_highlighted[:, :, 2][x2, y2] = 1
# ax2 = plt.subplot(gs[1, :])
# ax2.imshow(img_highlighted)
# ax2.set_ylim([1000, 2000])
# ax2.set_xlim([0, 3500])
#
#
