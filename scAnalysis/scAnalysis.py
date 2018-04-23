from matplotlib import cm
from sm_analytics_python.sm_annotation_utils import sm_annotation_utils as smau
import glob
from scipy.spatial.distance import correlation
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import basinhopping
import tifffile as tiff
import pandas as pd
import os
import itertools
import seaborn as sns
import tqdm
import codecs


def defMOLfeatures(MF,
                   tf_obj=ion2fluoTF,
                   CDs=[0.75],
                   tol_fact=-0.2,
                   filter = 'correlation',
                   norm_method='weighted_mean_sampling_area_MarkCell_overlap_ratio_sampling_area',
                   area_prop=0.3,
                   fetch_ann='online',
                   hdf5_path='C:/Users\Luca\Google Drive\A-Team\projects/1c\hepatocytes, '
                             'DKFZ\datasets/Molecular images/2017-09-12-luca-mz-images.hdf5'):

    """Defines molecular intensities of individual cells.

     Args:
         MF (str): path to the Main Folder.
         tf_obj (function): Image transformation to apply on ion image for registration.
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
         tol_fact (float): tolerance factor to use for the filter 'mean'.
         CDs (list): correlation distance tresholds used for filtering background annotation images, only used when
            filter is 'correlation'. Default value is 0.75.
         norm_method (str): normalization method to use:
            'Mean_intensity': mean Intensity of AM touching the cell OI,
            'weighted_mean_sampling_area_MarkCell_overlap_int_power1,7': mean Intensity of AM dvidied by its
                sampling area weighted by the number of pixels of the AM sampling area touching the cell OI,
            'weighted_mean_sampling_area_MarkCell_overlap_int': mean Intensity of AM dvidied by its sampling area
                weighted by the number of pixels of the AM sampling area touching the cell OI,
            'weighted_mean_sampling_area_MarkCell_overlap_ratio_sampling_area': (default) mean Intensity of AM dvidied
                by its sampling area weighted by the ration between AM sampling area touching the cell OI and its
                sampling area,
            'squared_weighted_mean_sampling_area_MarkCell_overlap_ratio_sampling_area': mean Intensity of AM dvidied
                by its sampling area weighted by the ration between AM sampling area touching the cell OI and its
                sampling area,
            'mean_sampling_ratio': mean Intensity of AM dvidied by the proportion of its area sampling a cell,
            'mean_sampling_area': mean Intensity of AM dvidied by its sampling area in um.
         area_prop (float): lower threshold value of overlap between an ablation marks and cell below which the
            ablation mark are discarded.
         fetch_ann (str): method for fetching annotations:
            'online': (default) queries metaspace using the name of the .imzml data present in the MALDI input folder
                as dataset name,
            'offline': reads annotation images from a provided dataframe.
         hdf5_path (str): path to the .hdf5 file containing the annotation images. Only used when fetch_ann='offline'.
            The hdf file should be organized as follows:
            {
            'ds_name': dataset name (str),
            'image': annotation image (2D array),
            'mol_formula': MS1 sum formula (str) (example: C42H82NO6P)
             }

     """

    def getPixSize(MFI):
        """Reads the pixel size in um from the Nikon Ti E microscope (NIS elements software).

        Args:
            MFI (str): path to Main Folder Input.

        Returns:
            pix_size (float): pixel size in um.

        """
        txt_file = codecs.open(MFI + '/Microscopy/postMALDI/out.txt', 'r','utf-16')
        for row in txt_file:
            if row.startswith('Calibration'):
                pix_size = float(row.strip().split()[2].replace(',', '.'))
            else:
                pix_size = 0.73 #Pixel size used in all experiments
        return pix_size

    MFA = MF + 'Analysis/'
    marksMask = np.load(MFA + 'SURF/transformedMarksMask.npy')
    cellMask = tiff.imread(MFA + 'CellProfilerAnalysis/Labelled_cells.tif')
    fluo = tiff.imread(MFA + 'CellProfilerAnalysis/img_t1_z1_c2.tif')
    fluo_nucl = tiff.imread(MFA + 'CellProfilerAnalysis/img_t1_z1_c3.tif')
    bf = tiff.imread(MFA + 'CellProfilerAnalysis/img_t1_z1_c1.tif')
    window = 100
    coordX, coordY = np.load(MFA + 'SURF/transformedMarks.npy')
    os.chdir(MF + 'Input/MALDI/')
    # ds_name = glob.glob('*.imzML')[0].replace('.imzML', '')
    pixel_size = getPixSize(MF + 'Input/')

    Fname = MFA + 'scAnalysis/Molecular_features/'
    if not os.path.exists(Fname):
        os.makedirs(Fname)
    CD = CDs[0]
    # sm = smau.SMInstance()
    os.chdir(MF + 'Input/MALDI/')
    ds_name = glob.glob('*.imzML')[0].replace('.imzML', '')

    norm_MM = {}  # Normalized markMask -->
    for mark_ind, data in enumerate(marksMask):
        # print(i)

        norm_MM[str(mark_ind)] = {}
        norm_MM[str(mark_ind)]['x'] = np.array(marksMask[mark_ind, 0] - np.min(coordX) + window).astype(np.int64)
        norm_MM[str(mark_ind)]['y'] = np.array(marksMask[mark_ind, 1] - np.min(coordY) + window).astype(np.int64)

    # marks_cell = {}
    # Record the cell indexes for the different marks
    # if not os.path.exists(MFA + 'tSNE/' + 'cells_per_marks.npy'):
    #     ind = 1
    #     for key, value in norm_MM.items():
    #         labels_underMM = np.unique(cellMask[norm_MM[key]['x'], norm_MM[key]['y']])
    #         # Record which cells are under which ablation marks
    #         marks_cell[key] = [np.where(np.unique(cellMask) == labeli)[0][0]
    #                            for labeli in labels_underMM[labels_underMM > 0]]
    #         print(key, ind)
    #         ind = ind +1
    #     np.save(MFA + 'tSNE/' + 'cells_per_marks.npy', [marks_cell])
    # else:
    #     marks_cell = np.load(MFA + 'tSNE/' + 'cells_per_marks.npy')

    # nucl_fluo = {}
    # for cell_ind in tqdm.tqdm(np.unique(cellMask)):
    #     if cell_ind > 0:
    #         nucl_fluo[str(cell_ind)] = []  # mean Raw intensity from second channel image
    #         cellMask_bw = cellMask == cell_ind
    #         for mark_ind, value in norm_MM.items():
    #             if True in np.unique(cellMask_bw[norm_MM[mark_ind]['x'], norm_MM[mark_ind]['y']]):
    #                 nucl_fluo[str(cell_ind)] = np.append(nucl_fluo[str(cell_ind)],
    #                                                      np.mean(fluo_nucl[cellMask == cell_ind]))

    if not os.path.exists(Fname + 'marks_flitered_fluo.npy'):

        # plt.imshow(cellMask)
        # plt.imshow(fluo)
        marks_cellLabels = {}
        marks_samplingArea = {}
        mark_area = {}
        #Record the cell indexes for the different marks
        for i, [mark_ind, value] in enumerate(norm_MM.items()):
            marks_cellLabels[mark_ind] = np.unique(cellMask[norm_MM[mark_ind]['x'], norm_MM[mark_ind]['y']])
            # Record which cells are under which ablation marks
            # marks_cell[key] = [np.where(np.unique(cellMask) == labeli)[0][0]
            #                    for labeli in labels_underMM[labels_underMM > 0]]
            marks_samplingArea[mark_ind] = len(np.where(cellMask[norm_MM[mark_ind]['x'], norm_MM[mark_ind]['y']]>0)[0])
            mark_area[mark_ind] = len(norm_MM[mark_ind]['x'])
            print('{}/{}'.format(i, len(norm_MM.keys())))

        nucl_fluo = {}
        cell_marks = {}
        cell_fluo = {}
        marks_fluo = {}
        marks_cell_overlap = {}
        # marks_area = {}
        marks_fluo_overlap = {}
        cell_area = {}
        marks_cell_overlap_indexes = {}

        for cell_ind in tqdm.tqdm(np.unique(cellMask)):
            #cell_ind is the object number given by CellProfiler
            if cell_ind > 0:
                # label
                # i = label
                nucl_fluo[str(cell_ind)] = []  # mean Raw intensity from third channel image
                cell_area[str(cell_ind)] = []  #number of pixels in each cell
                cell_marks[str(cell_ind)] = [] #marks index (or linear ion image pixel index) for each cell
                cell_fluo[str(cell_ind)] = [] #mean Raw intensity from second channel image
                marks_cell_overlap[str(cell_ind)] = [] #number of pixels shared by a cell and its touching ablation mark(s)
                # marks_area[str(cell_ind)] = [] #number of pixels in each a
                marks_fluo_overlap[str(cell_ind)] = [] #mean Raw intensity from second channel image of pixels shared between a cell and its touching ablation mark(s)
                marks_fluo[str(cell_ind)] = [] #mean Raw intensity from second channel image of pixels under touching ablation mark(s) of a cell
                marks_cell_overlap_indexes[str(cell_ind)] = [] #Indexes of touching ablation mark(s) of the cell
                cellMask_bw = cellMask == cell_ind

                # print(np.mean(fluo[cellMask == label]))
                for mark_ind, value in norm_MM.items():

                    if True in np.unique(cellMask_bw[norm_MM[mark_ind]['x'], norm_MM[mark_ind]['y']]):
                        # print(label, i, key)

                        # plt.show()
                        nucl_fluo[str(cell_ind)] = np.append(nucl_fluo[str(cell_ind)],
                                                             np.mean(fluo_nucl[cellMask == cell_ind]))
                        overlap_indices = np.where(cellMask_bw[norm_MM[mark_ind]['x'], norm_MM[mark_ind]['y']] == 1)[0]
                        marks_fluo[str(cell_ind)] = np.append(cell_fluo[str(cell_ind)],
                                                       np.mean(fluo[norm_MM[mark_ind]['x'], norm_MM[mark_ind]['y']]))  # TODO: remove np.fliplr
                        marks_fluo_overlap[str(cell_ind)] = np.append(marks_fluo_overlap[str(cell_ind)], np.mean(
                                                               fluo[norm_MM[mark_ind]['x'][overlap_indices],
                                                                    norm_MM[mark_ind]['y'][overlap_indices]]))
                        # Record which ablation marks are over which cells
                        cell_marks[str(cell_ind)] = np.append(cell_marks[str(cell_ind)], mark_ind)
                        cell_fluo[str(cell_ind)] = np.append(cell_fluo[str(cell_ind)],
                                                      np.mean(fluo[cellMask == cell_ind]))
                        marks_cell_overlap[str(cell_ind)] = np.append(marks_cell_overlap[str(cell_ind)],
                                                       len(overlap_indices))
                        marks_cell_overlap_indexes[str(cell_ind)] = np.append(marks_cell_overlap_indexes[str(cell_ind)],
                                                       mark_ind)
                        # marks_area[str(cell_ind)] = np.append(marks_area[str(cell_ind)], len(norm_MM[mark_ind]['x']))
                        cell_area[str(cell_ind)] = np.append(cell_area[str(cell_ind)], len(np.where(cellMask_bw == 1)[0]))
                        #
                        # plt.scatter(norm_MM[key]['y'], norm_MM[key]['x'], np.mean(fluo[cellMask == label]) / 200)

        cellMask_bw_all = cellMask > 0
        # plt.scatter(Y-np.min(Y),X-np.min(X))
        pmi = []  # Positive Mark Index
        mbi = []  # Mark brightfield intensity
        overLaps = []
        for i in tqdm.tqdm(range(np.shape(marksMask)[0])):
            status = 0
            cell_mark_OL = 0
            bi = []
            for j in range(np.shape(marksMask[i][0])[0]):
                # bi = np.append(bi, bf[
                #     int(marksMask[i][0][j] - np.min(coordX) + window), int(marksMask[i][1][j] - np.min(coordY) + window)])
                if cellMask_bw_all[
                    int(marksMask[i][0][j] - np.min(coordX) + window), int(
                                        marksMask[i][1][j] - np.min(coordY) + window)]:
                    status = 1
                    cell_mark_OL += 1
                    # if status == 1:
            pmi = np.append(pmi, status)
            # mbi = np.append(mbi, np.mean(bi))
            overLaps = np.append(overLaps, cell_mark_OL)

        np.save(Fname + 'marks_flitered_fluo.npy', [norm_MM, cell_marks, nucl_fluo, cell_fluo, marks_fluo, marks_cell_overlap,
                                                    mark_area, overlap_indices, marks_fluo_overlap, cell_area,
                                                    marks_cell_overlap_indexes, marks_cellLabels, marks_samplingArea, pmi, overLaps])
    else:
        norm_MM, cell_marks, nucl_fluo, cell_fluo, marks_fluo, marks_cell_overlap, mark_area, overlap_indices, marks_fluo_overlap, cell_area, \
        marks_cell_overlap_indexes, marks_cellLabels, marks_samplingArea, pmi, overLaps = np.load(Fname + 'marks_flitered_fluo.npy')

    if fetch_ann == 'online':
        if not os.path.exists(Fname + 'filter_results.npy'): #TODO uncomment this
            pmi = np.reshape(pmi, [int(np.sqrt(coordX.shape[0])), int(np.sqrt(coordX.shape[0]))]).ravel()
            # mbi = np.reshape(mbi, [int(np.sqrt(coordX.shape[0])), int(np.sqrt(coordX.shape[0]))]).ravel()

            def err_func(params, pmi, an_vec255):
                threshold = params
                corr = correlation(pmi, an_vec255 > threshold)
                return corr

            def running_mean(l, N):
                sum = 0
                result = list(0 for x in l)

                for i in range(0, N):
                    sum = sum + l[i]
                    result[i] = sum / (i + 1)

                for i in range(N, len(l)):
                    sum = sum - l[i - N] + l[i]
                    result[i] = sum / N

                return result
            # for CD in np.linspace(0, 1, 50):
            # pmi = np.flipud(d.isotope_images('C13H8O6', '+H')[0]).ravel() #TODO small hack, remove it after

            sm = smau.SMInstance()
            d = sm.dataset(ds_name)

            if filter == 'correlation':
                fdr = 0.5
                results1 = sm.msm_scores([d], d.annotations(fdr, database='HMDB'), db_name='HMDB').T
                results2 = sm.msm_scores([d], d.annotations(fdr, database='ChEBI'), db_name='ChEBI').T
                results3 = sm.msm_scores([d], d.annotations(fdr, database='LIPID_MAPS'), db_name='ChEBI').T
                results4 = sm.msm_scores([d], d.annotations(fdr, database='SwissLipids'), db_name='ChEBI').T
                results = pd.concat([results1, results2, results3, results4]).drop_duplicates()
                filter_results = []
                for i in tqdm.tqdm(range(results.shape[0])):
                    row = results.reset_index().iloc[i,:]
                    images = d.isotope_images(row.formula, row.adduct)
                    an_vec255 = tf_obj(images[0]).ravel()  # TODO: np.fliplr
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
                    filter_results = np.append(filter_results, best_corr)
                    # print('Molecule: {}, Corr. dist. = '.format(row.formula) + str(best_corr)[:5])

            elif filter == 'mean':
                fdr = 0.2
                # tol_fact = 0.2 #tolerence factor: number of std the mean of on-sample pixels have to be superior to the
                pmi_on = [i == 1.0 for i in pmi]
                pmi_off = [i == 0.0 for i in pmi]
                # mean of off_samples pixels for the annotation image to be considered as on_sample, higher means more severe
                results1 = sm.msm_scores([d], d.annotations(fdr, database='HMDB'), db_name='HMDB').T
                results2 = sm.msm_scores([d], d.annotations(fdr, database='ChEBI'), db_name='ChEBI').T
                results3 = sm.msm_scores([d], d.annotations(fdr, database='LIPID_MAPS'), db_name='LIPID_MAPS').T
                results4 = sm.msm_scores([d], d.annotations(fdr, database='SwissLipids'), db_name='SwissLipids').T
                results = pd.concat([results1, results2, results3, results4]).drop_duplicates()
                filter_results = []
                for i in range(results.shape[0]):
                    row = results.reset_index().iloc[i,:]
                    images = d.isotope_images(row.formula, row.adduct)
                    an_vec255 = tf_obj(images[0]).ravel()  # TODO: np.fliplr
                    result = np.mean(an_vec255[pmi_on]) > np.mean(an_vec255[pmi_off]) + (tol_fact*np.std(an_vec255[pmi_off]))
                    if result == True:
                        print(row.formula)
                    filter_results = np.append(filter_results, result)
                print('{} annotations passed the filter'.format(np.unique(filter_results, return_counts=True)[1][1]))

            np.save(Fname + 'filter_results.npy', [filter_results, pmi])
        else:
            filter_results, pmi = np.load(Fname + 'filter_results.npy')

        CD_fname = Fname + 'CD={}/'.format(CDs[0])
        Tol_fname = Fname + 'tol_fact={}/'.format(tol_fact)
        if not os.path.exists(Tol_fname + 'tsne_inputs_Tol=' + str(tol_fact) + '.npy'):

                tsne_data = {}
                tsne_data_norm = {}
                tsne_wholeCell_fluo = {}
                tsne_sf = []
                tsne_adduct = []
                tsne_cellMarks_fluo = {}
                tsne_cellMarks_fluo_norm = {}
                tsne_overlap = {}
                tsne_cellArea = {}

                if not np.shape(np.where(filter_results < CDs[0])[0])[0] == 0 or 1.0 in filter_results:

                    for key, value in cell_marks.items():
                        if not np.shape(cell_fluo[key])[0] == 0:
                            tsne_data[key] = []
                            tsne_data_norm[key] = []
                            tsne_cellArea[key] = np.mean([int(i) for i in cell_area[key]])
                            tsne_overlap[key] = np.mean([int(i) for i in marks_cell_overlap[key]])
                            tsne_wholeCell_fluo[key] = np.mean([int(i) for i in cell_fluo[key]])
                            tsne_cellMarks_fluo[key] = np.mean([int(i) for i in marks_fluo[key]])
                            tsne_cellMarks_fluo_norm[key] = np.mean([int(i) for i in marks_fluo_overlap[key]])

                    if filter == 'correlation':
                        if not os.path.exists(CD_fname):
                            os.makedirs(CD_fname)
                    elif filter == 'mean':
                        if not os.path.exists(Tol_fname):
                            os.makedirs(Tol_fname)

                    for i, result in enumerate(filter_results):
                        if result < CDs[0] and filter == 'correlation' or filter == 'mean' and result == 1.0:
                            sf = results.reset_index().as_matrix()[i, 0]
                            adduct = results.reset_index().as_matrix()[i, 1]
                            an_vec255 = tf_obj(d.isotope_images(sf, adduct)[0]).ravel()
                            for key, value in cell_marks.items():
                                if not np.shape(cell_fluo[key])[0] == 0:

                                    tsne_data[key] = np.append(tsne_data[key], np.mean(an_vec255[[int(k) for k in cell_marks[key]]]))
                                    if norm_method == 'sum':
                                        tsne_data_norm[key] = np.append(tsne_data_norm[key], np.sum(
                                            an_vec255[[int(k) for k in cell_marks[key]]] / [int(k) for k in
                                                                                            marks_cell_overlap[key]]))
                                    if norm_method == 'mean':
                                        tsne_data_norm[key] = np.append(tsne_data_norm[key], np.mean(
                                            an_vec255[[int(k) for k in cell_marks[key]]] / [int(k) for k in
                                                                                            marks_cell_overlap[key]]))
                            tsne_sf = np.append(tsne_sf, sf)
                            tsne_adduct = np.append(tsne_adduct, adduct)
                    key = np.array(list(tsne_data.keys()))[0]
                    print(np.shape(tsne_data[key]))

                    if np.shape(tsne_data[key])[0] > 0:
                        tsne_input = np.array([val for key, val in tsne_data.items()])
                        tsne_input_norm = np.array([val for key, val in tsne_data_norm.items()])
                        overlapCellMarks_input = np.array([val for key, val in tsne_overlap.items()])
                        fluoCell_input = np.array([val for key, val in tsne_wholeCell_fluo.items()])
                        fluoMarks_input = np.array([val for key, val in tsne_cellMarks_fluo.items()])
                        fluoMarks_input_norm = np.array([val for key, val in tsne_cellMarks_fluo_norm.items()])
                        keys_input = np.array([key for key, val in tsne_data.items()])
                        cell_area_input = np.array([val for key, val in tsne_cellArea.items()])

                        ObjN = np.array([int(i) for i in keys_input])
                        MOLdf = pd.concat([pd.DataFrame({'ObjectNumber_lu': ObjN,
                                                           'fluoCellMean_lu': fluoCell_input,
                                                           'fluoMarksMean_lu': fluoMarks_input_norm,
                                                           'CellAreaPixel_lu': cell_area_input,
                                                           'OverlapCellMarks_lu': overlapCellMarks_input}),
                                             pd.DataFrame(data=tsne_input, columns=tsne_sf)],
                                            axis=1)

                        if filter == 'correlation':
                            np.save(CD_fname + 'MOLdata_CD=' + str(CD) + '_nAnnot = ' + str(np.shape(tsne_data[key])[0]) + ' .npy',
                                    [tsne_input, fluoCell_input, fluoMarks_input, keys_input, tsne_sf, tsne_adduct, tsne_data,
                                     filter_results, CDs, pmi, tsne_input_norm, fluoMarks_input_norm, overlapCellMarks_input, cell_area_input])

                            MOLdf.set_index('ObjectNumber_lu').to_csv(CD_fname + 'MOLallData.csv')
                            MOLdf[[i for i in itertools.chain(['ObjectNumber_lu'], tsne_sf[:])]].\
                                set_index('ObjectNumber_lu').to_csv(CD_fname + 'MOLonlyData.csv')

                        elif filter == 'mean':
                            np.save(Tol_fname + 'MOLdata_Tol=' + str(tol_fact) + '_nAnnot = ' + str(
                                np.shape(tsne_data[key])[0]) + ' .npy',
                                    [tsne_input, fluoCell_input, fluoMarks_input, keys_input, tsne_sf, tsne_adduct,
                                     tsne_data,
                                     filter_results, CDs, pmi, tsne_input_norm, fluoMarks_input_norm,
                                     overlapCellMarks_input, cell_area_input])

                            MOLdf.set_index('ObjectNumber_lu').to_csv(Tol_fname + 'MOLallData.csv')
                            MOLdf[[i for i in itertools.chain(['ObjectNumber_lu'], tsne_sf[:])]]. \
                                set_index('ObjectNumber_lu').to_csv(Tol_fname + 'MOLonlyData.csv')

    if fetch_ann == 'offline':

        CD_fname = Fname + '/offline/'
        if not os.path.exists(CD_fname):
            os.makedirs(CD_fname)

        df_im0 = pd.read_hdf(hdf5_path)
        df_im = df_im0[df_im0['ds_name'] == ds_name]

        tsne_data = {}
        tsne_data_norm = {}
        tsne_wholeCell_fluo = {}
        tsne_nucl_fluo = {}
        tsne_sf = []
        tsne_adduct = []
        tsne_cellMarks_fluo = {}
        tsne_cellMarks_fluo_norm = {}
        tsne_overlap = {}
        tsne_cellArea = {}
        MIA = {} #'mark_indexed_area',
        MISA = {} # 'marks_indexed_SamplingArea',
        MSR = {} # 'marks_sampling_ratio',
        MCOI = {} #'marks_cell_overlap_int',
        MCORWA = {} # 'marks_cell_overlap_ratio_whole_area',
        MCORSA = {} # 'marks_cell_overlap_ratio_sampling_area'
        nMarks = {}
        # mark_indexes_d = {} # indexes of ablation marks touching the indexed cell (cell_ind)
        # marks_intensities_d = {} # raw intensities of ablation marks for the indexed annotation
        # mark_indexed_area_d = {}# number of pixels (from microscopy present within each
        # # the indexed ablation marks
        # marks_indexed_SamplingArea_d = {} # the number of pixels from the indexed
        # # ablation marks in common with any other cells
        # marks_sampling_ratio_d = {}  # The proportion of the indexed ablation marks's pixels touching any cells
        # marks_cell_overlap_int_d = {}
        # marks_cell_overlap_ratio_whole_area_d = {}
        # marks_cell_overlap_ratio_sampling_area_d = {}

        for cell_ind, value in cell_marks.items():
            if not np.shape(cell_fluo[cell_ind])[0] == 0:
                tsne_data[cell_ind] = []
                tsne_data_norm[cell_ind] = []
                MIA[cell_ind] = []
                MISA[cell_ind] = []
                MSR[cell_ind] = []
                MCOI[cell_ind] = []
                MCORWA[cell_ind] = []
                MCORSA[cell_ind] = []
                nMarks[cell_ind] = []
                tsne_cellArea[cell_ind] = np.mean([int(i) for i in cell_area[cell_ind]])
                tsne_overlap[cell_ind] = np.mean([int(i) for i in marks_cell_overlap[cell_ind]])
                tsne_nucl_fluo[cell_ind] = np.mean([int(i) for i in nucl_fluo[cell_ind]])
                tsne_wholeCell_fluo[cell_ind] = np.mean([int(i) for i in cell_fluo[cell_ind]])
                tsne_cellMarks_fluo[cell_ind] = np.mean([int(i) for i in marks_fluo[cell_ind]])
                tsne_cellMarks_fluo_norm[cell_ind] = np.mean([int(i) for i in marks_fluo_overlap[cell_ind]])

        if not os.path.exists(CD_fname):
            os.makedirs(CD_fname)
        mark_areas_array = np.array([k for k in mark_area.values()])
        low_t = np.percentile(np.array([k for k in mark_area.values()]), 5)
        high_t = np.percentile(np.array([k for k in mark_area.values()]), 99)

        # for i, best_corr in enumerate(filter_results):
        #     if best_corr < CD:
        #         sf = results.reset_index().as_matrix()[i, 0]
        #         adduct = results.reset_index().as_matrix()[i, 1]
        #         an_vec255 = np.rot90(np.flipud(d.isotope_images(sf, adduct)[0]), 3).ravel()
        #
        # marks_df = pd.DataFrame(columns=['mark_indexed_area', 'marks_indexed_SamplingArea', 'marks_sampling_ratio',
        #                                  'marks_cell_overlap_int', 'marks_cell_overlap_ratio_whole_area', 'marks_cell_overlap_ratio_sampling_area'])
        # marks_df['mark_indexed_area'] = marks_df['mark_indexed_area'].astype(object)
        # marks_df['marks_indexed_SamplingArea'] = marks_df['marks_indexed_SamplingArea'].astype(object)
        # marks_df['marks_sampling_ratio'] = marks_df['marks_sampling_ratio'].astype(object)
        # marks_df['marks_cell_overlap_int'] = marks_df['marks_cell_overlap_int'].astype(object)
        # marks_df['marks_cell_overlap_ratio_whole_area'] = marks_df['marks_cell_overlap_ratio_whole_area'].astype(object)
        # marks_df['marks_cell_overlap_ratio_sampling_area'] = marks_df['marks_cell_overlap_ratio_sampling_area'].astype(object)

        for j in tqdm.tqdm(range(df_im.shape[0])):
            an_vec255 = tf_obj(df_im.iloc[j, 1]).ravel()
            sf = df_im.iloc[j, 2]
            adduct = '-H'
            for cell_ind, value in cell_marks.items():
                if not np.shape(cell_fluo[cell_ind])[0] == 0:

                    marks_indexes = [int(mark_ind) for mark_ind in cell_marks[cell_ind]] #indexes of ablation marks touching the indexed cell (cell_ind)
                    marks_intensities = an_vec255[marks_indexes] #raw intensities of ablation marks for the indexed annotation
                    mark_indexed_area = np.array([mark_area[mark_ind] for mark_ind in cell_marks[cell_ind]]) #number of pixels (from microscopy present within each
                    # the indexed ablation marks
                    marks_indexed_SamplingArea = np.array([marks_samplingArea[mark_ind] for mark_ind in cell_marks[cell_ind]]) #the number of pixels from the indexed
                    # ablation marks in common with any other cells
                    marks_sampling_ratio = marks_indexed_SamplingArea/mark_indexed_area #The proportion of the indexed ablation marks's pixels touching any cells
                    tsne_data[cell_ind] = np.append(tsne_data[cell_ind],
                                               np.mean(marks_intensities))
                    marks_cell_overlap_int = np.array([int(k) for k in marks_cell_overlap[cell_ind]])
                    marks_cell_overlap_ratio_whole_area = marks_cell_overlap_int/mark_indexed_area
                    marks_cell_overlap_ratio_sampling_area = marks_cell_overlap_int / marks_indexed_SamplingArea

                    # Critical step: set two treshold to discard ablation marks:
                    # 1 - The proportion of the AM touching the cell OI (default 30%)
                    # 1 - The proportion of the AM sampling any cell (default 30%)

                    marks_ind_filter = np.array(marks_cell_overlap_ratio_sampling_area>area_prop) & \
                                       np.array(marks_sampling_ratio>area_prop) & \
                                       np.array(low_t<mark_indexed_area) & \
                                       np.array(mark_indexed_area<high_t)

                    if norm_method == 'Mean_intensity':
                        # Mean Intensity of AM touching the cell OI
                        tsne_data_norm[cell_ind] = np.append(tsne_data_norm[cell_ind], np.mean(marks_intensities[marks_ind_filter]))

                    if norm_method == 'weighted_mean_sampling_area_MarkCell_overlap_int_power1,7':
                        # Mean Intensity of AM dvidied by its sampling area weighted by the number of pixels of the
                        # AM sampling area touching the cell OI
                        cell_ion_intensity = np.sum(
                            marks_intensities[marks_ind_filter] * marks_cell_overlap_int[marks_ind_filter]**1.7 /
                            marks_indexed_SamplingArea[marks_ind_filter]) / np.sum(
                            marks_cell_overlap_int[marks_ind_filter])  # /marks_sampling_ratio[marks_ind_filter])
                        tsne_data_norm[cell_ind] = np.append(tsne_data_norm[cell_ind], cell_ion_intensity)

                    if norm_method == 'weighted_mean_sampling_area_MarkCell_overlap_int':
                        # Mean Intensity of AM dvidied by its sampling area weighted by the number of pixels of the
                        # AM sampling area touching the cell OI
                        cell_ion_intensity = np.sum(
                            marks_intensities[marks_ind_filter] * marks_cell_overlap_int[marks_ind_filter] /
                            marks_indexed_SamplingArea[marks_ind_filter]) / np.sum(
                            marks_cell_overlap_int[marks_ind_filter])  # /marks_sampling_ratio[marks_ind_filter])
                        tsne_data_norm[cell_ind] = np.append(tsne_data_norm[cell_ind], cell_ion_intensity)

                    if norm_method == 'weighted_mean_sampling_area_MarkCell_overlap_ratio_sampling_area':
                        # Mean Intensity of AM dvidied by its sampling area weighted by the ration between
                        # AM sampling area touching the cell OI and its sampling area
                        cell_ion_intensity = np.sum(
                            marks_intensities[marks_ind_filter] * marks_cell_overlap_ratio_sampling_area[marks_ind_filter] /
                            marks_sampling_ratio[marks_ind_filter]) / np.sum(
                            marks_cell_overlap_ratio_sampling_area[marks_ind_filter])  # /marks_sampling_ratio[marks_ind_filter])
                        tsne_data_norm[cell_ind] = np.append(tsne_data_norm[cell_ind], cell_ion_intensity)

                    if norm_method == 'squared_weighted_mean_sampling_area_MarkCell_overlap_ratio_sampling_area':
                        # Mean Intensity of AM dvidied by its sampling area weighted by the ration between
                        # AM sampling area touching the cell OI and its sampling area
                        cell_ion_intensity = np.sum(
                            marks_intensities[marks_ind_filter] * marks_cell_overlap_ratio_sampling_area[marks_ind_filter]**2 /
                            marks_sampling_ratio[marks_ind_filter]) / np.sum(
                            marks_cell_overlap_ratio_sampling_area[marks_ind_filter]**2)  # /marks_sampling_ratio[marks_ind_filter])
                        tsne_data_norm[cell_ind] = np.append(tsne_data_norm[cell_ind], cell_ion_intensity)

                    if norm_method == 'mean_sampling_ratio':
                        # Mean Intensity of AM dvidied by the proportion of its area sampling a cell
                        cell_ion_intensity = np.mean(
                            marks_intensities[marks_ind_filter] / marks_sampling_ratio[marks_ind_filter])
                        tsne_data_norm[cell_ind] = np.append(tsne_data_norm[cell_ind], cell_ion_intensity)

                    if norm_method == 'mean_sampling_area':
                        # Mean Intensity of AM dvidied by its sampling area in um
                        cell_ion_intensity = np.mean(
                            marks_intensities[marks_ind_filter] / (marks_indexed_SamplingArea[marks_ind_filter]*pixel_size))
                        tsne_data_norm[cell_ind] = np.append(tsne_data_norm[cell_ind], cell_ion_intensity)

                    MIA[cell_ind] = np.append(MIA[cell_ind], np.mean(mark_indexed_area))
                    MISA[cell_ind] = np.append(MISA[cell_ind], np.mean(marks_indexed_SamplingArea))
                    MSR[cell_ind] = np.append(MSR[cell_ind], np.mean(marks_sampling_ratio))
                    MCOI[cell_ind] = np.append(MCOI[cell_ind], np.mean(marks_cell_overlap_int))
                    MCORWA[cell_ind] = np.append(MCORWA[cell_ind], np.mean(marks_cell_overlap_ratio_whole_area))
                    MCORSA[cell_ind] = np.append(MCORSA[cell_ind], np.mean(marks_cell_overlap_ratio_sampling_area))
                    nMarks[cell_ind] = np.append(nMarks[cell_ind], np.shape(marks_cell_overlap_ratio_sampling_area)[0])
            tsne_sf = np.append(tsne_sf, sf)
            tsne_adduct = np.append(tsne_adduct, adduct)
        key = np.array(list(tsne_data.keys()))[0]

        tsne_input = np.array([val for key, val in tsne_data.items()])
        tsne_input_norm = np.array([val for key, val in tsne_data_norm.items()])
        overlapCellMarks_input = np.array([val for key, val in tsne_overlap.items()])
        fluoNucl_input = np.array([val for key, val in tsne_nucl_fluo.items()])
        fluoCell_input = np.array([val for key, val in tsne_wholeCell_fluo.items()])
        fluoMarks_input = np.array([val for key, val in tsne_cellMarks_fluo.items()])
        fluoMarks_input_norm = np.array([val for key, val in tsne_cellMarks_fluo_norm.items()])
        keys_input = np.array([key for key, val in tsne_data.items()])
        cell_area_input = np.array([val for key, val in tsne_cellArea.items()])

        MIA_input = np.array([np.mean(val) for key, val in MIA.items()])
        MISA_input = np.array([np.mean(val) for key, val in MISA.items()])
        MSR_input = np.array([np.mean(val) for key, val in MSR.items()])
        MCOI_input = np.array([np.mean(val) for key, val in MCOI.items()])
        MCORWA_input = np.array([np.mean(val) for key, val in MCORWA.items()])
        MCORSA_input = np.array([np.mean(val) for key, val in MCORSA.items()])
        nMarks_input = np.array([np.mean(val) for key, val in nMarks.items()])

        ObjN = np.array([int(i) for i in keys_input])
        mask = ~np.any(np.isnan(tsne_input_norm), axis=1)
        MOLdf_final = pd.concat([pd.DataFrame({'ObjectNumber_lu': ObjN[mask],
                                               'fluoNuclMean_lu': fluoNucl_input[mask],
                                               'fluoCellMean_lu': fluoCell_input[mask],
                                               'fluoMarksMean_lu': fluoMarks_input_norm[mask],
                                               'CellAreaPixel_lu': cell_area_input[mask],
                                               'OverlapCellMarks_lu': overlapCellMarks_input[mask],
                                               'mark_indexed_area': MIA_input[mask],
                                               'marks_indexed_SamplingArea': MISA_input[mask],
                                               'marks_sampling_ratio': MSR_input[mask],
                                               'marks_cell_overlap_int': MCOI_input[mask],
                                               'marks_cell_overlap_ratio_whole_area': MCORWA_input[mask],
                                               'marks_cell_overlap_ratio_sampling_area': MCORSA_input[mask],
                                               'number_of_marks': nMarks_input[mask]}),
                                 pd.DataFrame(data=tsne_input_norm[mask], columns=tsne_sf)], axis=1)
        MOLdf_final.set_index('ObjectNumber_lu').to_csv(CD_fname + 'MOLallData.csv')
        MOLdf_final[[i for i in itertools.chain(['ObjectNumber_lu'], tsne_sf[:])]].set_index(
            'ObjectNumber_lu').to_csv(
            CD_fname + 'MOLonlyData.csv')

def defMORPHfeatures(MF):
    """Reads morphological features of interest from CellProfiler csv output, and save them as a a new csv.
    Documentation of the quantified features:
    # Main page: http://cellprofiler.org/manuals/current/
    # http://cellprofiler.org/manuals/current/MeasureObjectSizeShape.html
    # http://cellprofiler.org/manuals/current/MeasureObjectIntensity.html
    # http://cellprofiler.org/manuals/current/MeasureObjectNeighbors.html

    Args:
        MF (str): path to Main Folder.

    """

    MFA = MF + 'Analysis/'
    Fname = MFA + 'scAnalysis/Morphological_features/'
    if not os.path.exists(Fname):
        os.makedirs(Fname)
    feat_df = pd.read_csv(MFA + 'CellProfilerAnalysis/Cells.csv')

    features_OI = ['ObjectNumber',
                   'AreaShape_Area',
                   'AreaShape_Compactness',
                   'AreaShape_Eccentricity',
                   'AreaShape_EulerNumber',
                   'Intensity_IntegratedIntensity_OrigRed_highDR',
                   'Intensity_MADIntensity_OrigRed_highDR',
                   'Intensity_MaxIntensity_OrigRed_highDR',
                   'Intensity_MeanIntensity_OrigRed_highDR',
                   'Intensity_MedianIntensity_OrigRed_highDR',
                   'Intensity_MinIntensity_OrigRed_highDR',
                   'Intensity_StdIntensity_OrigRed_highDR',
                   'Location_Center_X',
                   'Location_Center_Y',
                   'Neighbors_FirstClosestDistance_Adjacent',
                   'Neighbors_NumberOfNeighbors_Adjacent',
                   'Neighbors_PercentTouching_Adjacent',
                   'Neighbors_SecondClosestDistance_Adjacent',
                   'Intensity_IntegratedIntensity_OrigBlue',
                   'Intensity_MADIntensity_OrigBlue',
                   'Intensity_MaxIntensity_OrigBlue',
                   'Intensity_MeanIntensity_OrigBlue',
                   'Intensity_MedianIntensity_OrigBlue',
                   'Intensity_MinIntensity_OrigBlue',
                   'Intensity_StdIntensity_OrigBlue',
                   ]

    features_OI_names = ['ObjectNumber',
                         'Area',
                         'Compactness',
                         'Eccentricity',
                         'EulerNumber',
                         'Intensity_SUM',
                         'Intensity_MAD',
                         'Intensity_MAX',
                         'Intensity_MEAN',
                         'Intensity_MEDIAN',
                         'Intensity_MIN',
                         'Intensity_STD',
                         'Location_X',
                         'Location_Y',
                         'FirstClosestDistance',
                         'NumberOfNeighbors',
                         'PercentTouching',
                         'SecondClosestDistance',
                         'Intensity_SUM_nucl',
                         'Intensity_MAD_nucl',
                         'Intensity_MAX_nucl',
                         'Intensity_MEAN_nucl',
                         'Intensity_MEDIAN_nucl',
                         'Intensity_MIN_nucl',
                         'Intensity_STD_nucl',
                         ]

    feat_df_selected = feat_df[features_OI]
    feat_df_selected.columns = features_OI_names
    feat_df_selected.set_index('ObjectNumber').to_csv(Fname + 'MORPHselectedData.csv')
    feat_df.set_index('ObjectNumber').to_csv(Fname + 'MORPHallData.csv')

def mergeMORPHnMOL4cyt(MF, CDs, fetch_ann='offline',  tol_fact=0.2, filter = 'mean'):
    """Merge molecular data from the cells analyzed with SpaceM with their morphological features measured
    by CellProfiler. The matching is done using the values from Objectnumber'.

    Args:
        MF (str): path to Main Folder.
        CDs (list): correlation distance tresholds used for filtering background annotation images, only used when
            filter is 'correlation'. Default value is 0.75.
        fetch_ann (str): method for fetching annotations:
            'online': (default) queries metaspace using the name of the .imzml data present in the MALDI input folder
                as dataset name,
            'offline': reads annotation images from a provided dataframe..
        tol_fact (int): tolerance factor to use for the filter 'mean'.
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
         tol_fact (float): tolerance factor to use for the filter 'mean'.

    """

    if fetch_ann == 'online' and filter == 'correlation':
        MOLcsv_p = MF + 'Analysis/scAnalysis/Molecular_features/CD={}/MOLallData.csv'.format(CDs[0])
    elif fetch_ann == 'online' and filter == 'mean':
        MOLcsv_p = MF + 'Analysis/scAnalysis/Molecular_features/tol_fact={}/MOLallData.csv'.format(tol_fact)
    if fetch_ann == 'offline':
        MOLcsv_p = MF + 'Analysis/scAnalysis/Molecular_features/offline/MOLallData.csv'
    MORPHcsv_p = MF + 'Analysis/scAnalysis/Morphological_features/MORPHselectedData.csv'

    MOLdf = pd.read_csv(MOLcsv_p)
    # MOLdf_log = pd.DataFrame(columns= np.array(MOLdf.columns),data=np.nan_to_num(np.log10(MOLdf)))
    MORPHdf = pd.read_csv(MORPHcsv_p)
    MORPHnMOL_df = pd.concat([MORPHdf.iloc[MOLdf.ObjectNumber_lu-1,:].reset_index(), MOLdf.reset_index()], axis=1)
    MORPHnMOL_df = MORPHnMOL_df.set_index('ObjectNumber').drop(['index'], axis=1)
    MORPHnMOL_df.to_csv(MF + 'Analysis/scAnalysis/MORPHnMOL.csv')
    # quick check
    sns.regplot(MORPHnMOL_df.CellAreaPixel_lu, MORPHnMOL_df.Area, scatter_kws={"s": 80})
    plt.axis('equal')
    plt.xlabel('Cell Area computed from custom Pipeline')
    plt.ylabel('Cell Area computed from CellProfiler')
    plt.savefig(MF + 'Analysis/scAnalysis/AREA_index_check.png', psi=100)
    plt.close('all')

    sns.regplot(MORPHnMOL_df.fluoCellMean_lu, MORPHnMOL_df.Intensity_MEAN, scatter_kws={"s": 80})
    # plt.axis('equal')
    plt.xlabel('Fluo intensity computed from custom Pipeline')
    plt.ylabel('Fluo intensity computed from CellProfiler')
    plt.savefig(MF + 'Analysis/scAnalysis/FLUO_index_check.png', psi=100)
    plt.close('all')

def mapAnn2microCells(MF, MFA, csv_p, tf_obj,
                      labelled_cells_path, ds_index=10, draw_AM=False,
                      coloring_field='NumberOfNeighbors',
                      clip_percentile=100, cmap=cm.jet, log10=False):
    """Create an image using the label image from CellProfiler where cells are colored based on their intensity for a
    given metabolite.

    Args:
        MF (str): path to the Main Folder.
        MFA (str): path to the Main Folder Analysis.
        csv_p (str): path to the csv containing the molecular and morphological features of the cells.
        tf_obj (function): Image transformation to apply on ion image for registration.
        labelled_cells_path (str): path to the label image from CellProfiler.
        ds_index (int): index of the dataset. Stored in the csv under the field 'ds_index'.
        draw_AM (bool): whether drawing the ablation marks colored with their metabolite intensity on top of the cells.
        coloring_field (str): field from which the intensity will be used to color the cells/ablation marks.
        clip_percentile (float): percentile value to clip the intensities (hot/cold spot removal). The data are clipped
            in both direction using that value (ex: a clip_percentile values of 2.5 will result in 95% of the value range)
        cmap (matplotlib.cm): colormap to use to color the cells.
        log10 (bool): whether log10 transform the intensities from the csv.

    Returns:
        color_mask (array): the resulting labeled image in which each pixel from each cells have their corresponding
            value from the given coloring_field (2D).

    """

    data_i = pd.read_csv(csv_p, sep='\s*,\s*', header=0, encoding='ascii', engine='python')
    if log10:
        data = np.log10(data_i+1)
        data.ObjectNumber_lu = data_i.ObjectNumber_lu
    else:
        data = data_i
    cell_mask = tiff.imread(labelled_cells_path)
    color_mask = np.zeros(cell_mask.shape)
    sf_intensity_nz = np.copy(data[coloring_field].as_matrix())
    sf_intensity_nz = np.clip(sf_intensity_nz, np.percentile(sf_intensity_nz, 100-clip_percentile), np.percentile(sf_intensity_nz, clip_percentile))
    sf_min = np.nanmin(sf_intensity_nz)
    sf_max = np.nanmax(sf_intensity_nz - sf_min)
    if ds_index == None:
        for ObjN in tqdm.tqdm(data.ObjectNumber_lu.as_matrix()):
            color_mask[cell_mask == ObjN] = (data[data.ObjectNumber_lu == ObjN][coloring_field] - sf_min) / sf_max
    else:
        for ObjN in tqdm.tqdm(data[data.id_rep == ds_index].ObjectNumber_lu.as_matrix()):
            color_mask[cell_mask == ObjN] =(data[data.ObjectNumber_lu == ObjN][data.id_rep == ds_index][coloring_field]- sf_min) / sf_max

    out1 = np.array(cmap(color_mask[100:-100,100:-100])*255, dtype=np.uint8)
    # out2 = np.array(color_mask[100:-100,100:-100]*255, dtype=np.uint8)
    tiff.imsave(MFA + 'Mapping/' + coloring_field + '_cmap.tif', out1[:,:,:-1])
    np.save(MFA + 'Mapping/' + coloring_field + '_cmap.npy', color_mask[100:-100,100:-100])

    # tiff.imsave(MFA + 'Mapping/' + coloring_field + '_gray.tif', color_mask[100:-100,100:-100])
    if draw_AM:
        marksMask = np.load(MFA + 'SURF/transformedMarksMask.npy')
        coordX, coordY = np.load(MFA + 'SURF/transformedMarks.npy')
        images = pd.read_hdf('C:/Users\Luca\Google Drive\A-Team\projects/1c\hepatocytes, DKFZ\datasets\Molecular images/2017-09-12-luca-mz-images.hdf5')
        os.chdir(MF + 'Input/MALDI/')
        ds_name = glob.glob('*.imzML')[0].replace('.imzML', '')
        ion_img = tf_obj(images[np.array(images.mol_formula == coloring_field) & np.array(images.ds_name == ds_name)].image.as_matrix()[0])
        window = 0
        sns.set_style("whitegrid", {'axes.grid': False})
        plt.figure()
        plt.switch_backend('TkAgg')
        plt.get_backend()
        mng = plt.get_current_fig_manager()
        mng.window.state('zoomed')
        plt.imshow(out1)
        dz = (ion_img-np.min(ion_img)) / (np.max(ion_img)-np.min(ion_img))
        colors = cm.viridis(dz.ravel())
        plt.show()
        plt.pause(0.05)
        for i in tqdm.tqdm(range(0, coordX.shape[0])):
            x = np.array(marksMask[i, 0] - np.min(coordX) + window).astype(np.int64)
            y = np.array(marksMask[i, 1] - np.min(coordY) + window).astype(np.int64)
            plt.scatter(y, x, 0.01, colors[i])
        plt.savefig(MFA + 'Mapping/' + coloring_field + '_AM.png', dpi=500)
        plt.close('all')

    return color_mask[100:-100, 100:-100]

def annotation2microscopyAblationMarks(MF, sf, adduct, clip_percentile, touch_cell_only, tf_obj):
    """Overlaps the segmented ablation marks on the merged fluorescence with bright microscopy and color them with their
    corresponding metabolite intensity.

    Args:
        MF (str): path to the Main Folder.
        sf (str): sum formula of the metabolite to use for coloring the ablation marks.
        adduct (str): adduct to consider.
        clip_percentile (float): percentile value to clip the intensities (hot/cold spot removal). The data are clipped
            in both direction using that value (ex: a clip_percentile values of 2.5 will result in 95% of the value range).
        touch_cell_only (bool): whether to show only the ablation marks which are touching the cells.
        tf_obj (function): Image transformation to apply on ion image for registration.

    """


    MFA = MF + 'Analysis/'
    img = plt.imread(MFA + 'CellProfilerAnalysis/Contour_cells_adjusted.png')
    cellMask = tiff.imread(MFA + 'CellProfilerAnalysis/Labelled_cells.tif')
    cellMask_bw = cellMask>0
    marksMask = np.load(MFA + 'SURF/transformedMarksMask.npy')
    coordX, coordY = np.load(MFA + 'SURF/transformedMarks.npy')


    sm = smau.SMInstance()
    os.chdir(MF + 'Input/MALDI/')
    ds_name = glob.glob('*.imzML')[0].replace('.imzML', '')
    d = sm.dataset(ds_name)
    # results = sm.msm_scores([d], d.annotations(fdr, database='HMDB'), db_name='HMDB').T
    ion_img = np.log10(tf_obj(d.isotope_images(sf, adduct)[0])+1)
    # plt.imshow(ion_img, interpolation='none')
    window = 0
    sns.set_style("whitegrid", {'axes.grid': False})
    plt.figure()
    plt.switch_backend('TkAgg')
    plt.get_backend()
    mng = plt.get_current_fig_manager()
    mng.window.state('zoomed')
    plt.imshow(img)
    # ax.grid(False)
    dz = np.clip(ion_img, np.percentile(ion_img, clip_percentile), np.percentile(ion_img, 100 - clip_percentile)).ravel()
    colors = plt.cm.viridis((dz - dz.min()) / (dz.max() - dz.min()))
    plt.show()
    plt.pause(0.05)
    if touch_cell_only:
        for i in tqdm.tqdm(range(0, coordX.shape[0])):
            # norm_MM[str(i)] = {}
            x = np.array(marksMask[i, 0] - np.min(coordX) + window).astype(np.int64)
            y = np.array(marksMask[i, 1] - np.min(coordY) + window).astype(np.int64)
            if True in np.unique(cellMask_bw[x+window, y+window]):
                # print(i)
                plt.scatter(y, x, 1, colors[i], vmin=2.05)
    else:
        for i in tqdm.tqdm(range(0, coordX.shape[0])):
            # print(i)
            # norm_MM[str(i)] = {}
            x = np.array(marksMask[i, 0] - np.min(coordX) + window).astype(np.int64)
            y = np.array(marksMask[i, 1] - np.min(coordY) + window).astype(np.int64)
            plt.scatter(y, x, 0.01, colors[i], vmin=2.67)
            # plt.pause(0.05)

    if not os.path.exists(MFA + 'CellProfilerAnalysis/overlaps'):
        os.makedirs(MFA + 'CellProfilerAnalysis/overlaps')
    plt.savefig(MFA + 'CellProfilerAnalysis//overlaps/' + sf + '.png', dpi=1200, bbox_inches='tight', frameon='false')
    print('Overlap finished')
    plt.close('all')