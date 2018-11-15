from subprocess import call
import matplotlib.pyplot as plt
import numpy as np
import tqdm
from scipy import ndimage

def callCP(MFA, cp_p, cppipe_p):
    """Call CellProfiler (http://cellprofiler.org/) to perform cell segmentation. CellProfiler segmentation pipeline
    is in the spaceM folder with the '.cppipe' extension.

     Args:
         MFA (str): path to Main Folder Analysis.
         cp_p (str): path to CellProfiler path
         cppipe_p (str): path to CellProfiler pipeline file

     """
    # CP headless info https://github.com/CellProfiler/CellProfiler/wiki/Adapting-CellProfiler-to-a-LIMS-environment
    call([cp_p,
          '-r',
          '-c',
          '-p',
          cppipe_p,
          '-o',
          MFA + 'CellProfilerAnalysis\\', '--file-list',
          MFA + 'CellProfilerAnalysis\input_files.txt'])

def  cellOutlines(FluoBrightfield_p, fluo_window, label_p, save_p, clusters=[], cluster_col=[], labels_OI = []):
    """Visualize the cell segmentation results from CellProfiler by drawing a black outline around the estimated cell
    boundaries.

     Args:
         FluoBrightfield_p (str): path to image to draw the cells outlines on.
         fluo_window (int): number of pixels surrounding the frame of interest
         label_p (str): path to the label image created by CellProfiler
         save_p (str): path to the generated image with cells outlines

     """
    if fluo_window > 0 :
        labelI = plt.imread(label_p)[fluo_window:-fluo_window]
        fluoI = plt.imread(FluoBrightfield_p)[fluo_window:-fluo_window]
    else:
        labelI = plt.imread(label_p)
        fluoI = plt.imread(FluoBrightfield_p)
    values = np.unique(labelI)
    perimAll = np.zeros(np.shape(labelI))
    struct = ndimage.generate_binary_structure(2, 1)
    FIC = fluoI#[fluo_window:-fluo_window]
    if np.shape(labels_OI)[0] > 0:
        label_list = labels_OI
    else:
        label_list = np.unique(labelI)

    for seed in tqdm.tqdm(values):
        BW = (labelI==seed)*1
        if seed in label_list and seed >0:
            perim = BW - ndimage.binary_erosion(BW, structure=struct, iterations=1).astype(BW.dtype)
            perim = ndimage.binary_dilation(perim, structure=struct, iterations=1).astype(BW.dtype)
            # perimAll = perimAll + (BW-erode)
            if np.shape(clusters)[0] > 0:
                if seed in clusters[0]:
                    color = cluster_col[clusters[1][clusters[0] == seed][0]]
                else:
                    color = [0,0,0]
                FIC[perim == 1, :] = color
            else:
                color = [0,0,0]
                FIC[perim == 1, :] = color
    plt.imshow(FIC)


    PAC = perimAll#[fluo_window:-fluo_window]
    PAC_d = ndimage.binary_dilation(PAC, structure=struct, iterations=1).astype(BW.dtype)

    CC = FIC*np.dstack([np.invert(PAC_d.astype('bool'))] * 3)
    plt.imsave(save_p, CC)

def cellDistribution_MALDI(MF):
    """Maps the distribution of the cells over the sampled area by MALDI as a binary matrix. Can also be called an On/Off
        sample mask, where pixels with a value of 1 are off sample (the corresponding ablation mark of that MALDI pixel
        does not overlap with a cell) and a value of 1 are ON sample (there is overlap).

     Args:
         MF (str): path to the Main Folder.

     """
    MFA = MF + 'Analysis/'
    cellMask = tiff.imread(MFA + 'CellProfilerAnalysis/Labelled_cells.tif')
    marksMask = np.load(MFA + 'Fiducials/transformedMarksMask.npy')
    coordX, coordY = np.load(MFA + 'Fiducials/transformedMarks.npy')
    window = 100

    cellMask_bw_all = cellMask > 0
    pmi = []  # Positive Mark Index
    overLaps = []
    for i in tqdm.tqdm(range(np.shape(marksMask)[0])):
        status = 0
        cell_mark_OL = 0
        bi = []
        for j in range(np.shape(marksMask[i][0])[0]):
            if cellMask_bw_all[
                int(marksMask[i][0][j] - np.min(coordX) + window), int(
                    marksMask[i][1][j] - np.min(coordY) + window)]:
                status = 1
                cell_mark_OL += 1
                # if status == 1:
        pmi = np.append(pmi, status)
        overLaps = np.append(overLaps, cell_mark_OL)

    np.save(MFA + 'CellProfilerAnalysis/cellDistribution_MALDI.npy', [pmi, overLaps])


