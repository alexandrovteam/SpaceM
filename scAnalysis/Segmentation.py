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
    file = open(MFA + 'CellProfilerAnalysis\input_files.txt', 'w')
    file.write(MFA + 'CellProfilerAnalysis\img_t1_z1_c1.tif\n')
    file.write(MFA + 'CellProfilerAnalysis\img_t1_z1_c2.tif\n')
    file.write(MFA + 'CellProfilerAnalysis\img_t1_z1_c2_adjusted.tif\n')
    file.write(MFA + 'CellProfilerAnalysis\img_t1_z1_c3_adjusted.tif')
    file.close()
    call([cp_p,
          '-r',
          '-c',
          '-p',
          cppipe_p,
          '-o',
          MFA + 'CellProfilerAnalysis\\', '--file-list',
          MFA + 'CellProfilerAnalysis\input_files.txt'])

def  cellOutlines(FluoBrightfield_p, fluo_window, label_p, save_p):
    """Visualize the cell segmentation results from CellProfiler by drawing a black outline around the estimated cell
    boundaries.

     Args:
         FluoBrightfield_p (str): path to image to draw the cells outlines on.
         fluo_window (int): number of pixels surrounding the frame of interest
         label_p (str): path to the label image created by CellProfiler
         save_p (str): path to the generated image with cells outlines

     """
    labelI = plt.imread(label_p)
    fluoI = plt.imread(FluoBrightfield_p)
    values = np.unique(labelI)
    perimAll = np.zeros(np.shape(labelI))
    struct = ndimage.generate_binary_structure(2, 1)

    for seed in tqdm.tqdm(values):
        BW = (labelI==seed)*1
        if seed in labelI and seed >0:
            erode = ndimage.binary_erosion(BW, structure=struct, iterations=1).astype(BW.dtype)
            perimAll = perimAll + (BW-erode)

    PAC = perimAll[fluo_window:-fluo_window]
    PAC_d = ndimage.binary_dilation(PAC, structure=struct, iterations=1).astype(BW.dtype)
    FIC = fluoI[fluo_window:-fluo_window]
    CC = FIC*np.dstack([np.invert(PAC_d.astype('bool'))] * 3)
    plt.imsave(save_p, CC)
