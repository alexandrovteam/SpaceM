from spaceM.Pipeline import getPath, stitchMicroscopy, \
    ablationMarksFinder_old, \
    fiducialsFinder, \
    ablationMarksFilter, \
    registration, \
    cellSegmentation, \
    spatioMolecularMatrix, \
    curator
import spaceM.CellProfilerProjectFiles.CP_upstream as CP_upstream
import numpy as np

import spaceM
import pandas as pd
import os, gc
from subprocess import call
import GUI_maldi_helper

def img_tf(img):
    """Image transformation to apply to the tile images prior to stitching
    Args:
        img (array): the image to transform (2D).
    Returns:
        out (array): the transformed ion image.
    """
    # return img #20171206_CoCulture
    return np.fliplr(img) #20171106_Hepa_Nov, 20170622_Hepa_June_DAN_Untreated_FA_LPS_TNFa

def ion2fluoTF(ion_img):
    """Image transformation to apply on ion image for registration.
    Args:
        ion_img (ndarray): the ion image to transform (2D).
    Returns:
        out (array): the transformed ion image.
    """
    # return ion_img
    # return ion_img.T  # --> TF1 HepaJune
    # return np.fliplr(ion_img) #--> TF2 HepaJune, 20171206_CoCulture\M5
    return np.flipud(ion_img) # --> TF for HepaNov17, 20180514_Coculture

def prepCP(MF):
    CP_upstream.hepatocytes(MF) #--> HepaJune, HepaNov datasets
    # CP_upstream.coculture(MF) #--> 20180514_Coculture

channels = ['gray', 'red', 'blue']

MF = getPath('MF')
print('Main folder: {}'.format(MF))

stitchMicroscopy(MF,
                 preMALDI=False,
                 postMALDI=True,
                 tf=img_tf,
                 merge_colors=['gray', 'red', 'blue'],
                 merge_filenames=['img_t1_z1_c1.tif', 'img_t2_z1_c1.tif', 'img_t3_z1_c1.tif'])

# curator(MF + 'Analysis/StitchedMicroscopy/postMALDI_FLR/img_t1_z1_c1')
# ablationMarksFinder_old(MF)
# fiducialsFinder(MF)

# curator(MF + 'Analysis/gridFit/ablation_marks_XY.npy', 'Ablation marks\nSelect ablation marks')
# curator(MF + 'Analysis/Fiducials/postXYpenmarks.npy', 'Post-MALDI\nRemove noisy detections')
# curator(MF + 'Analysis/Fiducials/preXYpenmarks.npy', 'Pre-MALDI\nRemove noisy detections')

ablationMarksFilter(MF)
registration(MF, tf_obj=ion2fluoTF, ili_fdr=0.2)

cellSegmentation(MF,
                 merge_colors=['gray', 'green', 'red'],
                 merge_filenames=['img_t2_z1_c1.tif', 'img_t3_z1_c1.tif', 'img_t4_z1_c1.tif'],
                 prepCP_fun=prepCP)

# curator(MF + 'Analysis/CellProfilerAnalysis/Labelled_cells.tiff', 'Label image, remove wrong segmentations')

spatioMolecularMatrix(MF, tf_obj=ion2fluoTF, fetch_ann = 'offline')

