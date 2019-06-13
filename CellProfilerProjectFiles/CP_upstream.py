import spaceM
import matplotlib.pyplot as plt
import numpy as np
import tifffile as tif
import scipy.ndimage as scim
from skimage.morphology import ball

def scale(input):
    """Scale array between 0 and 1"""
    return (input - np.min(input)) / ((np.max(input) - np.min(input)))
def contrast(arr, min, max=1.0):
    """Clip array between min and max values"""
    return np.clip(arr, np.percentile(arr, min*100), np.percentile(arr, max*100))

def hepatocytes(MF):

    spaceM.ImageFileManipulation.manipulations.imAdjQuantiles(pc=[1],
                                                              im_p=MF + 'Analysis/CellProfilerAnalysis/img_t1_z1_c3.tif',
                                                              adj_p=MF + 'Analysis/CellProfilerAnalysis/img_t1_z1_c3_adjusted.tif')

    spaceM.ImageFileManipulation.manipulations.imAdjQuantiles(pc=[1],
                                                              im_p=MF + 'Analysis/CellProfilerAnalysis/img_t1_z1_c2.tif',
                                                              adj_p=MF + 'Analysis/CellProfilerAnalysis/img_t1_z1_c2_adjusted.tif')

    file = open(MF + '/Analysis/CellProfilerAnalysis\input_files.txt', 'w')
    file.write(MF + '/Analysis/CellProfilerAnalysis\img_t1_z1_c1.tif\n')
    file.write(MF + '/Analysis/CellProfilerAnalysis\img_t1_z1_c2.tif\n')
    file.write(MF + '/Analysis/CellProfilerAnalysis\img_t1_z1_c2_adjusted.tif\n')
    file.write(MF + '/Analysis/CellProfilerAnalysis\img_t1_z1_c3_adjusted.tif')
    file.close()

def coculture(MF):

    GFP_clip = spaceM.ImageFileManipulation.manipulations.imAdjQuantiles(pc=[92,98],
                                                              im_p=MF + 'Analysis/CellProfilerAnalysis/img_t3_z1_c1.tif',
                                                              adj_p=[])

    mCherry_clip = spaceM.ImageFileManipulation.manipulations.imAdjQuantiles(pc=[50,98],
                                                              im_p=MF + 'Analysis/CellProfilerAnalysis/img_t4_z1_c1.tif',
                                                              adj_p=[])

    DAPI_clip = spaceM.ImageFileManipulation.manipulations.imAdjQuantiles(pc=[90,97.5],
                                                              im_p=MF + 'Analysis/CellProfilerAnalysis/img_t5_z1_c1.tif',
                                                              adj_p=[])

    max_zstack = scale(np.max(np.dstack((GFP_clip, mCherry_clip, DAPI_clip)), axis=2))*65535
    tif.imsave(MF + 'Analysis/CellProfilerAnalysis/DAPI-mCherry-GFP-MaxStack.tif', max_zstack.astype('uint16'))



    Phase_clip = spaceM.ImageFileManipulation.manipulations.imAdjQuantiles(pc=[1, 96],
                                                              im_p=MF + 'Analysis/CellProfilerAnalysis/img_t2_z1_c1.tif',
                                                              adj_p=[])*-1
    # plt.imshow(Phase_clip, cmap='gray')
    DAPI_clip2 = spaceM.ImageFileManipulation.manipulations.imAdjQuantiles(pc=[5],
                                                                           im_p=MF + 'Analysis/CellProfilerAnalysis/img_t5_z1_c1.tif',
                                                                           adj_p=[])
    s = ball(50)
    h = int((s.shape[1] + 1) / 2)
    s = s[:h, :, :].sum(axis=0)
    s = (255 * (s - s.min())) / (s.max() - s.min())
    Phase_clip_RB = scim.white_tophat(Phase_clip, structure=s)
    DAPI_clip_RB = scim.white_tophat(DAPI_clip2, structure=s)

    max_zstack2 = scale(np.max(np.dstack((DAPI_clip_RB, Phase_clip_RB)), axis=2)) * 65535
    tif.imsave(MF + 'Analysis/CellProfilerAnalysis/DAPI-Phase-MaxStack.tif', max_zstack2.astype('uint16'))

    file = open(MF + 'Analysis/CellProfilerAnalysis\input_files.txt', 'w')
    file.write(MF + 'Analysis/CellProfilerAnalysis/DAPI-mCherry-GFP-MaxStack.tif\n')
    file.write(MF + 'Analysis/CellProfilerAnalysis/DAPI-Phase-MaxStack.tif\n')
    file.write(MF + 'Analysis/CellProfilerAnalysis/img_t2_z1_c1.tif\n')
    file.write(MF + 'Analysis/CellProfilerAnalysis/img_t3_z1_c1.tif\n')
    file.write(MF + 'Analysis/CellProfilerAnalysis/img_t4_z1_c1.tif\n')
    file.write(MF + 'Analysis/CellProfilerAnalysis/img_t5_z1_c1.tif\n')
    file.close()

def coculture_NIH(MF):


    from skimage import data, color
    from skimage.transform import hough_circle, hough_circle_peaks
    from skimage.feature import canny
    from skimage.draw import circle_perimeter
    from skimage.util import img_as_ubyte

    GFP_clip = spaceM.ImageFileManipulation.manipulations.imAdjQuantiles(pc=[0,99.9],
                                                              im_p=MF + 'Analysis/CellProfilerAnalysis/img_t3_z1_c1.tif',
                                                              adj_p=[])
    # plt.imshow(scale(GFP_clip*-1), cmap='gray')
    GFP_inv = scale(GFP_clip*-1) * 65535
    GFP_clip = scale(GFP_clip) * 65535

    tif.imsave(MF + 'Analysis/CellProfilerAnalysis/GFP_inverted.tif', GFP_inv.astype('uint16'))
    tif.imsave(MF + 'Analysis/CellProfilerAnalysis/GFP_clip.tif', GFP_clip.astype('uint16'))



    hough_radii = np.arange(5, 20, 2)
    hough_res = hough_circle(GFP_clip, hough_radii)
    accums, cx, cy, radii = hough_circle_peaks(hough_res, hough_radii)
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10, 4))
    for center_y, center_x, radius in zip(cy, cx, radii):
        plt.scatter(center_y, center_x, radius, c=[0,1,0])





    max_zstack2 = scale(np.max(np.dstack((DAPI_clip_RB, Phase_clip_RB)), axis=2)) * 65535
    tif.imsave(MF + 'Analysis/CellProfilerAnalysis/DAPI-Phase-MaxStack.tif', max_zstack2.astype('uint16'))

    file = open(MF + 'Analysis/CellProfilerAnalysis\input_files.txt', 'w')
    file.write(MF + 'Analysis/CellProfilerAnalysis/DAPI-mCherry-GFP-MaxStack.tif\n')
    file.write(MF + 'Analysis/CellProfilerAnalysis/DAPI-Phase-MaxStack.tif\n')
    file.write(MF + 'Analysis/CellProfilerAnalysis/img_t2_z1_c1.tif\n')
    file.write(MF + 'Analysis/CellProfilerAnalysis/img_t3_z1_c1.tif\n')
    file.write(MF + 'Analysis/CellProfilerAnalysis/img_t4_z1_c1.tif\n')
    file.write(MF + 'Analysis/CellProfilerAnalysis/img_t5_z1_c1.tif\n')
    file.close()

def coculture_MH(MF):

    spaceM.ImageFileManipulation.manipulations.imAdjQuantiles(pc=[1],
                                                              im_p=MF + 'Analysis/CellProfilerAnalysis\img_t1_z1_c3.tif',
                                                              adj_p=MF + 'Analysis/CellProfilerAnalysis\img_t1_z1_c3_adjusted.tif')
    spaceM.ImageFileManipulation.manipulations.imAdjQuantiles(pc=[1],
                                                              im_p=MF + 'Analysis/CellProfilerAnalysis\img_t1_z1_c0.tif',
                                                              adj_p=MF + 'Analysis/CellProfilerAnalysis\img_t1_z1_c0_adjusted.tif')

    #Make Probability map in ILASTIK from img_t1_z1_c3_adjusted.tif and save as prob.tif in the same CellProfilerAnalysis folder

    bf = tif.imread(MF + 'Analysis/CellProfilerAnalysis\img_t1_z1_c0.tif')
    prob = tif.imread(MF + 'Analysis/CellProfilerAnalysis\prob.tif')
    ctfr = tif.imread(MF + 'Analysis/CellProfilerAnalysis\img_t1_z1_c3_adjusted.tif')

    bf = bf.astype('int32')
    dx = scim.sobel(bf, 1)  # horizontal derivative
    dy = scim.sobel(bf, 0)  # vertical derivative
    mag = np.hypot(dx, dy)

    segm_in = scale(scale(ctfr) * scale(prob) + scale(mag)) * 255
    tif.imsave(MF + 'Analysis/CellProfilerAnalysis\segmentation.tif', segm_in.astype('uint8'))

    file = open(MF + 'Analysis/CellProfilerAnalysis/input_files.txt', 'w')
    file.write(MF + 'Analysis/CellProfilerAnalysis/img_t1_z1_c1.tif\n')
    file.write(MF + 'Analysis/CellProfilerAnalysis/img_t1_z1_c2.tif\n')
    file.write(MF + 'Analysis/CellProfilerAnalysis/img_t1_z1_c3.tif\n')
    file.write(MF + 'Analysis/CellProfilerAnalysis/segmentation.tif\n')
    file.close()

def BrdU(stitch_folder = r'Y:\rappez\20190523_BrdU_Hoesch_Hepa\FI1\stitch/'):

    spaceM.ImageFileManipulation.manipulations.imAdjQuantiles(pc=[0.99],
                                                              im_p=stitch_folder + 'img_t1_z1_c1.tif',
                                                              adj_p=stitch_folder + 'img_t1_z1_c1_adj.tif')

    spaceM.ImageFileManipulation.manipulations.imAdjQuantiles(pc=[0.99],
                                                              im_p=stitch_folder + 'img_t2_z1_c1.tif',
                                                              adj_p=stitch_folder + 'img_t2_z1_c1_adj.tif')

    file = open(stitch_folder + 'input_files.txt', 'w')
    file.write(stitch_folder + 'img_t1_z1_c1.tif\n')
    file.write(stitch_folder + 'img_t2_z1_c1.tif\n')
    file.write(stitch_folder + 'img_t1_z1_c1_adj.tif\n')
    file.write(stitch_folder + 'img_t2_z1_c1_adj.tif')
    file.close()


