import SpaceM.ImageFileManipulation.FIJIcalls as fc
import SpaceM.ImageFileManipulation.manipulations as manip
import SpaceM.Registration.AblationMarkFinder as amf
import SpaceM.Registration.ImageRegistration as ir
from SpaceM import WriteILIinput
import SpaceM.scAnalysis.Segmentation as scSg
import SpaceM.scAnalysis.scAnalysis as scSc
import SpaceM.scAnalysis.manifoldAnalysis as scMa
import numpy as np
import os, gc, matlab.engine

#Defines paths of main analysis directroy
MF = 'E:/Experiments/manual_cleaning_demo/' #MF: Main Folder
matrix = 'DHB' #or 'DAN'
CDs = [0.75] #correlation distances

def ion2fluoTF(ion_img):
    """Image transformation to apply on ion image for registration.

    Args:
        ion_img (array): the ion image to transform (2D).

    Returns:
        out (array): the transformed ion image.

    """
    # return ion_img.T  # --> TF1 HepaJune dataset batches FASTER
    # return np.fliplr(ion_img) #--> TF2 HepaJune dataset batch
    return np.flipud(ion_img)  # --> TF for HepaNov17

MFA = MF + 'Analysis/'
MFI = MF + 'Input/'
MFA_Spom = MFA + 'StitchedMicroscopy/postMALDI_FLR/'
MFA_Sprm = MFA + 'StitchedMicroscopy/preMALDI_FLR/'
MFA_gF = MFA + 'gridFit/'

#Generate subdirectory for analysis and for stitched microscopy
if not os.path.exists(MFA):
    os.makedirs(MFA)
    os.mkdir(MFA + 'StitchedMicroscopy')

# STEP2 -> process image preMALDI
# In Analysis/StitchedMicroscopy, make a new folder containing all individual images
# flipped left to right (required for FIJI stitching)
# create formatted coordinate text file of all images (also required for stitching)
if not os.path.exists(MFA_Sprm):
    os.makedirs(MFA_Sprm)
tif_files = manip.PixFliplr(MFI + 'Microscopy/preMALDI/', MFA_Sprm) #flip images left to right and copy them to destination fodler
fc.TileConfFormat(MFI + 'Microscopy/preMALDI/', MFA_Sprm, tif_files) #reformat metadafile from microscope into coordinate text file for FIJI stitching
gc.collect() #free the RAM
fc.callFIJIstitch(MFA_Sprm) #Call Fiji to perform stitching
# Generate a composite image with the two first channels
fc.callFIJImergeRedGray(base_path = MFA_Sprm,
                        red_filename = 'img_t1_z1_c2',
                        gray_filename = 'img_t1_z1_c1',
                        save_filename= 'Composite.png')
gc.collect()
print('Pre-MALDI Stitching finished')

#Perform the same for postMALDI microscopy dataset
if not os.path.exists(MFA_Spom):
    os.makedirs(MFA_Spom)
tif_files = manip.PixFliplr(MFI + 'Microscopy/postMALDI/', MFA_Spom)
fc.TileConfFormat(MFI + 'Microscopy/postMALDI/', MFA_Spom, tif_files)
gc.collect()
fc.callFIJIstitch(MFA_Spom)
print('Pre-MALDI Stitching finished')

# STEP3 -> detect ablation marks
if not os.path.exists(MFA_gF):
    os.makedirs(MFA_gF )
#reads X and Y coordinates of all images and store them into variables, for both pre- and post-MALDI
[POST_picXcoord, POST_picYcoord] = fc.readTileConfReg(MFA_Spom)
[PRE_picXcoord, PRE_picYcoord] = fc.readTileConfReg(MFA_Sprm)
#Loop over all post-MALDI images and detect ablation marks. Put coordinates from individual images into the global
#stitched image dimensions using images coordinates
amf.MarkFinderFT(MFA, MFA_Spom, POST_picXcoord, POST_picYcoord, matrix)
#Should run the script pyCLIMS/manualCleaning.py to manually select points. Follow instruction in the script

#Fit a grid on the ablation marks in order to:
# - filter extra detection
# - sort the indexes to match the pixel indexes of MALDI ion images
amf.GridFit(MFA, MFI, optimization=False, manual_cleaning=True, MarkFinderFT=True) #Look in the function definition for details

# Quick preview of the ablation mark detections in 'ili:
# Crops the stitched post-MALDI around the ablation marks detections
# Generate ili .csv input files
if not os.path.exists(MFA_gF + 'marks_check/'):
    os.makedirs(MFA_gF + 'marks_check/')

#Crops input images around coordinates
#On brightfield
manip.crop2coords(MFA + 'gridFit/xye_clean2.npy',
                             MFA_Spom + 'img_t1_z1_c1', MFA_gF + 'marks_check/PHASE_crop_bin1x1_window100.png', window = 100)
#On second channel (in the case of DHB, autofluorescence is used to detect ablation marks)
manip.crop2coords(MFA + 'gridFit/xye_clean2.npy',
                             MFA_Spom + 'img_t1_z1_c3', MFA_gF + 'marks_check/FLUO_crop_bin1x1_window100.png', window = 100)

manip.crop2coords(MFA + 'gridFit/xye_clean2.npy',
                             MFA_Spom + 'img_t1_z1_c1', MFA_gF + 'marks_check/PHASE_crop_bin1x1.png', window = 0)
#write ili .csv input files
nbin = fc.imbin4ili(MFA_gF + 'marks_check/PHASE_crop_bin1x1.png', maxsize=50e6)
predata = WriteILIinput.preCSVdatagen(MFA + 'gridFit/xye_clean2.npy', radius=10, nbin=nbin, PlainFirst=True)
WriteILIinput.writeCSV(path=MFA_gF + 'marks_check/ablation_marks_checkDETECTIONS.csv', data=predata)

predata = WriteILIinput.preCSVdatagen(MFA + 'gridFit/xye_grid.npy', radius=10, nbin=nbin, PlainFirst=True)
WriteILIinput.writeCSV(path=MFA_gF + 'marks_check/ablation_marks_checkTHEORETICAL.csv', data=predata)

# STEP4 -> detect SURF features around fiducial penmarks on both pre- and post-MALDI microscopy dataset
if not os.path.exists(MFA + 'SURF/'):
    os.makedirs(MFA + 'SURF/')
ir.SURF(MFA, folder=MFA_Spom, picXcoord=POST_picXcoord, picYcoord=POST_picYcoord, prefix='post', int_treshold = [])
ir.SURF(MFA, folder=MFA_Sprm, picXcoord=PRE_picXcoord, picYcoord=PRE_picYcoord, prefix='pre', int_treshold = [])

ir.SURF_Alignment(MFA) #performs registration of SURF features coordinates from both datasets.
# Can run manualCleaning.py here

#Note: this can be run already on the ablation amrks coordinates in order to have an early preview of the data in 'ili
# later, the area of each ablation marks will be defined. The same function will then have to be run a second time to apply the
# transformation on the areas
ir.TransformMarks(MFA) #applies the registration function defined by ir.SURF_Alignment() on the ablation marks coordinates

# STEP5 -> Create input for ili of registered ablation marks coordinates over the pre-MALDI stitched microscopy dataset
# each ablation marks is colored by its corresponding ion image pixel intensity.
if not os.path.exists(MFA + 'ili/'):
    os.makedirs(MFA + 'ili/')
manip.crop2coords(MFA + 'SURF/transformedMarks.npy',
                  MFA_Sprm + 'Composite.png',
                  MFA + 'ili/FLUO_crop_bin1x1.png',
                  window = 0)
manip.crop2coords(MFA + 'SURF/transformedMarks.npy',
                  MFA_Sprm + 'Composite.png',
                  MFA + 'ili/FLUO_crop_bin1x1_window100.png',
                  window = 100)
gc.collect()
nbin = fc.imbin4ili(MFA + 'ili/FLUO_crop_bin1x1.png', maxsize=50e6)

WriteILIinput.annotationSM2CSV(MFA, MFI, fdr=0.5, nbin = nbin, radius=20, tf_obj=ion2fluoTF) #fetches ion images from METASPACE
# and write intensities in ili csv inut file

amf.marksSegmentedMask(MFA) #detects the coordinate of all pixels (segments) within each ablation marks
ir.TransformMarks(MFA) #applies registration function on the masks

#STEP6 -> Calls Cellprofiler to segment cells from the cropped pre-MALDI stitched microsocpy
if not os.path.exists(MFA + 'CellProfilerAnalysis/'):
    os.makedirs(MFA + 'CellProfilerAnalysis/')
CP_window = 100
manip.crop2coords4CP(MFA + 'SURF/transformedMarks.npy',
                             MFA_Sprm, MFA + 'CellProfilerAnalysis/', window = CP_window)
fc.callFIJImergeRedGray(base_path = MFA + 'CellProfilerAnalysis/',
                        red_filename = 'img_t1_z1_c2.tif',
                        gray_filename = 'img_t1_z1_c1.tif',
                        save_filename = 'Composite_window100_adjusted.png')
gc.collect()
eng = matlab.engine.start_matlab()
dummy = eng.imAdjQuantiles(0.01, MFA + 'CellProfilerAnalysis/img_t1_z1_c2.tif',
                           MFA + 'CellProfilerAnalysis/img_t1_z1_c2_adjusted.tif')
print('Start CellProfiler Anlalysis')
scSg.callCP(MFA)
print('Finished CellProfiler Anlalysis')
#Draw cells outlines to check the cell segmentation results
dummy = eng.cellOutlines(MFA + 'CellProfilerAnalysis/Composite_window100_adjusted.png',
                         CP_window,
                         MFA + 'CellProfilerAnalysis/Labelled_cells.tif',
                         MFA + 'CellProfilerAnalysis/Contour_cells_adjusted.png')
eng.quit()

CDs = [0.75]
#STEP7 -> Start the single cell analysis
#Objective here is to assign a single intensity value from each ion for each sampled cell
if not os.path.exists(MFA + 'scAnalysis/'):
    os.makedirs(MFA + 'scAnalysis/')
scSc.defMORPHfeatures(MF)
fetch_ann ='online'
filter = 'mean' # either 'mean' or 'correlation'
tol_fact = -0.2
# for offline fetching: need to specify annotations, datasets ids
scSc.defMOLfeatures(MF, tf_obj=ion2fluoTF, CDs=CDs, norm_method='weighted_mean_sampling_area_MarkCell_overlap_ratio_sampling_area',
                    fetch_ann=fetch_ann, tol_fact=tol_fact, filter = filter)
scSc.mergeMORPHnMOL4cyt(MF, CDs=CDs, fetch_ann=fetch_ann, tol_fact=tol_fact, filter = filter)
scMa.tSNEgen(MF, CDs=CDs, metric='chebyshev', fetch_ann=fetch_ann,
             tol_fact=tol_fact, filter = filter)


# #Generate tSNE on 1c data
# scMa.tSNEgen(MF, CDs=CDs, metric='chebyshev', n_iter=1, fetch_ann=fetch_ann)
# #Generate cyt3 (MATLAB package) input file (.csv)
# scMa.genCYTinput(MF, iter_OI=0)
# if not os.path.exists(MFA + 'Mapping/'):
#     os.makedirs(MFA + 'Mapping/')
# ColoringFields = ['ObjectNumber', 'Area', 'Compactness', 'Eccentricity',
#        'EulerNumber', 'Intensity_SUM', 'Intensity_MAD', 'Intensity_MAX',
#        'Intensity_MEAN', 'Intensity_MEDIAN', 'Intensity_MIN',
#        'Intensity_STD', 'Location_X', 'Location_Y', 'FirstClosestDistance',
#        'NumberOfNeighbors', 'PercentTouching', 'SecondClosestDistance',
#        'ObjectNumber_lu', 'CellAreaPixel_lu', 'OverlapCellMarks_lu',
#        'fluoCellMean_lu', 'fluoMarksMean_lu', 'C10H14N5O7P', 'C10H17N3O6S',
#        'C18H13NO4', 'C18H32O2', 'C18H34O2', 'C19H24N2O2', 'C20H17FO4S',
#        'C20H30O2', 'C20H32O2', 'C20H34O2', 'C21H37O6P', 'C21H39O6P',
#        'C21H39O7P', 'C21H41O6P', 'C21H41O7P', 'C21H44NO7P', 'C22H30N2O5',
#        'C22H32O2', 'C22H36N2O5S', 'C23H44NO7P', 'C23H46NO7P', 'C23H48NO7P',
#        'C25H44NO7P', 'C25H48NO7P', 'C25H50NO7P', 'C26H32O12', 'C26H34O12',
#        'C26H36O12', 'C35H67O8P', 'C36H34O17', 'C37H69O8P', 'C37H71O8P',
#        'C39H69O8P', 'C39H73O8P', 'C39H74NO8P', 'C39H76NO8P', 'C40H77O10P',
#        'C41H74NO8P', 'C41H76NO8P', 'C41H78NO8P', 'C42H78NO10P',
#        'C42H82NO6P', 'C43H79O13P', 'C43H81O13P', 'C44H74O13P2',
#        'C44H76NO10P', 'C44H78NO10P', 'C45H66O5', 'C45H68O5', 'C45H78NO8P',
#        'C47H68O5', 'C47H70O5', 'C47H82O16P2', 'C47H84O16P2', 'C6H11O8P',
#        'C10H13N5O9P2', 'C16H23INO', 'C17H13O7', 'C22H28N2O5', 'C23H39O6P',
#        'C23H39O7P', 'C23H41O6P', 'C23H45O7P', 'C26H35N7O2', 'C27H51O12P',
#        'C29H49O12P', 'C36H57O9', 'C36H68NO6P', 'C36H69O8P', 'C38H71O8P',
#        'C38H73O8P', 'C38H75N2O6P', 'C39H67N5O6', 'C39H71O8P', 'C40H71O8P',
#        'C40H72O19', 'C40H75O8P', 'C41H71O8P', 'C41H73O8P', 'C42H80O10P',
#        'C42H81NO10P', 'C43H73O8P', 'C43H82O13P', 'C43H83NO8P',
#        'C44H78NO9P', 'C45H86O13P', 'C46H66O6', 'C46H77NO10P',
#        'C46H79NO10P', 'C46H86NO11P', 'C46H89N2O6P', 'C47H86O13P',
#        'C47H88O13P', 'C48H70O6', 'C48H81NO10P', 'C48H83NO10P',
#        'C50H85NO10P', 'C55H72N4O6', 'C6H8N3O7P2']
# replicates = ['U21', 'U22', 'U23', 'F11', 'F13', 'F21', 'L12', 'L13', 'L31', 'T22', 'T23', 'T31']
# replicates_int = [0,  1,      2,    3,      4,    5,      6,     7,     8,     9,     10,    11]
# ColoringFields = ['PercentTouching', 'PC1', 'PC2', 'PC3']
# ColoringFields = ['C47H84O16P2', 'C18H34O2', 'C43H73O8P', 'C47H82O16P2']
# csv_p = 'C:/Users/Luca\Google Drive\A-Team\projects/1c\hepatocytes, DKFZ\datasets\CSV/all_annotations/tresh=0.3\mean_sampling_area/NS=3/all_DS_cutoff=5_combat_fluoNormalized.csv'
# for CF in ColoringFields:
#     print(ColoringFields)
#     scSc.mapAnn2microCells(MF, MFA, csv_p=csv_p, ds_index=10,
#                            labelled_cells_path=MF + 'Analysis/CellProfilerAnalysis/Labelled_cells.tif',
#                            coloring_field=CF, scale='linear', draw_AM=False,
#                            clip_percentile=98, tf_obj=ion2fluoTF)
#     dummy = eng.overlapCellsCMAP(MFA + 'CellProfilerAnalysis/Contour_cells_adjusted.png', MFA + 'Mapping/' + CF + '_cmap.tif',
#                                MFA + 'Mapping/'+ CF +'.png' )
#
# eng.quit()


