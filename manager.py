import spaceM
import numpy as np
import pandas as pd
import os, gc, matlab.engine

def getPaths():
    return pd.read_json(os.path.dirname(spaceM.__file__) + '/paths.json')

def stitchMicroscopy(tf, preMALDI=True, postMALDI=True, composite=False):

    if not os.path.exists(MF + 'Analysis/'):
        os.makedirs(MF + 'Analysis/')
        os.mkdir(MF + 'Analysis/StitchedMicroscopy/')

    if preMALDI:

        if not os.path.exists(MF + 'Analysis/StitchedMicroscopy/preMALDI_FLR/'):
            os.makedirs(MF + 'Analysis/StitchedMicroscopy/preMALDI_FLR/')

        tif_files = spaceM.ImageFileManipulation.manipulations.PixFliplr(
            tf,
            MF + 'Input/Microscopy/preMALDI/',
            MF + 'Analysis/StitchedMicroscopy/preMALDI_FLR/')
        spaceM.ImageFileManipulation.FIJIcalls.TileConfFormat(MF + 'Input/Microscopy/preMALDI/',
                                                              MF + 'Analysis/StitchedMicroscopy/preMALDI_FLR/',
                                                              tif_files)
        gc.collect()
        spaceM.ImageFileManipulation.FIJIcalls.callFIJIstitch(MF + 'Analysis/StitchedMicroscopy/preMALDI_FLR/')
        print('Pre-MALDI Stitching finished')

    if postMALDI:

        if not os.path.exists(MF + 'Analysis/StitchedMicroscopy/postMALDI_FLR/'):
            os.makedirs(MF + 'Analysis/StitchedMicroscopy/postMALDI_FLR/')

        tif_files = spaceM.ImageFileManipulation.manipulations.PixFliplr(
            tf,
            MF + 'Input/Microscopy/postMALDI/',
            MF + 'Analysis/StitchedMicroscopy/postMALDI_FLR/')
        spaceM.ImageFileManipulation.FIJIcalls.TileConfFormat(MF + 'Input/Microscopy/postMALDI/',
                                                              MF + 'Analysis/StitchedMicroscopy/postMALDI_FLR/',
                                                              tif_files)
        gc.collect()
        spaceM.ImageFileManipulation.FIJIcalls.callFIJIstitch(MF + 'Analysis/StitchedMicroscopy/postMALDI_FLR/')
        print('Pre-MALDI Stitching finished')

    if composite:
        spaceM.ImageFileManipulation.FIJIcalls.callFIJImergeRedGray(
            base_path=MF + 'Analysis/StitchedMicroscopy/preMALDI_FLR/',
            red_filename='img_t1_z1_c2',
            gray_filename='img_t1_z1_c1',
            save_filename='Composite.png')

def ablationMarksFinder():

    if not os.path.exists(MF + 'Analysis/gridFit/'):
        os.makedirs(MF + 'Analysis/gridFit/')

    spaceM.Registration.AblationMarkFinder.MarkFinderFT(MF)

def fiducialsFinder():

    if not os.path.exists(MF + 'Analysis/Fiducials/'):
        os.makedirs(MF + 'Analysis/Fiducials/')

    spaceM.Registration.ImageRegistration.penMarksFeatures(MF,  prefix='post', int_treshold=[])
    spaceM.Registration.ImageRegistration.penMarksFeatures(MF,  prefix='pre' , int_treshold=[])

def ablationMarksFilter(marks_check=True):

    spaceM.Registration.AblationMarkFinder.GridFit(MF)

    if marks_check:

        if not os.path.exists(MF + 'Analysis/gridFit/marks_check/'):
            os.makedirs(MF + 'Analysis/gridFit/marks_check/')

        spaceM.ImageFileManipulation.manipulations.crop2coords(
            MF + 'Analysis/gridFit/xye_clean2.npy',
            MF + 'Analysis/StitchedMicroscopy/postMALDI_FLR/img_t1_z1_c1',
            MF + 'Analysis/gridFit/marks_check/PHASE_crop_bin1x1_window100.png',
            window=100)
        spaceM.ImageFileManipulation.manipulations.crop2coords(
            MF + 'Analysis/gridFit/xye_clean2.npy',
            MF + 'Analysis/StitchedMicroscopy/postMALDI_FLR/img_t1_z1_c1',
            MF + 'Analysis/gridFit/marks_check/PHASE_crop_bin1x1.png',
            window=0)

        nbin = spaceM.ImageFileManipulation.FIJIcalls.imbin4ili(
            MF + 'Analysis/gridFit/marks_check/PHASE_crop_bin1x1.png',
            maxsize=50e6)
        predata = spaceM.WriteILIinput.preCSVdatagen(
            MF + 'Analysis/gridFit/xye_clean2.npy',
            radius=10,
            nbin=nbin,
            PlainFirst=True)
        spaceM.WriteILIinput.writeCSV(
            path=MF + 'Analysis/gridFit/marks_check/ablation_marks_checkDETECTIONS.csv',
            data=predata)

        predata = spaceM.WriteILIinput.preCSVdatagen(
            MF + 'Analysis/gridFit/xye_grid.npy',
            radius=10,
            nbin=nbin,
            PlainFirst=True)
        spaceM.WriteILIinput.writeCSV(
            path=MF + 'Analysis/gridFit/marks_check/ablation_marks_checkTHEORETICAL.csv',
            data=predata)

def registration(ili_only=False):

    if not os.path.exists(MF + 'Analysis/Fiducials/optimized_params.npy'):
        spaceM.Registration.ImageRegistration.fiducialsAlignment(MF + 'Analysis/')

    if ili_only==True:
        if not os.path.exists(MF + 'Analysis/Fiducials/transformedMarksMask.npy'):
            spaceM.Registration.AblationMarkFinder.marksSegmentedMask(MF + 'Analysis/')

    if not os.path.exists(MF + 'Analysis/Fiducials/transformedMarksMask.npy'):
        spaceM.Registration.ImageRegistration.TransformMarks(MF + 'Analysis/')

    if not os.path.exists(MF + 'Analysis/ili/'):
        os.makedirs(MF + 'Analysis/ili/')

    spaceM.ImageFileManipulation.manipulations.crop2coords(
        MF + 'Analysis/Fiducials/transformedMarks.npy',
        MF + 'Analysis/StitchedMicroscopy/preMALDI_FLR/Composite.png',
        MF + 'Analysis/ili/FLUO_crop_bin1x1.png',
                      window=0)
    spaceM.ImageFileManipulation.manipulations.crop2coords(
        MF + 'Analysis/Fiducials/transformedMarks.npy',
        MF + 'Analysis/StitchedMicroscopy/preMALDI_FLR/Composite.png',
        MF + 'Analysis/ili/FLUO_crop_bin1x1_window100.png',
                      window=100)
    gc.collect()
    nbin = spaceM.ImageFileManipulation.FIJIcalls.imbin4ili(MF + 'Analysis/ili/FLUO_crop_bin1x1.png', maxsize=50e6)
    spaceM.WriteILIinput.annotationSM2CSV(
        MF + 'Analysis/',
        MF + 'Input/',
        fdr=0.5,
        nbin=nbin,
        radius=20,
        tf_obj=ion2fluoTF)

def cellSegmentation():

    if not os.path.exists(MF + 'Analysis/CellProfilerAnalysis/'):
        os.makedirs(MF + 'Analysis/CellProfilerAnalysis/')

    CP_window = 100
    spaceM.ImageFileManipulation.manipulations.crop2coords4CP(
        MF + 'Analysis/Fiducials/transformedMarks.npy',
        MF + 'Analysis/StitchedMicroscopy/preMALDI_FLR/',
        MF + 'Analysis/CellProfilerAnalysis/',
        window=CP_window)
    spaceM.ImageFileManipulation.FIJIcalls.callFIJImergeRedGray(
        base_path=MF + 'Analysis/CellProfilerAnalysis/',
        red_filename='img_t1_z1_c2.tif',
        gray_filename='img_t1_z1_c1.tif',
        save_filename='Composite_window100_adjusted.png')
    gc.collect()

    eng = matlab.engine.start_matlab()
    dummy = eng.imAdjQuantiles(0.01,
                               MF + 'Analysis/CellProfilerAnalysis/img_t1_z1_c2.tif',
                               MF + 'Analysis/CellProfilerAnalysis/img_t1_z1_c2_adjusted.tif')

    print('Start CellProfiler Anlalysis')
    spaceM.scAnalysis.Segmentation.callCP(MF + 'Analysis/')
    print('Finished CellProfiler Anlalysis')

    dummy = eng.cellOutlines(MF + 'Analysis/CellProfilerAnalysis/Composite_window100_adjusted.png',
                             CP_window,
                             MF + 'Analysis/CellProfilerAnalysis/Labelled_cells.tif',
                             MF + 'Analysis/CellProfilerAnalysis/Contour_cells_adjusted.png')
    eng.quit()

def spatioMolecularMatrix(CDs=[0.75]):

    if not os.path.exists(MF + 'Analysis/scAnalysis/'):
        os.makedirs(MF + 'Analysis/scAnalysis/')

    spaceM.scAnalysis.scAnalysis.defMORPHfeatures(MF)
    fetch_ann = 'online' # either 'online' or 'offline'
    filter = 'mean'  # either 'mean' or 'correlation'
    tol_fact = -0.2

    spaceM.scAnalysis.scAnalysis.defMOLfeatures(
        MF,
        tf_obj=ion2fluoTF,
        CDs=CDs,
        norm_method='weighted_mean_sampling_area_MarkCell_overlap_ratio_sampling_area',
        fetch_ann=fetch_ann, tol_fact=tol_fact, filter=filter)

    spaceM.scAnalysis.scAnalysis.mergeMORPHnMOL4cyt(
        MF,
        CDs=CDs,
        fetch_ann=fetch_ann,
        tol_fact=tol_fact,
        filter=filter)

def ion2fluoTF(ion_img):
    """Image transformation to apply on ion image for registration.

    Args:
        ion_img (ndarray): the ion image to transform (2D).

    Returns:
        out (array): the transformed ion image.

    """
    # return ion_img.T  # --> TF1 HepaJune dataset batches FASTER
    # return np.fliplr(ion_img) #--> TF2 HepaJune dataset batch
    return np.flipud(ion_img)  # --> TF for HepaNov17

def img_tf(img):
    # return img #20171206_CoCulture
    return np.fliplr(img) #20171106_Hepa_Nov, 20170622_Hepa_June_DAN_Untreated_FA_LPS_TNFa

MF = getPaths()['MF'].as_matrix()[0]
stitchMicroscopy(tf=img_tf, postMALDI=False)
# ablationMarksFinder()
# fiducialsFinder()
# #Run GUI to curate ablation marks and fiducials
# ablationMarksFilter()
# registration()
# cellSegmentation()
# spatioMolecularMatrix()