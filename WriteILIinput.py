# from pyimzml.ImzMLParser import ImzMLParser, getionimage
from sm_analytics_python.sm_annotation_utils import sm_annotation_utils as smau
import numpy as np
import csv
import glob, os

def preCSVdatagen(xy_p, radius, nbin, PlainFirst):
    """Format the data before generating the csv input for ili'.

    Args:
        xy_p (str): path to the X and Y coordiantes of ablation marks .npy file.
        radius (int): displayed radius of the marks in ili'.
        nbin (int): bin factor used to bin the image for ili'.
        PlainFirst (bool): intensity values of each datapoints are equal to 1. Used to visualize the ablation mark
            coordinates on the postMALDI brighfield in ili'.
    Returns:
        data (list): formatted data (2D).

    """
    X, Y = np.load(xy_p)
    Xs = X /( nbin)  # todo check relevance of Y <-> X
    Ys = Y /( nbin)
    Ys = Ys - np.min(Ys)
    Xs = Xs - np.min(Xs)
    Rs = np.ones(np.shape(Xs)) * radius
    data = []
    data.append(list(np.append('Num', list(range(np.shape(Xs.ravel())[0])))))
    data.append(list(np.append('X', Ys.ravel())))
    data.append(list(np.append('Y', Xs.ravel())))
    data.append(list(np.append('Z', np.zeros(np.shape(Xs.ravel())))))
    data.append(list(np.append('R', Rs.ravel())))
    if PlainFirst:
        data.append(list(np.append('Flat', np.ones(np.shape(Xs.ravel())))))
    return data

def writeCSV(path, data):
    """Writes the formatted data in a csv file.

    Args:
        path (str): str of the csv file to write.
        data (list): data to write (2D).

    """
    data_csv = list(zip(*data))
    with open(path, 'w') as testfile:
        cw = csv.writer(testfile)
        for i in range(np.shape(data_csv)[0]):
            cw.writerow(data_csv[i])

def annotationSM2CSV(MFA, MFI, fdr, nbin, radius, tf_obj):
    """Fetches annotation images from METASPACE (http://metaspace2020.eu/#/about) and writes intensity values of
    each ablation marks in a csv input for ili' (https://ili.embl.de/). Used to visualize the ion signal on the
    preMALDI microsocpy after registration and validate the geometric transform to apply to the ion image.

    Args:
        MFA (str): path to Main Folder Analysis.
        MFI (str): path to Main Folder Input.
        fdr (float): fdr threshold value can only be 0.05, 0.1, 0.2 and 0.5.
        nbin (int): bin factor used to bin the image for ili'.
        radius (int): displayed radius of the marks in ili'.
        tf_obj (function): Image transformation to apply on ion image for registration.

    """
    def CSVdatagen(data, results, d):
        """Writes intensity values of each ablation marks in a csv input for ili'.
        TODO finish documentation
        Args:
            data (list): data to populate with ion intensities (2D).
            results (): .
            d (): .

        Returns:
            data (list): data to write in csv input for ili.

        """
        ind = 0
        for i, row in enumerate(results.reset_index().itertuples()):
            images = d.isotope_images(row.formula, row.adduct)
            print(row.formula)
            data.append(list(np.append(row[1], tf_obj(images[0]).ravel())))
            ind += 1
        return data

    sm = smau.SMInstance()
    os.chdir(MFI + 'MALDI/')
    ds_name = glob.glob('*.imzML')[0].replace('.imzML', '')
    d = sm.dataset(ds_name)
    results = sm.msm_scores([d], d.annotations(database='HMDB-v2.5', fdr=fdr)).T

    predata = preCSVdatagen(MFA + 'Fiducials/transformedMarks.npy', radius, nbin, PlainFirst=False)
    data_csv = CSVdatagen(predata, results, d)
    writeCSV(path = MFA + '/ili/sm_annotation_detections.csv', data = data_csv)

