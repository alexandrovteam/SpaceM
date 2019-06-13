# from pyimzml.ImzMLParser import ImzMLParser, getionimage
# from sm_analytics_python.sm_annotation_utils import sm_annotation_utils as smau
from metaspace import sm_annotation_utils as smau
import numpy as np
import csv
import glob, os, tqdm
import pandas as pd

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

def annotationSM2CSV(MFA, MFI, fdr, nbin, radius, tf_obj, db='HMDB-v4'):
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

    # config = {
    # 'graphql_url': 'http://staging.metaspace2020.eu/graphql',
    # 'moldb_url': 'http://staging.metaspace2020.eu/mol_db/v1',
    # 'jwt': None}
    sm = smau.SMInstance()
    sm.login(email='luca.rappez@embl.de', password='Zeppar12')

    os.chdir(MFI + 'MALDI/')
    ds_name = glob.glob('*.imzML')[0].replace('.imzML', '')
    d = sm.dataset(ds_name)
    results = sm.msm_scores([d], d.annotations(database=db, fdr=fdr), db_name=db).T

    predata = preCSVdatagen(MFA + 'Fiducials/transformedMarks.npy', radius, nbin, PlainFirst=False)
    data_csv = CSVdatagen(predata, results, d)
    writeCSV(path = MFA + '/ili/sm_annotation_detections.csv', data = data_csv)

def annotationSM2CSV_offline(MF,
                     tf_obj,
                     hdf5_path=r'F:\Google Drive\A-Team\projects\1c\hepatocytes_40samples, DKFZ\datasets/',
                     on_sample_list_path=r"F:\Google Drive\A-Team\projects\1c\hepatocytes_40samples, DKFZ\KATJAnMANUAL_ON_sample_annotations.csv"):

    MF = r'F:\Experiments\20171106_Hepa_Nov_ANALYSIS_PAPER\F3/'
    os.chdir(MF + 'Input/MALDI/')
    imzml_name = glob.glob('*.imzML')[0]
    ds_name = imzml_name.replace('.imzML', '')

    if os.path.isdir(hdf5_path):
        df_im0 = pd.concat([pd.read_hdf(p) for p in glob.glob(hdf5_path + '*.hdf5')])
    else:
        df_im0 = pd.read_hdf(hdf5_path)

    df_im = df_im0[df_im0['ds_name'] == ds_name].reset_index()
    on_mol_df = pd.read_csv(on_sample_list_path)
    Xs, Ys = np.load(MF + 'Analysis/Fiducials/transformedMarks.npy')
    Ys = Ys - np.min(Ys)
    Xs = Xs - np.min(Xs)

    ili_df = pd.DataFrame()
    ili_df['Num'] = list(range(len(Xs)))
    ili_df['X'] = Ys
    ili_df['Y'] = Xs
    ili_df['Z'] = np.ones(len(Xs)) * 0
    ili_df['R'] = np.ones(len(Xs)) * 20

    for i in tqdm.tqdm(df_im.index):
        mol_name = '{}, {}'.format(df_im.loc[i, 'mol_formula'], df_im.loc[i, 'adduct'])
        ili_df[mol_name] = tf_obj(df_im.loc[i, 'image']).ravel()

    ili_df.to_csv(MF + 'Analysis/ili/offline_on_sample.csv', index=False)





