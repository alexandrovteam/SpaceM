import os
import tifffile as tiff
import numpy as np
import matplotlib.pyplot as plt
import scipy.misc
import tqdm

def PixFliplr(file_path, MFAspm):
    """Transform images before stitching. This transformation depends on the microscope, the software and the canera used.

    Args:
        file_path (str): path of directory containing the tiled images to transform.
        MFAspm (str): path of directory where the transformed images will be saved.

    Returns:
        tif_files (array): array containing the names of the transformed images.

    """
    # if not os.path.exists(dir_fliplr):
    #     os.makedirs(dir_fliplr)
    tif_files = []
    ind = 1
    for item in tqdm.tqdm(os.listdir(file_path)):
        if item.endswith('.tif'):
            a = tiff.imread(file_path + item)
            # a = (a/float(np.max(a)))*16385.0
            if np.shape(a.shape)[0] == 2:
                b = np.zeros((np.max(a.shape), np.max(a.shape)),
                             dtype=np.uint16)  # lazy hack --> restricts to squared images
                b[:, :] = np.fliplr(a[:, :])
            else:
                n_chan = a.shape[0]
                b = np.zeros((n_chan, np.max(a.shape), np.max(a.shape)), dtype=np.uint16) # lazy hack --> restricts to squared images
                for i in range(n_chan):
                    b[i, :, :] = np.fliplr(a[i, :, :])
            tiff.imsave(MFAspm + item, b)
            tif_files.append(item)
            # print('Processed Image # {}'.format(ind))
            ind  = ind+1
    return tif_files

def crop2coords(coords_p, img_p, save_p, window):
    X, Y = np.load(coords_p)
    X = [int(x) for x in X]
    Y = [int(y) for y in Y]
    img = plt.imread(img_p)
    scipy.misc.imsave(save_p, img[np.min(X)-window:np.max(X)+window, np.min(Y)-window:np.max(Y)+window])

def crop2coords4CP(coords_p, imgF_p, saveF_p, window):
    X, Y = np.load(coords_p)
    X = [int(x) for x in X]
    Y = [int(y) for y in Y]
    for item in os.listdir(imgF_p):
        if item.startswith('img_t1_z1'):
            print(imgF_p + item)
            img = np.array(plt.imread(imgF_p + item), dtype=np.uint16)
            tiff.imsave(saveF_p + item + '.tif',
                               img[np.min(X) - window:np.max(X) + window, np.min(Y) - window:np.max(Y) + window])
