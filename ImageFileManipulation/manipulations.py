import os
import tifffile as tiff
import numpy as np
import matplotlib.pyplot as plt
import scipy.misc
import tqdm

def PixFliplr(tf, file_path, MFAspm):
    """Transform images before stitching. This transformation depends on the microscope, the software and the camera used.

    Args:
        tf (fun): transformation function (np.fliplr/np.flupud/.T) toi apply to the tile images prior to stitching.
        file_path (str): path of directory containing the tiled images to transform.
        MFAspm (str): path of directory where the transformed images will be saved.

    Returns:
        tif_files (array): array containing the names of the transformed images.

    """
    # if not os.path.exists(dir_fliplr):
    #     os.makedirs(dir_fliplr)
    tif_files = []
    ind = 1
    for item in os.listdir(file_path):
        if item.endswith('.tif'):
            a = tiff.imread(file_path + item)
            # a = (a/float(np.max(a)))*16385.0
            if np.shape(a.shape)[0] == 2:
                n_chan=1
                b = np.zeros((np.max(a.shape), np.max(a.shape)),
                             dtype=np.uint16)  # lazy hack --> restricts to squared images
                # b.shape = 1, 1, 1, a.shape[0], a.shape[1], 1
                b[:, :] = tf(a[:, :])
            else:
                n_chan = a.shape[0]
                b = np.zeros((n_chan, np.max(a.shape), np.max(a.shape)), dtype=np.uint16) # lazy hack --> restricts to squared images
                # b.shape =  n_chan, a.shape[1], a.shape[2], 1, 1 # dimensions in XYCZT order
                for i in range(n_chan):
                    b[i, :, :] = tf(a[i, :, :])

            tiff.imsave(MFAspm + item,  b)
            # plt.imsave(arr=b, fname=MFAspm + item)
            # fi.write_multipage(b, MFAspm + item)
            # io.imsave(MFAspm + item, b)
            tif_files.append(item)
            # print('Processed Image # {}'.format(ind))
            ind  = ind+1
    return tif_files

def crop2coords(coords_p, img_p, save_p, window):
    X, Y = np.load(coords_p)
    X = [int(x) for x in X]
    Y = [int(y) for y in Y]
    if img_p.split('/')[-1].split('.') == 'tif':
        img = tiff.imread(img_p)
    if len(img_p.split('/')[-1].split('.')) == 1:
        img = tiff.imread(img_p)
    else:
        img = plt.imread(img_p)
    tiff.imsave(save_p, img[np.min(X)-window:np.max(X)+window, np.min(Y)-window:np.max(Y)+window])
    if window == 0:
        scipy.misc.imsave(save_p, img[np.min(X)-window:np.max(X)+window, np.min(Y)-window:np.max(Y)+window])

def crop2coords4CP(coords_p, imgF_p, saveF_p, window):
    X, Y = np.load(coords_p)
    X = [int(x) for x in X]
    Y = [int(y) for y in Y]
    for item in os.listdir(imgF_p):
        if item.startswith('img_t'):
            print(imgF_p + item)
            img = np.array(tiff.imread(imgF_p + item), dtype=np.uint16)
            tiff.imsave(saveF_p + item + '.tif',
                               img[np.min(X) - window:np.max(X) + window, np.min(Y) - window:np.max(Y) + window])

def imAdjQuantiles(pc, im_p, adj_p):
    """Increase contrast of an image by clipping it at its 'pc' lowest and 'pc' highest percentile values.

    Args:
        pc (list): list percentile values to clip. If the list contains only one value, the clip will be symmetrical.
            If two values are provided, the bottom clipping is performed using a percentile value equal to the value
            at the first position in the list and the top clipping using the value at the second. The list must have a
            maximum of two values.
        im_p (str): path of the image to clip
        adj_p (str or list): save path of the clipped image. If an empty list is provided, the function returns the
            adjusted image variable.

    """

    def scale(input):
        """Scale array between 0 and 1"""
        return (input - np.min(input)) / ((np.max(input) - np.min(input)))

    pc_low = pc_high = 0

    if np.shape(pc)[0] == 1:
        pc_low = pc[0]
        pc_high = 100 - pc[0]
    elif np.shape(pc)[0] == 2:
        pc_low = pc[0]
        pc_high = pc[1]
    else:
        print('Clipping percentile format is wrong, no clipping is performed \non image {}'.format(im_p))

    img = scale(plt.imread(im_p))
    low_in = np.percentile(img, pc_low)
    high_in = np.percentile(img, pc_high)
    adjusted = scale(np.clip(img, low_in, high_in)) * 255
    # plt.imshow(adjusted)
    if adj_p == []:
        return adjusted
    else:
        tiff.imsave(adj_p, adjusted.astype('uint8'))

