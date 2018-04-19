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

# def flr(path_in, path_out):
#     tif_files = []
#     ind = 1
#     for item in os.listdir(path_in):
#         if item.endswith('.tif'):
#             print('Processing Image # {}'.format(item))
#             a = tiff.imread(path_in + item)
#             a = (a / float(np.max(a))) * 16385.0
#             phase_lr = np.fliplr(a[:, :])  # TODO:1->0
#             # fluo_lr = np.fliplr(a[0, :, :])  # TODO:0->1
#             b = np.zeros((np.shape(phase_lr)[0], np.shape(phase_lr)[1]), dtype=np.uint16)
#             b[:, :] = phase_lr
#             # b[1, :, :] = fluo_lr
#             tiff.imsave(path_out + item, b)
#
#             tif_files.append(item)
#     return tif_files

def crop2coords(coords_p, img_p, save_p, window):
    X, Y = np.load(coords_p)
    # X = np.int(X)#/0.73
    # Y = np.int(Y)#/0.73
    X = [int(x) for x in X]
    Y = [int(y) for y in Y]
    img = plt.imread(img_p)
    scipy.misc.imsave(save_p, img[np.min(X)-window:np.max(X)+window, np.min(Y)-window:np.max(Y)+window])

def crop2coords4CP(coords_p, imgF_p, saveF_p, window):
    X, Y = np.load(coords_p)
    # X = X  # /0.73
    # Y = Y  # /0.73
    X = [int(x) for x in X]
    Y = [int(y) for y in Y]
    for item in os.listdir(imgF_p):
        if item.startswith('img_t1_z1'):
            print(imgF_p + item)
            img = np.array(plt.imread(imgF_p + item), dtype=np.uint16)
            # print(np.max(img))
            # img = img/np.max(img) * 65535
            # tiff = TIFF.open(saveF_p + item + '.tif', mode='w')
            # # scipy.misc.imsave(saveF_p + item + '.tif',
            # #                   img[np.min(X) - window:np.max(X) + window, np.min(Y) - window:np.max(Y) + window])
            # tiff.write_image(img[np.min(X) - window:np.max(X) + window, np.min(Y) - window:np.max(Y) + window])
            tiff.imsave(saveF_p + item + '.tif',
                               img[np.min(X) - window:np.max(X) + window, np.min(Y) - window:np.max(Y) + window])


# import time
# import os
# import tifffile as tiff
# import numpy as np
# import multiprocessing
# import pyCLIMS.ImageFileManipulation.PixFlipLR as pflr
#
# file_path = 'D:/Experiments/20160810_HeLa_Stm_WT_mCHerry/code_dev/input_all/Input/Microscopy/preMALDI/'
#
# start = time.time()
# for item in os.listdir(file_path):
#     pflr.flr(file_path + item)
# elapsed = time.time() - start
# print elapsed
#
# start = time.time()
# pool = multiprocessing.Pool()
# pool.map(pflr.flr, list(file_path + item for item in os.listdir(file_path)))
# elapsed = time.time() - start
# print elapsed
#
# from multiprocessing import Pool
# def f(x):
#     return x*x
#
# with Pool() as p:
#     p = Pool()
#
# p.map(f, [1, 2, 3])
#
#
# from multiprocessing import Process
# import sys
#
# # pool of 10 workers
# processes = []
# for i in range(4):
#     processes.append(Process(target=f, args=[1,2,3]))
#
# for p in processes:
#     p.start()
#
# for p in processes:
#     p.join()