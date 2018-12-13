from skimage import transform as tf
from skimage import filters
from skimage import exposure
from scipy import ndimage
from scipy.optimize import basinhopping
from scipy import spatial
import os, gc
import numpy as np
import matplotlib.pyplot as plt
import tifffile as tiff
import spaceM.ImageFileManipulation.FIJIcalls as fc

def penMarksFeatures(MF, prefix, whole_image=True):
    """Obtain coordinates of the pixels at the edge of the penmarks from the tile frames in both pre- and post-MALDI
    microscopy datasets using matlab implementation of the SURF algorithm.

    Args:
        MF (str): path to Main Folder.
        prefix (str): either 'pre' or 'post' for pre- or post-MALDI dataset, respectively.
        whole_image (bool): whether or not perform the fiducial detection on the stitched image. If False, performs
            on the tiled images, uses less RAM but much slower.
    """

    def fiducialFinder(im_p):
        # prefix= 'post'
        # im_p = MF + 'Analysis/StitchedMicroscopy/' + prefix + 'MALDI_FLR/' + 'img_t1_z1_c1'
        # im_p = 'E:\Experiments\TNFa_2.3_SELECTED\Analysis\StitchedMicroscopy\preMALDI_FLR\img_XY030.tif'

        im = tiff.imread(im_p)
        if len(np.shape(im)) > 2:
            im = im[0,:,:]
        # im = np.log2(im+1)
        val = filters.threshold_otsu(im, nbins=65536)
        val = np.mean(im) - np.std(im)/2

        hist, bins_center = exposure.histogram(im)

        plt.plot(bins_center, hist, lw=2)
        plt.axvline(val, color='w', ls='--', label='Threshold')
        plt.yscale('log')
        plt.xlabel('Pixel intensities')
        plt.ylabel('Log10(counts)')
        # plt.title('Mean- std/2')
        plt.legend()
        plt.savefig(MF + 'Analysis/Fiducials/' + prefix + '_histogram.png')
        plt.close('all')

        BW = np.zeros(np.shape(im))
        BW[im < val] = 1
        BW[im < 100] = 0
        im = []
        # opening
        struct2 = ndimage.generate_binary_structure(2,1)
        iteration = 10
        rec_o1 = ndimage.binary_erosion(BW, structure=struct2, iterations=iteration).astype(BW.dtype)
        BW = []
        rec_o2 = ndimage.binary_dilation(rec_o1, structure=struct2, iterations=iteration).astype(rec_o1.dtype)
        rec_o1 = []
        rec_c1 = ndimage.binary_dilation(rec_o2, structure=struct2, iterations=iteration).astype(rec_o2.dtype)
        rec_o2 = []
        rec_c2 = ndimage.binary_dilation(rec_c1, structure=struct2, iterations=iteration).astype(rec_c1.dtype)
        gc.collect()
        edge = filters.sobel(rec_c2)
        rec_o2 = []
        # plt.imshow(edge)
        x,y = np.where(edge > 0)
        return x,y

    folder = MF + 'Analysis/StitchedMicroscopy/' + prefix + 'MALDI_FLR/'
    if whole_image:
        X,Y = fiducialFinder(im_p=folder + 'img_t1_z1_c1')
    else:
        [picXcoord, picYcoord] = fc.readTileConfReg(folder)
        X = []
        Y = []
        for item in os.listdir(folder):
            if item.endswith('.tif'):
                x_coord, y_coord = fiducialFinder(folder + item)
                picInd = int(item[len(item) - 7:len(item) - 4])
                for i in range(len(x_coord)):
                    if x_coord[i] < (1608 - 1608*0.1) and y_coord[i] < (1608 - 1608*0.1):
                        xScaled = x_coord[i] + picXcoord[picInd - 1]
                        yScaled = y_coord[i] + picYcoord[picInd - 1]
                        X = np.append(X, xScaled)
                        Y = np.append(Y, yScaled)
                print(item)
    np.save(MF + 'Analysis/Fiducials/' + prefix + 'XYpenmarks.npy', [X, Y])

    plt.figure()
    plt.scatter(X, Y, 1)
    plt.xlabel('X dimension', fontsize=20)
    plt.ylabel('Y dimension', fontsize=20)
    plt.title('Fiducials detection ' + prefix +'MALDI', fontsize=25)
    plt.axis('equal')
    plt.savefig(MF + 'Analysis/Fiducials/' + prefix + 'CHECK.png', dpi = 500)
    plt.close('all')

def transform(postX, postY, transX, transY, rot):
    """Coordinate transform function.

    Args:
        postX (list): X coordinates to transform (1D).
        postY (list): Y coordinates to transform (1D).
        transX (float): Translation value in X dimension.
        transY (float): Translation value in Y dimension.
        rot (float): Rotation value in degree.

    Returns:
        transformed (list): Transformed coordinates (2D).

    """
    tform = tf.SimilarityTransform(scale=1, rotation=rot, translation=(transX, transY))
    transformed = tf.matrix_transform(np.transpose([postX, postY]), tform.params)
    return transformed

def fiducialsAlignment(MFA):
    """Define the coordinate transform parameters leading to the optimal overlap between the pre and post-MALDI
    fiducials.

     Args:
         MFA (str): path to Main Folder Analysis.

     """

    def get_distance(x_spots,y_spots,xe,ye,n_neighbor):
        """Measure the euclidean distance between each point of an array to its n nearest neighbor in a second array using
        kd-tree algorithm.

        Args:
            x_spots (list): X coordinates of the array to query (1D).
            y_spots (list): Y coordinates of the array to query (1D).
            xe (list): X coordinates of the array to index (1D).
            ye (list): Y coordinates of the array to index(1D).
            n_neighbor (int): The number of nearest neighbor to consider.

        Returns:
            distances (list): Distances of the indexed points n nearest queried neighbors (2D).

        """
        data = list(zip(xe.ravel(), ye.ravel()))
        tree = spatial.KDTree(data)
        return tree.query(list(zip(x_spots.ravel(), y_spots.ravel())),n_neighbor)

    def err_func(params, preX, preY, postX, postY):
        """Error function passed in the optimizer. Transforms coordinates of target frame and returns the mean nearest neighbor
        distance to the 1st frame fiducials.

        Args:
            params (list): Array of coordinate transformation function:
                           [Translation in X(float), Translation in Y(float), rotation(float)] (1D).
            preX (list): X coordinates of the 1st frame fiducials (1D).
            preY (list):  Y coordinates of the 1st frame fiducials (1D).
            postX (list): X coordinates of the target frame fiducials (1D).
            postY (list): Y coordinates of the target frame fiducials (1D).

        Returns:
            mean_distances (): Mean N nearest neighbor distance to the 1st frame fiducials.

        """
        transX, transY, rot = params
        transformed = transform(postX, postY, transX, transY, rot)
        distances = np.array(get_distance(transformed[:, 0], transformed[:, 1], preX, preY, 1)[0])
        return np.mean(distances)

    preX, preY = np.load(MFA + 'Fiducials/preXYpenmarks.npy')
    postX, postY = np.load(MFA + 'Fiducials/postXYpenmarks.npy')

    n_features = 2000
    post_den = int(np.round(np.shape(postX)[0] / n_features))
    postX_redu = postX[::post_den]#reduces optimizer computation time
    #TODO --> need to evaluate impact on alignment precision
    postY_redu = postY[::post_den]
    # print('post features length = {}'.format(np.shape(postX_redu)))

    pre_den = int(np.round(np.shape(preX)[0] / n_features))
    preX_redu = preX[::pre_den]
    preY_redu = preY[::pre_den]
    # print('pre features length = {}'.format(np.shape(preX_redu)))

    minF = basinhopping(err_func, x0=(0, 0, 0 ),
                        niter=1, T=1.0, stepsize=10,
                        minimizer_kwargs={'args': ((preX_redu, preY_redu, postX_redu, postY_redu))},
                        take_step=None, accept_test=None, callback=None, interval=50, disp=True,
                        niter_success=1)
    # print(minF)

    # #Case of 180deg rotation:
    # deg = 180
    # rad =  -deg / 180. * np.pi
    # transformed = transform(postX_redu, postY_redu, 0, 0, rad)
    #
    # minF_rot = basinhopping(err_func, x0=((np.mean(postX) - np.mean(transformed[:, 0])), (np.mean(postY) - np.mean(transformed[:, 1])), rad),
    #                     niter=1, T=1.0, stepsize=10,
    #                     minimizer_kwargs={'args': ((preX_redu, preY_redu, postX_redu, postY_redu))},
    #                     take_step=None, accept_test=None, callback=None, interval=50, disp=True,
    #                     niter_success=1)
    # print(minF_rot)
    # if minF.fun < minF_rot.fun:

    transX = minF.x[0]
    transY = minF.x[1]
    rot = minF.x[2]
    # scale = minF.x[3]
    #     print('rot = 0 deg')
    # else:
    #     transX = minF_rot.x[0]
    #     transY = minF_rot.x[1]
    #     rot = minF_rot.x[2]
    #     # scale = minF.x[3]
    #     print('rot = 180 deg')

    np.save(MFA + '/Fiducials/optimized_params.npy', [transX, transY, rot])
    transformed = transform(postX, postY, transX, transY, rot)
    plt.figure()
    plt.scatter(transformed[:, 0], transformed[:, 1],1)
    plt.scatter(preX, preY, 1, 'r')
    plt.axis('equal')
    plt.savefig(MFA + '/Fiducials/surfRegResults.png', dpi = 500)
    plt.close('all')

def TransformMarks(MFA):
    """Transform the ablation mark coordinates from the post-MALDI dataset using the geometric transform parameters
    defined in SURF_Alignment() function to estimate their position in the pre-MALDI dataset.

     Args:
         MFA (str): path to Main Folder Analysis.

     """
    xe_clean2, ye_clean2 = np.load(MFA + '/gridFit/xye_clean2.npy')
    x_spots, y_spots = np.load(MFA + '/gridFit/xye_grid.npy')
    transX, transY, rot = np.load(MFA + '/Fiducials/optimized_params.npy')
    scale = 1
    shape, pix_size = np.load(MFA + '/gridFit/metadata.npy')
    xye_tf = transform(xe_clean2, ye_clean2, transX, transY, rot)
    X = xye_tf[:, 0]
    Y = xye_tf[:, 1]
    np.save(MFA + '/Fiducials/transformedMarks.npy', [X, Y])
    xyg_tf = transform(x_spots.ravel(), y_spots.ravel(), transX, transY, rot)
    Xg = xyg_tf[:, 0]
    Yg = xyg_tf[:, 1]
    np.save(MFA + '/Fiducials/transformedGrid.npy', [Xg, Yg])

    if os.path.exists(MFA + 'gridFit/marksMask.npy'):
        marksMask = np.load(MFA + 'gridFit/marksMask.npy')
        tfMarksMask = []
        for i in range(np.shape(marksMask)[0]):
            if np.shape(marksMask[i][0].T)[0] > 1:
                tfMask = transform(marksMask[i][0].T, marksMask[i][1].T, transX, transY, rot)
                tfMarksMask.append([tfMask[:, 0], tfMask[:, 1]])
            else:
                tfMarksMask.append([[], []])
                print('empty')
        np.save(MFA + '/Fiducials/transformedMarksMask.npy', tfMarksMask)

    tfMarksMask = np.array(tfMarksMask)
    penmarks = np.load(MFA + 'Fiducials/preXYpenmarks.npy')
    plt.figure(figsize=[50,25])
    plt.scatter(penmarks[0,:], penmarks[1,:], 1, c='k')
    plt.axis('equal')
    for i in range(np.shape(tfMarksMask)[0]):
        if np.shape(tfMarksMask[i][0])[0] > 1:
            plt.scatter(tfMarksMask[i][0], tfMarksMask[i][1], 1, 'r')
    plt.scatter(X, Y, 1, c='g')
    plt.savefig(MFA + '/Fiducials/registration_result.png', dpi=100)
    plt.close('all')
