import tifffile as tiff
import numpy as np
import glob, os, spaceM
from matplotlib import pyplot as plt
from scipy import spatial
import codecs
from scipy.optimize import basinhopping
import re
import itertools
import scipy.io
import spaceM.ImageFileManipulation.manipulations as manip
import tqdm
import seaborn as sns
from scipy import ndimage
from skimage.measure import label, regionprops
import spaceM.ImageFileManipulation.FIJIcalls as fc
from PIL import Image, ImageFile
import tifffile as tiff
ImageFile.LOAD_TRUNCATED_IMAGES = True
Image.MAX_IMAGE_PIXELS = None
from scipy.optimize import curve_fit
from scipy import exp


def crop_img(MF, im_p, #=MF+'Analysis/StitchedMicroscopy/postMALDI_FLR/img_t1_z1_c1',
             coords_p):#=MF+'Analysis/gridFit/AM_cropped_cropCoords.npy'):
    im = plt.imread(im_p)
    coords = np.load(coords_p, allow_pickle=True).item()
    im_crop = im[int(np.round(coords['topLeft'][1])): int(np.round(coords['bottomRight'][1])),
              int(np.round(coords['topLeft'][0])): int(np.round(coords['bottomRight'][0]))]
    # plt.imshow(im[10465:15560, 5011:10357])
    tiff.imsave(file=MF+'Analysis/gridFit/AM_cropped_2.tif', data=im_crop)


def DHB_prep(MF, path1, path2):
    im1 = contrast(plt.imread(path1), 0.01, 0.999)
    im2 = contrast(plt.imread(path2), 0.01, 0.999)
    im3 = scale(contrast(scale(scale(im2) * scale(im1) * -1), 0.9, 1))* 65535
    tiff.imsave(MF + 'Analysis/gridFit/AM_cropped_3.png', im3.astype(np.uint16))
    return im3


def scale(arr):
    """Scale array between 0 and 1"""
    return (arr - np.min(arr)) / (np.max(arr) - np.min(arr))


def contrast(arr, min, max=1.0):
    """Clip array between min and max values"""
    return np.clip(arr, np.percentile(arr, min*100), np.percentile(arr, max*100))


def spotFinder(MF, path1, path2=None, matrix='DHB', show_results=True):
    """Detect ablation marks on the tiled images. For details, see the paper.
    Args:
        path1 (str): path of the image to detect ablation marks on.

    """
    if matrix == 'DHB':
        img_i = DHB_prep(MF, path1, path2)
        img_2 = scale(contrast(scale(ndimage.gaussian_filter(img_i, sigma=6)), 0, 0.99))
        # img_segmentation = scale(contrast(scale(ndimage.gaussian_filter(img_i, sigma=3)), 0, 0.99)) * 65535
    elif matrix == 'DAN':
        img_i = plt.imread(path1)
        img = scale(img_i)
        contrast_min = np.mean(img) + 1.5 * np.std(img)  # 2 for DAN, 1.5 for DHB
        if contrast_min >= 1: contrast_min = 0.8
        img_2 = contrast(img, min=contrast_min)

    ff = np.fft.fft2(img_2)
    F1 = 20*np.log10(np.abs(np.fft.fftshift(ff)))

    mask1 = F1 - ndimage.gaussian_filter(F1, sigma=15)
    mask1[mask1<0] = 0

    int_thresh = np.mean(mask1) + 3.5*np.std(mask1) #3.5 for DAN
    struct = ndimage.generate_binary_structure(2,1)
    mask2 = ndimage.binary_dilation(mask1>int_thresh, structure=struct, iterations=5).astype(mask1.dtype) #20 for DAN
    mask2_blur = scale(ndimage.gaussian_filter(mask2, sigma=2)) #lower sigma means more grid but more leaking
    ff_masked = np.fft.fftshift(mask2_blur)*ff #mask_2 for DAN, mask2_blur for DHB

    plt.imshow(mask2_blur)

    freq_up = 0.6 #0.6 for DAN
    freq_down = 0.0 #0 for DAN

    # freq_up = 0.2 #0.6 for DAN
    # freq_down = 0.0 #0 for DAN

    [N, M] = np.shape(img_2)
    dx = 1
    KX0 = (np.mod(1 / 2 + np.arange(0,M) / M, 1) - 1 / 2)
    KX1 = KX0 * (2 * np.pi / dx)
    KY0 = (np.mod(1 / 2 + np.arange(0,N) / N, 1) - 1 / 2)
    KY1 = KY0 * (2 * np.pi / dx)
    [KX, KY] = np.meshgrid(KX1, KY1)
    lpf = (KX * KX + KY * KY < freq_up ** 2)
    lpf2 = (KX * KX + KY * KY < freq_down ** 2)
    mix = lpf * ~lpf2
    mix_masked = mix * ff_masked
    F2 = 20 * np.log10(abs(np.fft.fftshift(mix_masked)) + 1)
    rec = np.real(np.fft.ifft2(mix_masked))
    rec = scale(rec) * 65535
    rec = contrast(rec, 0.1, 0.999)
    plt.imshow(rec)

    tiff.imsave(MF + 'Analysis/gridFit/FT_filtered.tiff', rec.astype(np.uint16))

    # plt.subplot(121), plt.imshow(F2)
    # plt.subplot(122), plt.imshow(rec)

    rec_bw = rec > np.mean(rec) + 1.5*np.std(rec) #2.5 for DAN
    # ax = plt.subplot(121)
    # plt.imshow(rec_bw)
    # plt.subplot(122, sharex=ax, sharey=ax), plt.imshow(img_2)

    label_img = label(rec_bw, connectivity=rec_bw.ndim)
    props = regionprops(label_img)
    centroids = np.array([[props[i].centroid[1],props[i].centroid[0]]  for i in range(len(props))])

    if show_results:
        def no_axis(ax):
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.get_xaxis().tick_bottom()
            ax.get_yaxis().tick_left()
            ax.set_xticks([])
            ax.set_yticks([])
            return ax

        plt.figure(figsize=(18,9))
        ax = no_axis(plt.subplot(241))
        plt.imshow(img_i, cmap='gray')
        no_axis(plt.subplot(245, sharex=ax, sharey=ax)), plt.imshow(img_2, cmap='gray')

        ax1 = no_axis(plt.subplot(242))
        plt.imshow(F1, cmap='gray')
        no_axis(plt.subplot(246, sharex=ax1, sharey=ax1)), plt.imshow(F2, cmap='gray')

        no_axis(plt.subplot(243, sharex=ax, sharey=ax)), plt.imshow(rec, cmap='gray')
        no_axis(plt.subplot(247, sharex=ax, sharey=ax)), plt.imshow(rec_bw, cmap='gray')
        x = [plt.scatter(x,y,100, color=[1,0,0]) for x,y in centroids]

        no_axis(plt.subplot(248, sharex=ax, sharey=ax)), plt.imshow(img_i, cmap='gray')
        x = [plt.scatter(x,y,100, color=[1,0,0]) for x,y in centroids]
        plt.show()
        plt.tight_layout()
        ax.set_ylim([1000,1200])
        ax.set_xlim([1000,1200])
        ax1.set_ylim([1500,3000])
        ax1.set_xlim([1500,3000])
        # mng = plt.get_current_fig_manager()
        # mng.window.showMaximized()
        plt.savefig(MF + 'Analysis/gridFit/spotFinder_diagnostic.png')
        plt.close('all')
    return centroids


def MarkFinderFT(MF):
    """Find center of mass of ablation marks on individual tile images from postMALDI microscopy dataset.

    Args:
        MF (str): path to Main Folder.

    """
    MFA = MF + 'Analysis/'
    MFA_Spom = MFA + 'StitchedMicroscopy/postMALDI_FLR/'

    [POST_picXcoord, POST_picYcoord] = fc.readTileConfReg(MFA_Spom)
    allxScaled = []
    allyScaled = []
    overlap_img = 0.10
    picNdetect = []
    picMean = []
    img = tiff.imread(MFA_Spom + os.listdir(MFA_Spom)[10])
    if np.shape(img.shape)[0] > 2:
        img = img[0, :, :]  # TODO: uncomment
    for item in os.listdir(MFA_Spom):
        if item.endswith('.tif'):
            ##NEW APPROACH
            # centroids = np.asarray(eng.spotFinder2(MFA_Spom + item, 3)) #TODO: remove matlab dependency
            centroids = spotFinder(MFA_Spom + item)
            picInd = int(item[len(item)-7:len(item)-4]) # 6-->7
            picNdetect = np.append(picNdetect, len(centroids))
            picMean = np.append(picMean, np.mean(img))
            for i in centroids:
                if i[0] < img.shape[0] - img.shape[0]*overlap_img and i[1] < img.shape[1] - img.shape[1]*overlap_img:
                    xScaled = i[1] + POST_picXcoord[picInd - 1] #had to invert X and Y coz of Matlab output
                    yScaled = i[0] + POST_picYcoord[picInd - 1]
                    allxScaled = np.append(allxScaled, xScaled)
                    allyScaled = np.append(allyScaled, yScaled)
            print('Detecting marks in image # {}; Total = {}'\
                .format(picInd, len(allxScaled)))
    np.save(MFA + 'gridFit/ablation_marks_XY.npy', \
            [allxScaled, allyScaled])
    np.save(MFA + 'gridFit/pics_nDetect_Mean.npy', \
            [picNdetect, picMean])
    # return allxScaled, allyScaled
    # eng.quit()

def getPixSize(MFI):
        """Reads the pixel size in um from the Nikon Ti E microscope (NIS elements software).

        Args:
            MFI (str): path to Main Folder Input.

        Returns:
            pix_size (float): pixel size in um.

        """
        txt_file = codecs.open(MFI + '/Microscopy/postMALDI/out.txt', 'r', 'utf-16')
        for row in txt_file:
            if row.startswith('Calibration'):
                # print(row)
                pix_size = float(row.strip().split()[2].replace(',', '.'))
            # else:
            #     pix_size = 0.73
        return pix_size

def GridFit(MF, optimization=False, manual_cleaning=True, MarkFinderFT=True, matrix='DAN'):
    """Fit a theoretical grid on the ablation marks coordinates to remove extra detections and re-index ablation marks
    uniformly.

    Args:
        MFA (str): path to Main Folder Analysis.
        MFI (str): path to Main Folder Input.
        optimization (bool): whether or not an optimization should be used to fit the grid on the ablation marks.
        manual_cleaning (bool): whether or not the ablation amrks coordinates have been cleaned with the
            script 'manualCleaning.py'.
        MarkFinderFT (bool): whether or not a second ablation mark detection is performed on the stitched microsocpy
            cropped around ablation marks. This should be done by default as it should improve the detections.

    """

    def create_grid(shape, rotation_deg, affine_lat, i, j):
        """Create node coordinates of a theoretical grid.

        Args:
            shape (array): number of rows and columns (1D).
            rotation_deg (float): rotation value un degree of the grid.
            affine_lat (float): spacing between the nodes.
            i (float): X coordinate of the center of the grid.
            j (float): Y coordinate of teh center of the grid.

        Returns:
            x_spots (array): X coordinates of the theoretical grid nodes (1D).
            y_spots (array): Y coordinates of the theoretical grid nodes (1D).

        """
        x_ther, y_ther = np.meshgrid((np.arange(shape[1]) - (shape[1] - 1) / 2.) * affine_lat, \
                                     (np.arange(shape[0]) - (shape[0] - 1) / 2.) * affine_lat)
        theta_t = -rotation_deg / 180. * np.pi
        x_spots = i + np.cos(theta_t) * (x_ther - i) - np.sin(theta_t) * (y_ther - j)
        y_spots = j + np.sin(theta_t) * (x_ther - i) + np.cos(theta_t) * (y_ther - j)
        center_xs = np.mean(x_spots)
        center_ys = np.mean(y_spots)
        x_spots = x_spots + (i - center_xs)
        y_spots = y_spots + (j - center_ys)
        return x_spots, y_spots

    def get_distance(x_spots, y_spots, xe, ye, n_neighbor):
        """Measure the euclidean distance between each point of an array to its n nearest neighbor in a second array using
        kd-tree algorithm.

        Args:
            x_spots (array): X coordinates of the array to query (1D).
            y_spots (array): Y coordinates of the array to query (1D).
            xe (array): X coordinates of the array to index (1D).
            ye (array): Y coordinates of the array to index(1D).
            n_neighbor (int): The number of nearest neighbor to consider.

        Returns:
            distances (array): Distances of the indexed points n nearest queried neighbors (2D).

        """
        data = list(zip(xe.ravel(), ye.ravel()))
        tree = spatial.KDTree(data)
        distances = tree.query(list(zip(x_spots.ravel(), y_spots.ravel())), n_neighbor)
        return distances

    def err_func(params,shape, xe, ye):
        """Error function passed in the optimizer. It defines the optimal grid parameters leading to the lowest
        mean distance between the grid nodes and their nearest ablation mark.

        Args:
            params (array): array of parameters to build the theoretical grid which are optimized (1D).
            shape (array): 'shape' parameter of the grid, which is non optimized (1D).
            xe (array): X coordinates of the ablation marks (1D).
            ye (array): Y coordinates of the ablation marks (1D).

        Returns:
            distance (float): mean distance between the grid nodes and theire nearest ablation marks summed to the
            squared distance between the extrema of the theoretical grid and the ablation mark coordinates.

        """
        rotation_deg, affine_lat, i, j = paramsw
        x_spots, y_spots = create_grid(shape,rotation_deg,affine_lat,i, j)
        distances, coord = get_distance(x_spots,y_spots,xe,ye,1)
        top1, bottom1, right1, left1 = getExtrema(shape, x_spots, y_spots)
        top2, bottom2, right2, left2 = getExtrema(shape, xe, ye)
        border_dist = abs(top1-top2) + abs(bottom1-bottom2) + abs(left1-left2) + abs(right1-right2)
        distance = np.mean(distances) + border_dist**2
        return distance

    def getShapeRes(MFI):
        """Read the number of rows, columns and the spacing (in um) between ablation marks of the MALDI acquisition in
        the UDP file.

        Args:
            MFI (str): path to Main Folder Input.

        Returns:
            shape (array): number of rows and columns of the MALDI acquisition (1D).
            resX (float): spacing in um between ablation marks in X dimension.
            resY (float): spacing in um between ablation marks in Y dimension.

        """
        UDP_path = glob.glob(MFI + '/MALDI/*.UDP')[0]
        UDP_file = codecs.open(UDP_path, 'r')
        shape = []
        for row in UDP_file:
            if row.startswith('    <MaxX>'):
                # print row.strip().replace(',','.')
                m = re.search('(?<=>)\d+', row.strip())
                shape.append(float(m.group(0)))
            elif row.startswith('    <MaxY>'):
                 # print row.strip().replace(',','.')
                m = re.search('(?<=>)\d+', row.strip())
                shape.append(float(m.group(0)))
            elif row.startswith('    <ResolutionX>'):
                # print row.strip().replace(',','.')
                m = re.search('(?<=>)\d+\.\d*', row.strip().replace(',', '.'))
                resX = float(m.group(0))
            elif row.startswith('    <ResolutionY>'):
                # print row.strip().replace(',','.')
                m = re.search('(?<=>)\d+\.\d*', row.strip().replace(',', '.'))
                resY = float(m.group(0))
        return shape, resX, resY

    def estimateAngle(xe,  ye, shape, MFA, Figures=True):
        """Estimate the relative angle between the ablation marks and the X axis.
        TODO: estimate in both X and Y and get average estimate. It might be more accurate
        Args:
            xe (array): X coordinates of the ablation marks (1D).
            ye (array): Y coordinates of the ablation marks (1D).
            shape (array): number of rows and columns of the MALDI acquisition (1D).
            MFA (str): path to Main Folder Analysis.
            Figures (bool): whether or not plot the results and save in analysis folder.

        Returns:
            rotation_deg (float): alignment rotation angle in degree to the X axis.

        """
        counts_x = []
        counts_y = []
        angle = []
        for i in np.linspace(-15, 15, 5000):
            rotation = i
            theta = rotation / 180. * np.pi
            x_spots = np.mean(xe) + np.cos(theta) * (xe - np.mean(xe)) \
                      - np.sin(theta) * (ye - np.mean(ye))
            y_spots = np.mean(ye) + np.sin(theta) * (xe - np.mean(xe)) \
                      + np.cos(theta) * (ye - np.mean(ye))
            a_x = np.histogram(x_spots, int(shape[0]) * 100)
            a_y = np.histogram(y_spots, int(shape[0]) * 100)
            count_x = np.asarray(a_x[0][a_x[0] > 0]).shape[0]
            count_y = np.asarray(a_y[0][a_y[0] > 0]).shape[0]
            # print count, theta
            angle.append(rotation)
            counts_x.append(count_x)
            counts_y.append(count_y)
            # print '{}/{} deg'.format(i, 15)
        rotation_deg_x = angle[np.where(counts_x == np.min(counts_x))[0][0]]
        # rotation_deg_y = angle[np.where(counts_y == np.min(counts_y))[0][0]]
        # rotation_deg = np.mean([rotation_deg_x,rotation_deg_y])
        print('Rough Rotation estimation is {} degree'.format(rotation_deg_x))
        if Figures == True:
            plt.figure()
            plt.plot(angle, counts_x, 'ko')
            plt.plot(rotation_deg_x, np.min(counts_x), 'ro', markersize=20, label='Global minimum')
            plt.legend()
            plt.xlabel('Rotation angles (rad)', fontsize=20)
            plt.ylabel('Number of non-zero bins \nof data 1D projection', fontsize=20)
            plt.title('Angle Estimation: angle=%.3f degree' % (rotation_deg_x), fontsize=25)
            plt.savefig(MFA + '/gridFit/estimateAngle.png')
            plt.close('all')
        return rotation_deg_x

    def cleanData(rotation, resX, resY, xe, ye, tolerance, manual_cleaning=True, Figures=True):
        """Remove outlier datapoints by estimating the ablation mark extrema in X and Y dimension.

        Args:
            rotation (float): alignment rotation angle in degree to the X axis.
            resX (float): spacing in um between ablation marks in X dimension.
            resY (float): spacing in um between ablation marks in Y dimension.
            xe (array): X coordinates of the ablation marks (1D).
            ye (array): Y coordinates of the ablation marks (1D).
            tolerance (int): cutoff tolerance.
            manual_cleaning (bool): whether or not the ablation amrks coordinates have been cleaned with the
                script 'manualCleaning.py'.
            Figures (bool): whether or not plot the results and save in analysis folder.

        Returns:
            xe (array): X coordinates of the ablation marks (1D).
            ye (array): Y coordinates of the ablation marks (1D).
            xe_r (array): X coordinates of the rotated ablation marks (1D).
            ye_r (array): Y coordinates of the rotated ablation marks (1D).

        """
        theta = rotation / 180. * np.pi
        xe_r = np.mean(xe) + np.cos(theta) * (xe - np.mean(xe)) \
               - np.sin(theta) * (ye - np.mean(ye))
        ye_r = np.mean(ye) + np.sin(theta) * (xe - np.mean(xe)) \
               + np.cos(theta) * (ye - np.mean(ye))
        if manual_cleaning:
            xr_clean=xe_r
            yr_clean=ye_r
            xe_clean=xe
            ye_clean=ye
        else:
            # remove outlier points to facilitate center estimation
            nbins_x = int((np.max(xe_r)-np.min(xe_r)) / resX)
            nbins_y = int((np.max(ye_r)-np.min(ye_r)) / resY)
            c_x, b_x = np.histogram(xe_r, nbins_x)
            c_y, b_y = np.histogram(ye_r, nbins_y)
            # treshc_x = tresh[0]
            # treshc_y = tresh[1]
            treshc_x = np.mean(c_x[c_x>0.0])
            treshc_y = np.mean(c_y[c_y>0.0])
            tx_low = np.min(b_x[:-1][c_x > treshc_x]) - tolerance
            tx_high = np.max(b_x[:-1][c_x > treshc_x]) + tolerance
            ty_low = np.min(b_y[:-1][c_y > treshc_y]) - tolerance
            ty_high = np.max(b_y[:-1][c_y > treshc_y]) + tolerance
            ind2d = np.ravel([np.where(xe_r < tx_low)[0], np.where(xe_r > tx_high)[0], np.where(ye_r < ty_low)[0],
                              np.where(ye_r > ty_high)[0]])
            ind_cut = np.unique(list(itertools.chain.from_iterable(ind2d)))

            xr_clean = np.delete(xe_r, ind_cut)
            yr_clean = np.delete(ye_r, ind_cut)
            xe_clean = np.delete(xe, ind_cut)
            ye_clean = np.delete(ye, ind_cut)

            plt.figure()
            x_line, = plt.plot(b_x[:-1], c_x, 'r', label = 'X_projections')
            thresh_x_line, = plt.plot(b_x[:-1], np.ones(np.shape(b_x[:-1])) * treshc_x, 'r', label='Threshold X')
            y_line, = plt.plot(b_y[:-1], c_y, 'k', label = 'Y_projections')
            thresh_y_line, = plt.plot(b_y[:-1], np.ones(np.shape(b_y[:-1])) * treshc_y, 'k', label='Threshold Y')
            plt.legend(handles = [x_line, y_line, thresh_x_line, thresh_y_line])#, ['X projections', 'Y projections'])
            plt.savefig(MFA + '/gridFit/XYprojections.png')
            # pics_ind_clean = np.delete(pics_ind, ind_cut)
        if Figures:

            crop = np.load(MFA + 'gridFit/AM_cropped_coords.npy', allow_pickle=True).item()
            minX = crop['topLeft'][0]
            minY = crop['topLeft'][1]
            xe_scaled = xe/pix_size - minY
            ye_scaled = ye/pix_size - minX
            xe_clean_scaled = xe_clean/pix_size - minY
            ye_clean_scaled = ye_clean/pix_size - minX
            plt.figure(figsize=(60,60))
            plt.imshow(plt.imread(MFA + 'gridFit/AM_cropped.tif'), cmap='gray')
            plt.scatter(ye_scaled, xe_scaled, 5, label='Input')
            plt.scatter(ye_clean_scaled, xe_clean_scaled, 2, color=[0,1,0], label='Output - cleaned')
            plt.savefig(MFA + '/gridFit/Cleaning_results.png', dpi=200)
            plt.tight_layout()
            plt.close('all')

        return xe_r, ye_r, xr_clean, yr_clean, xe_clean, ye_clean
        #, pics_ind_clean

    def getExtrema(shape, xe_r, ye_r):
        """Get extrema coordinates of the ablation mark grid.

        Args:
            shape (array): number of rows and columns of the MALDI acquisition (1D).
            xe_r (array): X coordinates of the rotated ablation marks (1D).
            ye_r (array): Y coordinates of the rotated ablation marks (1D).

        Returns:
            top (float): highest Y coordinate of the ablation marks.
            bottom (float): lowest Y coordinate of the ablation marks.
            right (float): highest X coordinate of the ablation marks.
            left (float): lowest X coordinate of the ablation marks.

        """
        #Get extreme coordinates of the acquisition matrix
        xe_r = xe_r.ravel()
        ye_r = ye_r.ravel()
        X1 = []
        Y1 = []
        X2 = []
        Y2 = []
        X3 = []
        Y3 = []
        X4 = []
        Y4 = []
        for i in range(int(shape[0])):
            ii = i + 1
            # horizontal up
            X1 = np.append(X1, xe_r[np.where(ye_r == np.sort(ye_r)[np.shape(ye_r)[0] - ii])])
            Y1 = np.append(Y1, ye_r[np.where(ye_r == np.sort(ye_r)[np.shape(ye_r)[0] - ii])])
            # horizontal low
            X2 = np.append(X2, xe_r[np.where(ye_r == np.sort(ye_r)[i])])
            Y2 = np.append(Y2, ye_r[np.where(ye_r == np.sort(ye_r)[i])])
            # vertical right
            X3 = np.append(X3, xe_r[np.where(xe_r == np.sort(xe_r)[np.shape(xe_r)[0] - ii])])
            Y3 = np.append(Y3, ye_r[np.where(xe_r == np.sort(xe_r)[np.shape(xe_r)[0] - ii])])
            # vertical left
            X4 = np.append(X4, xe_r[np.where(xe_r == np.sort(xe_r)[i])])
            Y4 = np.append(Y4, ye_r[np.where(xe_r == np.sort(xe_r)[i])])
        top = np.median(Y1)
        bottom = np.median(Y2)
        right = np.median(X3)
        left = np.median(X4)
        return top, bottom, right, left

    def estimateCenter(shape, xr_clean, yr_clean, MFA, Figures):
        """Estimate the X and Y coordinates of the center of the ablation mark grid.

        Args:
            shape (array): number of rows and columns of the MALDI acquisition (1D).
            xr_clean (array): X coordinates of the rotated ablation marks (1D).
            yr_clean (array): Y coordinates of the rotated ablation marks (1D).
            MFA (str): path to Main Folder Analysis.
            Figures (bool): whether or not plot the results and save in analysis folder.

        Returns:
            center_x (float): X coordinates of the center of the ablation mark grid.
            center_y (float): Y coordinates of the center of the ablation mark grid.
            top (float): highest Y coordinate of the ablation marks.
            bottom (float): lowest Y coordinate of the ablation marks.
            right (float): highest X coordinate of the ablation marks.
            left (float): lowest X coordinate of the ablation marks.

        """
        top, bottom, right, left = getExtrema(shape, xr_clean, yr_clean)
        center_x = left + abs(right - left) / 2
        center_y = bottom + abs(top - bottom) / 2
        if Figures == True:
            plt.figure()
            plt.plot([left, right, right, left, left], [bottom, bottom, top, top, bottom], 'r', linewidth=3.0 )
            # plt.scatter(xr, yr, 10, 'k')
            plt.scatter(xr_clean, yr_clean, 10)
            plt.scatter(center_x, center_y, 300, 'r', label='Estimated center')
            plt.xlabel('X dimension (um)')
            plt.ylabel('Y dimension (um)')
            plt.title('Rough center estimation\n')
            # plt.legend()
            plt.axis('equal')
            # print 'Rough center estimation: X = {}; Y = {}'.format(center_x, center_y)
            plt.savefig(MFA + '/gridFit/estimateCenter.png')
            plt.close('all')
        return center_x, center_y, top, bottom, right, left

    def cutData(xe_clean, ye_clean, xr_clean, yr_clean, bottom, top, left, right):
        """Removes outlier points (if outside of estimated acquisition extreme coordinates).
        This function is no longer used.

        Args:
            xe_clean (array): X coordinates of the ablation marks (1D).
            ye_clean (array): Y coordinates of the ablation marks (1D).
            xr_clean (array): X coordinates of the rotated ablation marks (1D).
            yr_clean (array): Y coordinates of the rotated ablation marks (1D).
            top (float): highest Y coordinate of the ablation marks.
            bottom (float): lowest Y coordinate of the ablation marks.
            right (float): highest X coordinate of the ablation marks.
            left (float): lowest X coordinate of the ablation marks.

        Returns:
            xe_cut (array): X coordinates of the ablation marks without outliers(1D).
            ye_cut (array): X coordinates of the ablation marks without outliers(1D).

        """
        inds = np.array(yr_clean < bottom) + np.array(yr_clean>top) + np.array(xr_clean>right) + np.array(xr_clean<left)
        xe_cut = np.delete(xe_clean, inds)
        ye_cut = np.delete(ye_clean, inds)

        plt.scatter(xe_clean[~inds], ye_clean[~inds], 10, 'b')
        plt.scatter(xe_clean[inds], ye_clean[inds], 10, 'r')

        return xe_cut, ye_cut

    def estimateLattice(shape, rotation_esti, center_x_esti, center_y_esti, top, bottom, right, left, xe, ye, MFA, Figures):
        """Estimate the spacing between ablation marks. Different theoretical grid is generated with increasing spacing
        between the nodes and the mean distance between each node and their closes ablation marks is reported.
        The function returns the spacing which led to the smallest distance.

        Args:
            shape (array): number of rows and columns of the MALDI acquisition (1D).
            rotation_esti (float): alignment rotation angle in degree to the X axis.
            center_x_esti (float): X coordinates of the center of the ablation mark grid.
            center_y_esti (float): Y coordinates of the center of the ablation mark grid.
            top (float): highest Y coordinate of the ablation marks.
            bottom (float): lowest Y coordinate of the ablation marks.
            right (float): highest X coordinate of the ablation marks.
            left (float): lowest X coordinate of the ablation marks.
            xe (array): X coordinates of the ablation marks (1D).
            ye (array): Y coordinates of the ablation marks (1D).
            MFA (str): path to Main Folder Analysis.
            Figures (bool): whether or not plot the results and save in analysis folder.

        Returns:
            affine_lat2 (float): affined estimation of the spacing between ablation marks.

        """
        mu = np.mean(
            [np.abs((np.median(bottom) - np.median(top)) / shape[1]), np.abs((np.median(left) - np.median(right)) / shape[0])])
        dist = []
        steps = []
        nsteps = 30
        for i in np.linspace(mu - 3, mu + 3, nsteps):
            # step = i
            x_spots, y_spots = create_grid(shape, rotation_esti, i, center_x_esti, center_y_esti)
            a,b = get_distance(x_spots, y_spots, xe, ye, 1)
            # print np.mean(a)
            dist.append(np.mean(a))
            steps.append(i)
            # print '{}/{}'.format(i, mu + 0.7)
        affine_lat1 = steps[np.where(dist == np.min(dist))[0][0]]
        if Figures == True:
            plt.figure()
            plt.plot(steps, dist)
            plt.plot(steps, dist, 'o', color=[0,1,0], markersize=5)
            plt.xlabel('Lattice size (um)', fontsize=20)
            plt.ylabel('Mean closest neighbor from\nTheoretical to Experimental grid', fontsize=20)
            plt.title('Rough lattice estimation', fontsize=25)
            plt.scatter(affine_lat1, np.min(dist), 200, 'r', label='Global minimum')
            plt.legend()
            plt.savefig(MFA + '/gridFit/estimateLattice1.png')
            plt.close('all')
        # return affine_lat

        if affine_lat1 < resX - 1 or affine_lat1 > resX + 1:
            affine_lat1 = resX

        dist = []
        steps = []
        nsteps = 100
        for i in np.linspace(affine_lat1 - 1, affine_lat1 + 1, nsteps):
            # step = i
            x_spots, y_spots = create_grid(shape, rotation_esti, i, center_x_esti, center_y_esti)
            a,b = get_distance(x_spots, y_spots, xe, ye, 1)
            # print np.mean(a)
            dist.append(np.mean(a))
            steps.append(i)
            # print '{}/{}'.format(i, mu + 0.7)
        affine_lat2 = steps[np.where(dist == np.min(dist))[0][0]]
        if Figures == True:
            plt.figure()
            plt.plot(steps, dist)
            plt.plot(steps, dist, 'o', color=[0,1,0], markersize=5)
            plt.xlabel('Lattice size (um)', fontsize=20)
            plt.ylabel('Mean closest neighbor from\nTheoretical to Experimental grid', fontsize=20)
            plt.title('Fine lattice estimation', fontsize=25)
            plt.scatter(affine_lat2, np.min(dist), 200, 'r', label='Global minimum')
            plt.legend()
            plt.savefig(MFA + '/gridFit/estimateLattice2.png')
            plt.close('all')
        return affine_lat2

    MFI = MF + 'Input/'
    MFA = MF + 'Analysis/'

    if MarkFinderFT:
        # Ablation mark detection efficiency can be affected by the number of ablation marks present in the image.
        # To improve the results, the same algorithm is run a seond time on the cropped stitched microsocpy around the
        # previous detections. That way all ablation marks are present in the same image, maximizing their associated
        # frequencies in the fourrier domain.

        # window = 200
        # if not os.path.exists(MFA + 'gridFit/SpotFinder_cropped_window200.png'):
        #     manip.crop2coords(MFA + 'gridFit/ablation_marks_XY.npy',
        #                       MFA + 'StitchedMicroscopy/postMALDI_FLR/img_t1_z1_c1',
        #                       MFA + 'gridFit/2nd_detection_cropped_window200.png', window=window)

        # if not os.path.exists(MFA + 'gridFit/ablation_marks_XY_2nd_detection.npy'):
        # X, Y = np.load(MFA + 'gridFit/ablation_marks_XY.npy')
        # minX = np.min(X) - window
        # minY = np.min(Y) - window
        # eng = matlab.engine.start_matlab()
        # centroids = np.asarray(eng.spotFinder2(MFA + 'gridFit/blue_window200.png', 1))
        # centroids_2nd = spotFinder(MFA + 'gridFit/2nd_detection_cropped_window200.png')
        # # centroids_curated = np.array(eng.filterOutMarks(MFA + '/gridFit/ablation_marks_XY_2ndFT.npy',
        # #                          MFA + 'gridFit/blue_window200.png')).T
        # centroids_2nd = np.reshape(centroids_2nd, np.shape(centroids_2nd))
        # dataPoints = np.zeros([i for i in reversed(np.shape(centroids_2nd))])
        # dataPoints[0,:] = centroids_2nd[:, 0] + minY
        # dataPoints[1,:] = centroids_2nd[:, 1] + minX
        #
        # # dataPoints = np.array(centroids).T

        centroids = spotFinder(MF,
                               path1=MFA + 'gridFit/AM_cropped.tif',
                               path2=MFA + 'gridFit/AM_cropped_2.tif',
                               matrix=matrix,
                               show_results=True)  # TODO provide option to swtch to DAN at higher level

        crop = np.load(MFA + 'gridFit/AM_cropped_coords.npy', allow_pickle=True).item()
        minX = crop['topLeft'][0]
        minY = crop['topLeft'][1]
        dataPoints = np.zeros([i for i in reversed(np.shape(centroids))])
        dataPoints[0, :] = centroids[:, 0] + minX
        dataPoints[1, :] = centroids[:, 1] + minY
        np.save(MFA + 'gridFit/ablation_marks_XY.npy', dataPoints)

    dataPoints = np.load(MFA + 'gridFit/ablation_marks_XY.npy', allow_pickle=True)
    shape, resX, resY = getShapeRes(MFI)
    pix_size = getPixSize(MFI)
    xe = dataPoints[1]*pix_size
    ye = dataPoints[0]*pix_size
    np.save(MFA + '/gridFit/metadata.npy', [shape, pix_size])
    rotation_esti = estimateAngle(xe, ye, shape, MFA, Figures = True)
    print('Angle estimation finished')  # \nEstimated rotation = {}deg\nData cleaning started ...'.format(rotation_esti)
    xe_r, ye_r, xr_clean, yr_clean, xe_clean, ye_clean = cleanData(rotation=rotation_esti,
                                                                   resX=resX,
                                                                   resY=resY,
                                                                   xe=xe,
                                                                   ye=ye,
                                                                   tolerance=100,
                                                                   manual_cleaning=True,
                                                                   Figures=True)
    print('Data cleaning finished')  # \nCenter estimation started ...'
    center_x_esti, center_y_esti, top, bottom, right, left = estimateCenter(shape, xr_clean, yr_clean, MFA, Figures = True)
    print('Center estimation finished')  # \nCenter X = {}\nCenter Y = {}\nLattice Estimation started ...'.format(center_x_esti, center_y_esti)
    #
    # center_x_esti = 12620
    # center_y_esti = 5368

    lattice_esti = estimateLattice(shape, 0, center_x_esti, center_y_esti, top, bottom, right, left, xr_clean, yr_clean, MFA, Figures = True)
    # lattice_esti = 49.3
    # if lattice_esti < resX - 0.5 or lattice_esti > resX + 0.5:
    #     lattice_esti = resX

    print('Lattice estimation finished')#\n Lattice  = {}um'.format(lattice_esti)
    np.save(MFA + '/gridFit/estimatedParams.npy', [rotation_esti, lattice_esti, center_x_esti, center_y_esti])
    print('Estimated parameters are \nLattice = {}\nRotation = {}\nCenterX = {}\nCenterY = {}'\
        .format(lattice_esti, rotation_esti, center_x_esti, center_y_esti))

    xr_spots, yr_spots = create_grid(shape, 0, lattice_esti, center_x_esti, center_y_esti)
    xe_spots, ye_spots = create_grid(shape, rotation_esti, lattice_esti, center_x_esti, center_y_esti)

    if optimization==True:
        minF = basinhopping(err_func, x0=np.array([rotation_esti, lattice_esti, center_x_esti, center_y_esti]),
                            niter=2, T=1.0, stepsize=0.5,
                            minimizer_kwargs={'args': ((shape, xr_clean, yr_clean))},
                            take_step=None, accept_test=None, callback=None, interval=50, disp=True,
                            niter_success=1)
        # print(minF)
        np.save(MFA + '/gridFit/optimizedParam.npy', [minF.x[0], minF.x[1], minF.x[2], minF.x[3]])
        np.save(MFA + '/gridFit/optimizationResults.npy', minF)
        rotation_opti, lattice_opti, center_x_opti, center_y_opti = np.load(MFA + '/gridFit/optimizedParam.npy')
        print('Optimization done. Parameters are now: \nLattice = {} --> {}\n'
              'Rotation = {} --> {}\n'
              'CenterX = {} --> {}\n'
              'CenterY = {} --> {}' \
              .format(lattice_esti, lattice_opti,
                      rotation_esti, rotation_opti,
                      center_x_esti, center_x_opti,
                      center_y_esti, center_y_opti))
        xr_spots, yr_spots = create_grid(shape, 0, lattice_opti, center_x_opti, center_y_opti)
        xe_spots, ye_spots = create_grid(shape, rotation_esti, lattice_opti, center_x_opti, center_y_opti)

    xr_spots = xr_spots/pix_size
    yr_spots = yr_spots/pix_size
    xe_spots = xe_spots/pix_size
    ye_spots = ye_spots/pix_size

    distances, coord = get_distance(xr_spots, yr_spots, xr_clean/pix_size, yr_clean/pix_size,1)
    xr_clean2 = xr_clean[coord]/pix_size
    yr_clean2 = yr_clean[coord]/pix_size

    xe_clean2 = xe_clean[coord]/pix_size
    ye_clean2 = ye_clean[coord]/pix_size

    crop = np.load(MFA + 'gridFit/AM_cropped_coords.npy', allow_pickle=True).item()
    minX = crop['topLeft'][0]
    minY = crop['topLeft'][1]
    xe_scaled = xe / pix_size - minY
    ye_scaled = ye / pix_size - minX
    xe_clean2_scaled = xe_clean2 - minY
    ye_clean2_scaled = ye_clean2 - minX
    xe_spots_scaled = xe_spots - minY
    ye_spots_scaled = ye_spots - minX

    plt.figure(figsize=(40, 40))
    plt.imshow(plt.imread(MFA + 'gridFit/AM_cropped.tif'), cmap='gray')
    plt.scatter(ye_scaled, xe_scaled, 5, label='Input')
    plt.scatter(ye_clean2_scaled, xe_clean2_scaled, 5, color=[0, 1, 0], label='Output - cleaned')
    plt.scatter(ye_spots_scaled, xe_spots_scaled, 5, color=[1, 0, 0], label='Output - cleaned')

    plt.savefig(MFA + '/gridFit/GridFit_results.png', dpi=100)
    plt.tight_layout()
    plt.close('all')

    np.save(MFA + '/gridFit/xye_clean2.npy', [xe_clean2, ye_clean2])
    np.save(MFA + '/gridFit/xyr_clean2.npy', [xr_clean2, yr_clean2])
    np.save(MFA + '/gridFit/xye_grid.npy', [xe_spots, ye_spots])
    np.save(MFA + '/gridFit/xyr_grid.npy', [xr_spots, yr_spots])


def regionGrowing(cIM, initPos, thresVal, maxDist):


    regVal = cIM[initPos[1], initPos[0]]
    nRow, nCol = np.shape(cIM)
    J = np.zeros([nRow, nCol], dtype=bool)
    queue = [[initPos[1], initPos[0]]]

    while len(queue) > 0:
        xv = queue[0][0]
        yv = queue[0][1]

        del queue[0]

        for i in [-1,0,1]:
            for j in [-1,0,1]:
                if xv+i>=0 and xv+1 <= nRow-1 and \
                        yv+j>=0 and yv+j <=nCol-1  and \
                        any([i, j]) and \
                        ~J[xv+i, yv+j] and \
                        np.sqrt((xv+i - initPos[1])**2 + (yv+j - initPos[0])**2) < maxDist and \
                        cIM[xv+i, yv+j] <= (regVal + thresVal*regVal) and cIM[xv+i, yv+j] >= (regVal - thresVal*regVal):
                            J[xv+i, yv+j] = True
                            queue.append([xv+i, yv+j])

    x, y = np.where(J == True)

    return x,y


def regionGrowingAblationMarks(MF, grow_thresh=0.35, blur=False, sigma=3, matrix='DAN', refine_seed=True, FT_filtered=False, maxDist=15):

    MFA = MF + 'Analysis/'
    # MFA = 'E:\Experiments\TNFa_2.3_SELECTED\Analysis/'
    marks = np.load(MFA + 'gridFit/xye_clean2.npy')
    window = 100

    if matrix == 'DAN':
        mark_check_p = MFA + 'gridFit/marks_check/PHASE_crop_bin1x1_window100.tiff'
        img0 = scale(plt.imread(mark_check_p))
        img_rgb = np.repeat(img0[:, :, np.newaxis], 3, axis=2)
        cut_thresh = 0.1 #was 0.8 before 20181126

    if matrix == 'DHB':
        mark_check_p = MFA + 'gridFit/AM_segmentation_source.tiff'
        path1 = MFA + 'gridFit/marks_check/Fluo_crop_bin1x1_window100.tiff'
        path2 = MFA + 'gridFit/marks_check/PHASE_crop_bin1x1_window100.tiff'
        img_i = DHB_prep(MF, path1, path2)
        img0 = scale(contrast(scale(ndimage.gaussian_filter(img_i, sigma=3)), 0, 0.99))
        img_results = scale(plt.imread(path2))
        img_rgb = np.repeat(img_results[:, :, np.newaxis], 3, axis=2)
        cut_thresh = 0.1

    if FT_filtered:
        big_FT = np.zeros(np.shape(plt.imread(MFA + 'StitchedMicroscopy/postMALDI_FLR/img_t1_z1_c1')))
        crop_coords = np.load(MFA + 'gridFit/AM_cropped_cropCoords.npy').item()
        small_FT = plt.imread(MFA + 'gridFit/FT_filtered.tiff')

        for key in crop_coords.keys():
            crop_coords[key] = [int(np.round(item)) for item in crop_coords[key]]

        big_FT[crop_coords['topLeft'][1] : crop_coords['bottomLeft'][1],
               crop_coords['bottomLeft'][0] : crop_coords['bottomRight'][0]] = small_FT

        tiff.imsave(MF + 'Analysis/gridFit/big_FT.tiff', big_FT.astype(np.uint16))
        spaceM.ImageFileManipulation.manipulations.crop2coords(
            MF + 'Analysis/gridFit/xye_clean2.npy',
            MF + 'Analysis/gridFit/big_FT.tiff',
            MF + 'Analysis/gridFit/marks_check/FT_crop_bin1x1_window100.tiff',
            window=100)
        img0 = scale(plt.imread(MF + 'Analysis/gridFit/marks_check/FT_crop_bin1x1_window100.tiff'))

    cut_window = 50
    marks_s = np.zeros(np.shape(marks))
    marks_s[0, :] = marks[0, :] - np.min(marks[0, :]) + window
    marks_s[1, :] = marks[1, :] - np.min(marks[1, :]) + window

    marksMask = []
    if blur:
        img = ndimage.gaussian_filter(img0, sigma=sigma)
    else:
        img = img0

    # img_rgb = np.repeat(img[:, :, np.newaxis], 3, axis=2)

    n_AM = np.shape(marks_s)[1]
    # n_AM = 300
    for i in range(n_AM):
        im_cut = img[int(marks_s[0, i]) - cut_window : int(marks_s[0, i]) + cut_window,
                 int(marks_s[1, i]) - cut_window:int(marks_s[1, i]) + cut_window]

        if refine_seed:
            thresh = np.percentile(im_cut, 96.5)
            posX, posY = np.where(im_cut > thresh)
            if len(posX) == 0:
                posX, posY = np.where(im_cut > cut_thresh) #TODO 0.8 for DAN
            tree = spatial.KDTree(list(zip(posX.ravel(), posY.ravel())))
            pts = np.array([[cut_window, cut_window]])
            dist, ind = tree.query(pts)
            seed = [posY[ind][0], posX[ind][0]]
        else:
            seed = [cut_window, cut_window]

        xi,yi = regionGrowing(cIM=im_cut,
                      initPos=seed,
                      thresVal=grow_thresh,
                      maxDist=maxDist)

        x = xi + marks[0, i] - cut_window
        y = yi + marks[1, i] - cut_window

        # marksMask[str(i)]['x'] = x
        # marksMask[str(i)]['y'] = y

        marksMask.append([x, y])

        x_int = [int(k) for k in xi + marks_s[0, i] - cut_window]
        y_int = [int(k) for k in yi + marks_s[1, i] - cut_window]

        img_rgb[x_int, y_int, 0] = 0
        img_rgb[x_int, y_int, 1] = 1
        img_rgb[x_int, y_int, 2] = 0

        width = 2

        if refine_seed:
            x_point = int(posX[ind][0] + marks_s[0, i] - cut_window)
            y_point = int(posY[ind][0] + marks_s[1, i] - cut_window)
            img_rgb[x_point - width:x_point + width, y_point - width:y_point + width, 0] = 1
            img_rgb[x_point - width:x_point + width, y_point - width:y_point + width, 1] = 0
            img_rgb[x_point - width:x_point + width, y_point - width:y_point + width, 2] = 0

        x_point_i = int(marks_s[0, i])
        y_point_i = int(marks_s[1, i])
        img_rgb[x_point_i - width:x_point_i + width, y_point_i - width:y_point_i + width, 0] = 0
        img_rgb[x_point_i - width:x_point_i + width, y_point_i - width:y_point_i + width, 1] = 0
        img_rgb[x_point_i - width:x_point_i + width, y_point_i - width:y_point_i + width, 2] = 1

    # plt.imshow(img_rgb)

    plt.imsave(MFA + 'gridFit/AM_segmentation.png', scale(img_rgb).astype(np.float))
    np.save(MFA + 'gridFit/marksMask.npy', marksMask)


def AM_filter(MF, n_std=4, manual_threshold=None):

    def gaus(x, a, x0, sigma):
        return a * exp(-(x - x0) ** 2 / (2 * sigma ** 2))

    MFA = MF + 'Analysis/'
    MFI = MF + 'Input/'
    pix_size = getPixSize(MFI)
    img_rgb_filter = plt.imread(MFA + 'gridFit/AM_segmentation.png')
    marksMask = np.load(MFA + 'gridFit/marksMask.npy', allow_pickle=True)
    marks = np.load(MFA + 'gridFit/xye_clean2.npy', allow_pickle=True)
    window = 100
    marks_s = np.zeros(np.shape(marks))
    marks_s[0, :] = marks[0, :] - np.min(marks[0, :]) + window
    marks_s[1, :] = marks[1, :] - np.min(marks[1, :]) + window
    n_AM = len(marksMask)

    # img_rgb_filter = np.copy(img_rgb)
    areas = np.array([len(x) for x, y in marksMask])
    indexes = np.arange(n_AM)

    y, x = np.histogram(areas, 100)
    x = x[1:]
    mean = np.mean(areas)
    sigma = np.std(areas)
    if not manual_threshold:
        popt, pcov = curve_fit(gaus, x, y, p0=[1, mean, sigma])
        AM_threshold_up = popt[1] + n_std * abs(popt[2])
        AM_threshold_down = popt[1] - n_std * abs(popt[2])
    else:
        AM_threshold_down = manual_threshold[0]
        AM_threshold_up = manual_threshold[1]

    AM_pass_indexes = []
    for i in indexes:
        if areas[i] < AM_threshold_up and areas[i] > AM_threshold_down:
            AM_pass_indexes.append(int(i))
    # AM_pass_indexes = [int(i) for i in indexes[areas < AM_threshold]]

    for i in range(n_AM):
        x_all, y_all = marksMask[i]

        x_int = [int(x) for x in x_all - marks[0, i] + marks_s[0, i]]
        y_int = [int(y) for y in y_all - marks[1, i] + marks_s[1, i]]

        if i not in AM_pass_indexes:
            img_rgb_filter[x_int, y_int, 0] = 1
            img_rgb_filter[x_int, y_int, 1] = 0
            img_rgb_filter[x_int, y_int, 2] = 0

    plt.imsave(MFA + 'gridFit/AM_segmentation_filtered.png', scale(img_rgb_filter).astype(np.float))
    np.save(MFA + 'gridFit/AM_pass_filter.npy', AM_pass_indexes)

    plt.figure()
    plt.plot(x*pix_size, y, label='data')
    if not manual_threshold:
        plt.plot(x*pix_size, gaus(x, *popt), label='fit')
    plt.axvline(AM_threshold_up*pix_size, label='threshold')
    plt.axvline(AM_threshold_down * pix_size, label='threshold')
    plt.legend()
    plt.xlabel('AM area (um**2)', fontsize=15)
    plt.ylabel('Counts', fontsize=15)
    plt.savefig(MF + 'Analysis/gridFit/AM_filter.png')
    plt.close('all')

def marksSegmentedMask(MFA):
    """Segment ablation marks using their center of mass as seed for a region growing algorithm.
    This function uses the open-source implementation of the region growing by Daniel Kellner
    available from the MatlabExchange (https://de.mathworks.com/matlabcentral/fileexchange/32532-region-growing--2d-3d-grayscale).

    Args:
        MFA (str): path to Main Folder Analysis.
    """
    eng = matlab.engine.start_matlab()
    window=100
    dummy = eng.regionGrowAblationMarks(MFA + 'gridFit/marks_check/PHASE_crop_bin1x1_window100.png',
                                            MFA + 'gridFit/xye_clean2.npy', window,
                                            MFA + 'gridFit/marksMask.mat')
    eng.quit()
    mat = scipy.io.loadmat(MFA + 'gridFit/marksMask.mat')
    marksMask = mat['out'][0]
    np.save(MFA + 'gridFit/marksMask.npy', marksMask)

    img = plt.imread(MFA + 'gridFit/marks_check/PHASE_crop_bin1x1_window100.png')
    coordX, coordY = np.load(MFA + 'gridFit/xye_clean2.npy')

    plt.figure()
    plt.switch_backend('TkAgg')
    plt.get_backend()
    mng = plt.get_current_fig_manager()
    mng.window.state('zoomed')
    sns.set_style("whitegrid", {'axes.grid': False})

    plt.imshow(img, cmap='gray')
    plt.scatter(coordY- np.min(coordY) + window, coordX- np.min(coordX) + window, 0.005, facecolors='g')
    plt.show()
    plt.pause(0.05)

    for i in tqdm.tqdm(range(np.shape(marksMask)[0])):
        if np.shape(marksMask[i][0].T)[1] > 1:
            plt.scatter(marksMask[i][1].T- np.min(coordY) + window, marksMask[i][0].T- np.min(coordX) + window, 0.5, facecolors='r')
            # tfMarksMask.append([tfMask[:, 0], tfMask[:, 1]])
    plt.savefig(MFA + 'gridFit/marksMask_vizu.png', dpi=1200, bbox_inches='tight', frameon='false')
    plt.close('all')









