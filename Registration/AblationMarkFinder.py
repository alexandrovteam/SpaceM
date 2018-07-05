import tifffile as tiff
import numpy as np
import glob, os
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



def spotFinder(path, layer=3):
    """Detect ablation marks on the tiled images. For details, see the paper.
    Args:
        path (str): path of the image to detect ablation marks on.

    """

    def scale(arr):
        """Scale array between 0 and 1"""
        return (arr - np.min(arr)) / (np.max(arr) - np.min(arr))

    def contrast(arr, min, max=1):
        """Clip array between min and max values"""
        return np.clip(arr, min, max)

    # img_i = plt.imread(path)
    im = tiff.imread(path)
    if len(np.shape(im)) > 2:
        img_i = im[0, :, :]
    img = scale(img_i)
    contrast_min = np.mean(img) + 2*np.std(img)
    if contrast_min >=1: contrast_min=0.8
    img_2 = contrast(img, min=contrast_min)
    # plt.imshow(img_2, cmap='gray')

    ff = np.fft.fft2(img_2)
    F1 = 20*np.log10(np.abs(np.fft.fftshift(ff)))
    # plt.imshow(F1, cmap='gray')

    mask1 = F1 - ndimage.gaussian_filter(F1, sigma=15)
    mask1[mask1<0] = 0
    # plt.imshow(mask1, cmap='gray')

    int_thresh = np.mean(mask1 + 3.5*np.std(mask1))
    struct = ndimage.generate_binary_structure(2,1)
    mask2 = ndimage.binary_dilation(mask1>int_thresh, structure=struct, iterations=20).astype(mask1.dtype)
    # plt.imshow(mask2, cmap='gray')

    ff_masked = np.fft.fftshift(mask2)*ff
    # F2 = 20 * np.log10(abs(np.fft.fftshift(ff_masked))+1)
    # plt.imshow(F2, cmap='gray')

    freq_up = 0.6
    freq_down = 0.0
    [N, M] = np.shape(img)
    dx = 1
    KX0 = (np.mod(1 / 2 + np.arange(0,M) / M, 1) - 1 / 2)
    KX1 = KX0 * (2 * np.pi / dx)
    KY0 = (np.mod(1 / 2 + np.arange(0,N) / N, 1) - 1 / 2)
    KY1 = KY0 * (2 * np.pi / dx)
    [KX, KY] = np.meshgrid(KX1, KY1)
    lpf = (KX * KX + KY * KY < freq_up ** 2)
    lpf2 = (KX * KX + KY * KY < freq_down ** 2)
    mix = lpf * ~lpf2
    rec = np.real(np.fft.ifft2(mix * ff_masked))
    rec = scale(rec)
    rec = contrast(rec, np.percentile(rec, 1), np.percentile(rec, 99))
    # plt.imshow(rec)

    rec_bw = rec > np.mean(rec) + 2.5*np.std(rec)
    label_img = label(rec_bw, connectivity=rec_bw.ndim)
    props = regionprops(label_img)
    centroids = [[props[i].centroid[1],props[i].centroid[0]]  for i in range(len(props))]

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


def GridFit(MF, optimization=False, manual_cleaning=True, MarkFinderFT=True):
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
        rotation_deg, affine_lat, i, j = params
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

    def getPixSize(MFI):
        """Reads the pixel size in um from the Nikon Ti E microscope (NIS elements software).

        Args:
            MFI (str): path to Main Folder Input.

        Returns:
            pix_size (float): pixel size in um.

        """
        txt_file = codecs.open(MFI + '/Microscopy/postMALDI/out.txt', 'r','utf-16')
        for row in txt_file:
            if row.startswith('Calibration'):
                pix_size = float(row.strip().split()[2].replace(',', '.'))
            else:
                pix_size = 0.73
        return pix_size

    def estimateAngle(xe,ye, shape, MFA, Figures=True):
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
        rotation_deg_y = angle[np.where(counts_y == np.min(counts_y))[0][0]]
        rotation_deg = np.mean([rotation_deg_x,rotation_deg_y])
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
        return rotation_deg

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
        if manual_cleaning==True:
            return xe_r, ye_r, xe, ye
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
            # pics_ind_clean = np.delete(pics_ind, ind_cut)
            if Figures == True:
                plt.figure()
                x_line, = plt.plot(b_x[:-1], c_x, 'r', label = 'X_projections')
                thresh_x_line, = plt.plot(b_x[:-1], np.ones(np.shape(b_x[:-1])) * treshc_x, 'r', label='Threshold X')
                y_line, = plt.plot(b_y[:-1], c_y, 'k', label = 'Y_projections')
                thresh_y_line, = plt.plot(b_y[:-1], np.ones(np.shape(b_y[:-1])) * treshc_y, 'k', label='Threshold Y')
                plt.legend(handles = [x_line, y_line, thresh_x_line, thresh_y_line])#, ['X projections', 'Y projections'])

                plt.savefig(MFA + '/gridFit/XYprojections.png')
            return xr_clean, yr_clean, xe_clean, ye_clean
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
            plt.scatter(xr_clean, yr_clean, 10, 'k')
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
            plt.plot(steps, dist, 'ko')
            plt.xlabel('Lattice size (um)', fontsize=20)
            plt.ylabel('Mean closest neighbor from\nTheoretical to Experimental grid', fontsize=20)
            plt.title('Finer lattice estimation', fontsize=25)
            plt.scatter(affine_lat1, np.min(dist), 200, 'r', label='Global minimum')
            plt.legend()
            plt.savefig(MFA + '/gridFit/estimateLattice1.png')
            plt.close('all')
        # return affine_lat

        dist = []
        steps = []
        for i in np.linspace(affine_lat1 - 0.15, affine_lat1 + 0.15, nsteps):
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
            plt.plot(steps, dist, 'ko')
            plt.xlabel('Lattice size (um)', fontsize=20)
            plt.ylabel('Mean closest neighbor from\nTheoretical to Experimental grid', fontsize=20)
            plt.title('Finer lattice estimation', fontsize=25)
            plt.scatter(affine_lat2, np.min(dist), 200, 'r', label='Global minimum')
            plt.legend()
            plt.savefig(MFA + '/gridFit/estimateLattice2.png')
            plt.close('all')
        return affine_lat2


    MFI = MF + 'Input/'
    MFA = MF + 'Analysis/'
    dataPoints = np.load(MFA + 'gridFit/ablation_marks_XY.npy')

    if MarkFinderFT:
        # Ablation mark detection efficiency can be affected by the number of ablation marks present in the image.
        # To improve the results, the same algorithm is run a seond time on the cropped stitched microsocpy around the previous detections.
        # That way all ablation marks are present in the same image, maximizing their associated frequencies in the fourrier domain.

        window = 200
        if not os.path.exists(MFA + 'gridFit/SpotFinder_cropped_window200.png'):
            manip.crop2coords(MFA + 'gridFit/ablation_marks_XY.npy',
                              MFA + 'StitchedMicroscopy/postMALDI_FLR/img_t1_z1_c1',
                              MFA + 'gridFit/2nd_detection_cropped_window200.png', window=window)

        # if not os.path.exists(MFA + 'gridFit/ablation_marks_XY_2nd_detection.npy'):
        X, Y = np.load(MFA + 'gridFit/ablation_marks_XY.npy')
        minX = np.min(X) - window
        minY = np.min(Y) - window
        # eng = matlab.engine.start_matlab()
        # centroids = np.asarray(eng.spotFinder2(MFA + 'gridFit/blue_window200.png', 1))
        centroids_2nd = spotFinder(MFA + 'gridFit/2nd_detection_cropped_window200.png')
        # centroids_curated = np.array(eng.filterOutMarks(MFA + '/gridFit/ablation_marks_XY_2ndFT.npy',
        #                          MFA + 'gridFit/blue_window200.png')).T
        centroids_2nd = np.reshape(centroids_2nd, np.shape(centroids_2nd))
        dataPoints = np.zeros([i for i in reversed(np.shape(centroids_2nd))])
        dataPoints[0,:] = centroids_2nd[:, 0] + minY
        dataPoints[1,:] = centroids_2nd[:, 1] + minX

        # dataPoints = np.array(centroids).T
        np.save(MFA + 'gridFit/ablation_marks_XY_2nd_detection.npy', dataPoints)

        dataPoints = np.load(MFA + 'gridFit/ablation_marks_XY_2nd_detection.npy')
        # eng.quit()

    shape, resX, resY = getShapeRes(MFI)
    pix_size = getPixSize(MFI)
    xe = dataPoints[1]*pix_size
    ye = dataPoints[0]*pix_size
    np.save(MFA + '/gridFit/metadata.npy', [shape, pix_size])
    rotation_esti = estimateAngle(xe, ye, shape, MFA, Figures = True)
    print('Angle estimation finished') #\nEstimated rotation = {}deg\nData cleaning started ...'.format(rotation_esti)
    xr_clean, yr_clean, xe_clean, ye_clean = cleanData(rotation_esti, resX, resY, xe, ye, tolerance = 100,
                                                       manual_cleaning = True, Figures = True)
    print('Data cleaning finished')#\nCenter estimation started ...'
    center_x_esti, center_y_esti, top, bottom, right, left = estimateCenter(shape, xr_clean, yr_clean, MFA, Figures = True)
    print('Center estimation finished')#\nCenter X = {}\nCenter Y = {}\nLattice Estimation started ...'.format(center_x_esti, center_y_esti)
    lattice_esti = estimateLattice(shape, 0, center_x_esti, center_y_esti, top, bottom, right, left, xr_clean, yr_clean, MFA, Figures = True)
    if lattice_esti < resX - 0.5 or lattice_esti > resX + 0.5:
        lattice_esti = resX
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

    # plt.figure()
    # plt.scatter(xe/pix_size, ye/pix_size, 30, 'k')
    # plt.scatter(xe_spots, ye_spots, 30, 'r')
    # plt.axis('equal')
    #
    # plt.figure()
    # plt.scatter(xe_spots, ye_spots, 10, 'r')
    # plt.scatter(xe_clean2, ye_clean2, 10, 'k')
    # plt.axis('equal')

    np.save(MFA + '/gridFit/xye_clean2.npy', [xe_clean2, ye_clean2])
    np.save(MFA + '/gridFit/xyr_clean2.npy', [xr_clean2, yr_clean2])
    np.save(MFA + '/gridFit/xye_grid.npy', [xe_spots, ye_spots])
    np.save(MFA + '/gridFit/xyr_grid.npy', [xr_spots, yr_spots])

    # for i in range(len(xe)):
    #     # print(i/len(xe))
    #     c = plt.cm.jet(i/len(xe))
    #     plt.scatter(xe_clean2[i],ye_clean2[i], 40, c)
    #     plt.axis('equal')


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
                        cIM[xv+i, yv+j] <= (regVal + thresVal) and cIM[xv+i, yv+j] >= (regVal - thresVal):
                            J[xv+i, yv+j] = True
                            queue.append([xv+i, yv+j])

    x, y = np.where(J == True)

    return x,y

def regionGrowingAblationMarks(MFA, grow_thresh=0.2, blur=True, sigma=3):

    # MFA = 'E:\Experiments\TNFa_2.3_SELECTED\Analysis/'
    marks = np.load(MFA + 'gridFit/xye_clean2.npy')
    mark_check_p = MFA + 'gridFit/marks_check/PHASE_crop_bin1x1_window100.png'
    img0 = plt.imread(mark_check_p)
    window=100
    cut_window = 50
    marks_s = np.zeros(np.shape(marks))
    marks_s[0, :] = marks[0, :] - np.min(marks[0, :]) + window
    marks_s[1, :] = marks[1, :] - np.min(marks[1, :]) + window
    img_rgb = np.repeat(img0[:, :, np.newaxis], 3, axis=2)
    # marksMask = {}
    # for i in range(np.shape(marks_s)[1]):
    #     marksMask[str(i)] = {}
    #     marksMask[str(i)]['x'] = []
    #     marksMask[str(i)]['y'] = []
    marksMask = []
    if blur:
        img = ndimage.gaussian_filter(img0, sigma=sigma)
    else:
        img = img0
    for i in tqdm.tqdm(range(np.shape(marks_s)[1])):
        im_cut = img[int(marks_s[0, i]) - cut_window : int(marks_s[0, i]) + cut_window,
                 int(marks_s[1, i]) - cut_window:int(marks_s[1, i]) + cut_window]
        thresh = np.percentile(im_cut, 98)
        posX, posY = np.where(im_cut > thresh)

        if len(posX) == 0:
            posX, posY = np.where(im_cut > 0.8)

        tree = spatial.KDTree(list(zip(posX.ravel(), posY.ravel())))
        pts = np.array([[cut_window, cut_window]])
        dist, ind = tree.query(pts)

        xi,yi = regionGrowing(cIM=im_cut,
                      initPos=[posY[ind][0], posX[ind][0]],
                      thresVal=grow_thresh,
                      maxDist=60)

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

    tiff.imsave(MFA + 'gridFit/AM_segmentation.tiff', img_rgb)
    np.save(MFA + 'gridFit/marksMask.npy',marksMask)

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