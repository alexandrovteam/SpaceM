import os
import codecs
import numpy as np
from subprocess import call
import math

def TileConfFormat(path, dir_fliplr, tif_files):
    """Extract the microscope motor stage coordinates at each frame from the metadata text file from the Nikon
    Ti-E microscope (NIS elements software) and reformat into readable format for the Priebisch software algorithm
    from FIJI.

    Args:
        path (str): path of the directory containing the microscope metadata file (named 'out.txt' by default).
        dir_fliplr (path): path of the directory containing the transformed tiled frames to stitch.
        tif_files (array): names of each tiled frames to stitch.

    """
    if 'out.txt' in os.listdir(path):
        txt_file = codecs.open(path + 'out.txt', 'r','utf-16')
        data = []
        out_file = open(dir_fliplr + 'TileConfiguration.txt','w')
        out_file.write('# Define the number of dimensions we are working on\ndim = 2\n\n# Define the image coordinates\n')
        i = 0
        for row in txt_file:
            if row.startswith('#'):
                #print(row.strip().split('\t'))
                if i <= np.shape(tif_files)[0]-1:
                    data.append(row.strip().split('\t'))
                    data[i][0] = str(i+1).zfill(3)
                    data[i][1] = float(data[i][1].replace(',','.'))
                    data[i][2] = float(data[i][2].replace(',','.'))
                    out_file.write('img_XY{}.tif; ; ({}, {})\n'.format(data[i][0],data[i][1],data[i][2]))
                    # print i
                    i = i+1
            elif row.startswith('Spectral Loop'):
                break
        out_file.close()

def callFIJImergeRedGray(base_path, red_filename, gray_filename,
                         save_filename):
    """Creates a FIJI macro and call FIJI executable to merge two channels (brighfield and Red) stored in independent files.

    Args:
        base_path (str): path of the directory containing the two images to merge.
        red_filename (str): name of the image file for red channel.
        gray_filename (str): name of the image file for gray channel.
        save_filename (str): name of the merged image containing both chanel. Saved as an RGB image

    """
    script_file_p = base_path + 'mergeRedGray_script.txt'
    base = os.path.splitext(script_file_p)[0]
    if not os.path.exists(base_path + 'mergeRedGray_script.ijm' ):
        out_file2 = open(script_file_p , 'w')
        out_file2.write('\
        open("{}")\
        \nopen("{}")\
        \nrun("Merge Channels...", "c1={} c4={} create")\
        \nsaveAs("PNG", "{}");\
        \nrun("Quit");'\
                        .format(base_path.replace('/', '\\\\') + red_filename,
                                base_path.replace('/', '\\\\') + gray_filename,
                                red_filename,
                                gray_filename,
                                base_path.replace('/', '\\\\') + save_filename))
        out_file2.close()

        os.rename(script_file_p, base + ".ijm")
    call(['C:\\Users\Luca\Documents\Fiji.app\ImageJ-win64.old.exe', '-macro', base.replace('/', '\\') + ".ijm"])#, stdout = PIPE)

def callFIJImergeGrayRedBlue(base_path, red_filename, gray_filename, blue_filename,
                         save_filename):
    """Creates a FIJI macro and call FIJI executable to merge three channels (brighfield, red and blue) stored in independent files.

    Args:
        base_path (str): path of the directory containing the two images to merge.
        red_filename (str): name of the image file for red channel.
        gray_filename (str): name of the image file for gray channel.
        blue_filename (str): name of the image file for blue channel.
        save_filename (str): name of the merged image containing both chanel. Saved as an RGB image

    """
    script_file_p = base_path + 'mergeRedGrayBlue_script.txt'
    base = os.path.splitext(script_file_p)[0]
    if not os.path.exists(base_path + 'mergeRedGrayBlue_script.ijm' ):
        out_file2 = open(script_file_p , 'w')
        out_file2.write('\
        open("{}")\
        \nopen("{}")\
        \nopen("{}")\
        \nrun("Merge Channels...", "c1={} c4={} c3={} create")\
        \nsaveAs("PNG", "{}");\
        \nrun("Quit");'\
                        .format(base_path.replace('/', '\\\\') + red_filename,
                                base_path.replace('/', '\\\\') + gray_filename,
                                base_path.replace('/', '\\\\') + blue_filename,
                                red_filename,
                                gray_filename,
                                blue_filename,
                                base_path.replace('/', '\\\\') + save_filename))
        out_file2.close()

        os.rename(script_file_p, base + ".ijm")
    call(['C:\\Users\Luca\Documents\Fiji.app\ImageJ-win64.old.exe', '-macro', base.replace('/', '\\') + ".ijm"])#, stdout = PIPE)

def callFIJIstitch(dir_fliplr):
    """Calls FIJI stitching algorithm on the transformed tiled frames using the reformatted metadata text file.
    The function creates a FIJI macro which is then ran by FIJI application.

    Args:
        dir_fliplr (str): path of the directory containing the tiled frames to stitch and the metadata text file.

    """
    script_file_p = dir_fliplr + 'stitch_script.txt'
    out_file2 = open(script_file_p , 'w')
    out_file2.write('run("Grid/Collection stitching", \
    "type=[Positions from file] order=[Defined by TileConfiguration] directory={} \
    layout_file=TileConfiguration.txt fusion_method=[Linear Blending] regression_threshold=0.30 \
    max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 \
    compute_overlap computation_parameters=[Save computation time (but use more RAM)] \
    image_output=[Write to disk] output_directory={}"); \
                    \nrun("Quit");'
                    .format(dir_fliplr.replace('/', '\\\\'),
                            dir_fliplr.replace('/', '\\\\')))
    out_file2.close()
    base = os.path.splitext(script_file_p)[0]
    os.rename(script_file_p, base + ".ijm")
    call(['C:\\Users\Luca\Documents\Fiji.app\ImageJ-win64.old.exe', '-macro', base.replace('/', '\\') + ".ijm"])#, stdout = PIPE)
    # os.remove('C:\\Users\Luca\AppData\Local\Temp\org.scijava.jython.shaded.jline_2_5_3.dll')

def readTileConfReg(dir_fliplr):
    """Extract registered tile image coordinates from the textfile generated by Priebisch stitching algorithm
    in FIJI.

    Args:
        dir_fliplr(str): path of the directory containing the registered coordinates.

    Returns:
        dataXscaled (array): registred tile images X coordinates.
        dataYscaled (array): registred tile images Y coordinates.

    """
    txt_file = codecs.open(dir_fliplr + 'TileConfiguration.registered.txt', 'r')
    dataX = []
    dataY = []
    for row in txt_file:
        if row.startswith('img'):
            dataY = np.append(dataY, float(row.split()[2].replace('(', '').replace(',', '')))
            dataX = np.append(dataX, float(row.split()[3].replace(')', '')))
    dataXscaled = dataX - np.min(dataX)
    dataYscaled = dataY - np.min(dataY)

    return dataXscaled, dataYscaled

def imbin4ili(file_p, maxsize):
    """Bin the images to make its size inferior to a given value. This is required for large images which have to
    be visualized in ili' which has an upper limit of 50MB.

    Args:
        file_p (str): path of the image to bin.
        maxsize (float): upper size threshold in bytes
    Returns:
        bin (int): defined bin factor.

    """
    def round_up_to_even(f):
        return math.ceil(f / 2.) * 2

    dirname = os.path.dirname(file_p) + '/'
    img_size = os.path.getsize(file_p)
    scale = np.sqrt(img_size/maxsize)
    if scale>1:
        bin = round_up_to_even(scale)

    # if bin > 0:
        script_file_p = dirname + 'bin_script.txt'
        out_file2 = open(script_file_p, 'w')
        out_file2.write('\
        open("{}");\
        \nrun("Bin...", "x={} y={} bin=Average");\
        \nsaveAs("PNG", "{}");\
        \nrun("Quit");' \
                        .format((file_p).replace('/', '\\\\'),
                                str(bin), str(bin),
                                (dirname).replace('/', '\\\\') + \
                                'FLUO_crop_bin' + str(bin) + 'x' + str(bin)+ '.png'))
        out_file2.close()
        base = os.path.splitext(script_file_p)[0]
        os.rename(script_file_p, base + ".ijm")
        call(['C:\\Users\Luca\Documents\Fiji.app\ImageJ-win64.exe', '-macro',
              base.replace('/', '\\') + ".ijm"])  # , stdout = PIPE)
    else: bin = 1
    return bin
