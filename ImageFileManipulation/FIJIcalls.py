import spaceM
import os, re
import codecs
import numpy as np
from subprocess import call
import math
import pandas as pd

fiji_path = pd.read_json(os.path.dirname(spaceM.__file__) + '\\paths.json')['Fiji path'].as_matrix()[0]

def TileConfFormat(path, dir_fliplr, tif_files):
    """Extract the microscope motor stage coordinates at each frame from the metadata text file from the Nikon
    Ti-E microscope (NIS elements software) and reformat into readable format for the Priebisch software algorithm
    from FIJI.

    Args:
        path (str): path of the directory containing the microscope metadata file (named 'out.txt' by default).
        dir_fliplr (str): path of the directory containing the transformed tiled frames to stitch.
        tif_files (array): names of each tiled frames to stitch.

    """
    if 'out.txt' in os.listdir(path):
        txt_file = codecs.open(path + 'out.txt', 'r', 'utf-16')
        data = []
        out_file = open(dir_fliplr + 'TileConfiguration.txt', 'w')
        out_file.write('# Define the number of dimensions we are working on\ndim = 2\n\n# Define the image coordinates\n')
        i = 0
        base = re.findall('^(.*)\d{3}.tif$', tif_files[0])[0]
        for row in txt_file:
            if row.startswith('#'):
                # print(row.strip().split('\t'))
                if i <= np.shape(tif_files)[0]-1:
                    data.append(row.strip().split('\t'))
                    data[i][0] = str(i+1).zfill(3)
                    data[i][1] = float(data[i][1].replace(',','.'))
                    data[i][2] = float(data[i][2].replace(',','.'))
                    out_file.write(base + '{}.tif; ; ({}, {})\n'.format(data[i][0],data[i][1],data[i][2]))
                    # re.findall('^(.*)(\d{3})$', 'seq000_XY120')
                    # print i
                    i = i+1
            elif row.startswith('Spectral Loop'):
                break
        out_file.close()

def callFIJImergeChannels(base_path,
                          colors,
                          filenames,
                          save_filename = 'Composite.png'):

    """Creates a FIJI macro and call FIJI executable to merge different channels stored in independent files into an RGB
    image (.png).

    Args:
        base_path (str): path of the directory containing the images to merge.
        colors (list): list of string of color names: 'red', 'green', 'blue', 'gray', 'cyan', 'magenta', 'yellow'.
        filename (str): list of string of image files names to merge. Their sequence in the list should match their
            respective color in the 'colors' argument.
        save_filename (str): name of the merged image containing all chanels. Saved as an RGB image

    """

    color_codex = {'red': 'c1',
                   'green': 'c2',
                   'blue': 'c3',
                   'gray': 'c4',
                   'cyan': 'c5',
                   'magenta': 'c6',
                   'yellow': 'c7'}

    string1 = ''
    for i in range(len(filenames)):
        string1 = string1 + 'open("{}")\n'.format(
            base_path.replace('/', '\\\\') + filenames[i]
        )

    string2 = ''
    for i in range(len(colors)):
        string2 = string2 + '{}={} '.format(
            color_codex[colors[i]],
            filenames[i])
    string2 = string2 + 'create'

    script_file_p = base_path + 'mergeRedGray_script.txt'
    base = os.path.splitext(script_file_p)[0]

    if not os.path.exists(base_path + 'mergeRedGray_script.ijm' ):
        out_file2 = open(script_file_p , 'w')
        out_file2.write(string1 + 'run("Merge Channels...", "{}")\
        \nsaveAs("PNG", "{}");\
        \nrun("Quit");'\
                        .format(string2,
                                base_path.replace('/', '\\\\') + save_filename))
        out_file2.close()

        if os.path.exists(base + ".ijm"):
            os.remove(base + ".ijm")
        os.rename(script_file_p, base + ".ijm")
    call([fiji_path, '-macro', base.replace('/', '\\') + ".ijm"])#, stdout = PIPE)

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
    if os.path.exists(base + ".ijm"):
        os.remove(base + ".ijm")
    os.rename(script_file_p, base + ".ijm")
    call([fiji_path, '-macro', base.replace('/', '\\') + ".ijm"])#, stdout = PIPE)
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
        call([fiji_path, '-macro',
              base.replace('/', '\\') + ".ijm"])  # , stdout = PIPE)
    else: bin = 1
    return bin
