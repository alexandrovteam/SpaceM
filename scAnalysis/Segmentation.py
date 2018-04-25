from subprocess import call

def callCP(MFA):
    """Call CellProfiler (http://cellprofiler.org/) to perform cell segmentation. CellProfiler segmentation pipeline
    is in the spaceM folder with the '.cppipe' extension.

     Args:
         MFA (str): path to Main Folder Analysis.

     """
    # CP headless info https://github.com/CellProfiler/CellProfiler/wiki/Adapting-CellProfiler-to-a-LIMS-environment
    file = open(MFA + 'CellProfilerAnalysis\input_files.txt', 'w')
    file.write(MFA + 'CellProfilerAnalysis\img_t1_z1_c1.tif\n')
    file.write(MFA + 'CellProfilerAnalysis\img_t1_z1_c2.tif\n')
    file.write(MFA + 'CellProfilerAnalysis\img_t1_z1_c2_adjusted.tif\n')
    file.write(MFA + 'CellProfilerAnalysis\img_t1_z1_c3.tif')
    file.close()
    call(['C:\\Program Files\CellProfiler\CellProfiler.exe', '-r', '-c', '-p',
          'C:\\Users\Luca\Documents\python_codebase\pyCLIMS\Hepatocytes_segmentation.cppipe',
          '-o', MFA + 'CellProfilerAnalysis\\', '--file-list', MFA + 'CellProfilerAnalysis\input_files.txt'])
