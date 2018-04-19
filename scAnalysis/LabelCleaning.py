import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import os
from functools import reduce
import tifffile as tiff
class Annotate(object):
    def __init__(self):
        self.ax = plt.gca()
        self.rect = Rectangle((0, 0), 1, 1, facecolor='None')
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.ax.add_patch(self.rect)
        self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)

    def on_press(self, event):
        # print 'press'
        self.x0 = event.xdata
        self.y0 = event.ydata

    def on_release(self, event):
        # print 'release'
        self.x1 = event.xdata
        self.y1 = event.ydata
        self.rect.set_width(self.x1 - self.x0)
        self.rect.set_height(self.y1 - self.y0)
        self.rect.set_xy((self.x0, self.y0))
        self.ax.figure.canvas.draw()
        return [self.x0, self.x1, self.y0, self.y1]

#CHOOSE

MF = 'Z:/rappez/20170622_Hepa_June_DAN_Untreated_FA_LPS_TNFa\Data/transformation_2/LPS_1.2_SELECTED/'

MFA = MF + 'Analysis/'
Path = MFA + 'CellProfilerAnalysis/Labelled_cells'
img = tiff.imread(Path + '.tif')
plt.imshow(img, interpolation = 'none')
print(np.shape(np.unique(img)))
lbl_correct = img
a = Annotate()

indexes = np.unique(img[int(a.y0):int(a.y1), int(a.x0):int(a.x1)])
for ind in indexes:
    lbl_correct[img == ind] = 0.0

plt.imshow(lbl_correct, interpolation = 'none')
a = Annotate()

# lbl_correct = lbl_correct.astype('uint16')

os.rename(Path + '.tif', Path + '_OLD.tif')
tiff.imsave(Path + '.tif', lbl_correct)

