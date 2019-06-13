import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import os
from functools import reduce

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
# MF = 'D:/Experiments/Hepa_June/LPS_1.1/'
MFA = MF + 'Analysis/'

Path = MFA + 'gridFit/ablation_marks_XY'
Path = MFA + 'Fiducials/postXYpenmarks'
Path = MFA + 'Fiducials/preXYpenmarks'

#CHOOSE
purpose = 'select'
purpose = 'delete'

oldX, oldY = np.load(Path + '.npy')


plt.scatter(oldX, oldY, 5)
plt.axis('equal')
a = Annotate()

#draw rectangle from upper left corner to the lower right one
#Close the image


#once image is closed, execute the rest of the script
if purpose == 'select':
    ind = np.unique(np.hstack([np.where(oldX < a.x0)[0], np.where(oldX > a.x1)[0], np.where(oldY > a.y0)[0], np.where(oldY < a.y1)[0]]))
    newX = np.delete(oldX, ind)
    newY = np.delete(oldY, ind)
elif purpose == 'delete':
    ind = reduce(np.intersect1d,(np.where(oldX > a.x0)[0], np.where(oldX < a.x1)[0], np.where(oldY < a.y0)[0], np.where(oldY > a.y1)[0]))
    newX = np.delete(oldX, ind)
    newY = np.delete(oldY, ind)

plt.scatter(newX, newY,5)
plt.axis('equal')

oldX = newX
oldY = newY
a = Annotate()

os.rename(Path + '.npy', Path + 'OLD.npy')
np.save(Path + '.npy', [newX, newY])
plt.scatter(newX, newY, 1)
plt.axis('equal')
plt.savefig(Path + '_manualCLEAN.png', dpi = 500)
plt.close('all')

