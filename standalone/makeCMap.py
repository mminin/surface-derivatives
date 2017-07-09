#This snippet makes a very simple colourmap for use in FD.py

import numpy as np
from osgeo import gdal
from scipy import ndimage
from PIL import Image
import matplotlib.pyplot as plt

q=np.linspace(0,255,60)
a=np.asarray([np.linspace(0,w,90) for w in q]).astype(int)
b=np.asarray([np.linspace(0,255,90) for w in q]).astype(int)
c=a[::-1].astype(int)
d=np.asarray([np.linspace(0,0,90) for w in q]).astype(int)

R=np.concatenate([b,c,d,d,a,b])
G=np.concatenate([a,b,b,c,d,d])
B=np.concatenate([d,d,a,b,b,c])

Comp=np.stack((R,G,B)).swapaxes(0,2)
CompB=np.uint8(Comp)

Image.fromarray(CompB).save('CMap.tif')
