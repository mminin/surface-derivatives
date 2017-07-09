#import scipy
from PIL import Image

import numpy as np
from scipy import ndimage
from osgeo import gdal
radToDegFactor = 57.2958

sourceFileName = 'Mikhail_test.tif'
kx=3
ky=3
spc=50

khx, khy=(kx-1)/2.,(ky-1)/2.
mx,my=np.meshgrid(np.linspace(-khx,khx,kx),np.linspace(-khy,khy,ky))

#Apply scaling to kernel, see M. Thesis for explanation:
mx/= spc*kx*khy*(khy+1)*ky 
my/= spc*ky*khx*(khx+1)*kx

sourceImg = gdal.Open(sourceFileName, gdal.GA_ReadOnly)
array = sourceImg.GetRasterBand(1).ReadAsArray()*1.0
FDx, FDy = ndimage.convolve(array,mx),ndimage.convolve(array,my)

Image.fromarray(FDx).save('FDx.tif')
Image.fromarray(FDy).save('FDy.tif')

Strike = np.arctan2(-FDy,FDx)*radToDegFactor
Dip = np.arctan(np.sqrt(FDx*FDx+FDy*FDy))*radToDegFactor

Image.fromarray(Strike).save('Strike.tif')
Image.fromarray(Dip).save('Dip.tif')

## LAST STEP HERE WOULD BE TO MAKE A COLOURMAP PLOTTER

