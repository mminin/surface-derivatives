import sys

sourceFileName = sys.argv[1] # 'Mikhail_test.tif'
if len(sys.argv)>=4: kx,ky,spc = [int(sys.argv[2]),int(sys.argv[3]),float(sys.argv[4])]
else:                kx,ky,spc = [3,3,50]
print(sys.argv[1])
cMapFileName   = 'CMap.tif'
FDxFileName    = 'FDx.tif';    FDyFileName    = 'FDy.tif'
StrikeFileName = 'Strike.tif'; DipFileName    = 'Dip.tif'
outputFileName = 'output.tif'

from PIL import Image
import numpy as np
from scipy import ndimage
from osgeo import gdal
import matplotlib.image as mpimg
radToDegFactor = 57.2958

khx, khy=(kx-1)/2.,(ky-1)/2.
mx,my=np.meshgrid(np.linspace(-khx,khx,kx),np.linspace(-khy,khy,ky))

mx/= spc*kx*khy*(khy+1)*ky # Apply scaling to kernel,
my/= spc*ky*khx*(khx+1)*kx # see M. Thesis for explanation

sourceImg = gdal.Open(sourceFileName, gdal.GA_ReadOnly)
array = sourceImg.GetRasterBand(1).ReadAsArray()*1.0
FDx, FDy = ndimage.convolve(array,mx),ndimage.convolve(array,my)
Image.fromarray(FDx).save(FDxFileName)
Image.fromarray(FDy).save(FDyFileName)

Strike = np.arctan2(-FDy,FDx)*radToDegFactor
Dip = np.arctan(np.sqrt(FDx*FDx+FDy*FDy))*radToDegFactor

Image.fromarray(Strike).save(StrikeFileName)
Image.fromarray(Dip).save(DipFileName)

output=mpimg.imread(cMapFileName)[Dip.astype(int),Strike.astype(int)]
Image.fromarray(output).save(outputFileName)
