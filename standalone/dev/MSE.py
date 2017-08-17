import sys

sourceFileName = sys.argv[1] # 'Mikhail_test.tif'
if len(sys.argv)>=4: kx,ky,spc = [int(sys.argv[2]),int(sys.argv[3]),float(sys.argv[4])]
else:                kx,ky,spc = [3,3,50]
print(sys.argv[1])
print(spc)
cMapFileName   = 'CMap.tif'
FDxFileName    = 'FDx.tif';    FDyFileName    = 'FDy.tif'
StrikeFileName = 'Strike.tif'; DipFileName    = 'Dip.tif'
outputFileName = 'output.tif'; outputEnhFileName='outputEnhanced.tif'

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

N = kx*ky

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

DipClip=Dip/15
DipEnh=DipClip*(DipClip<1)+(DipClip>1)
DipEnh*=89
outputEnh=mpimg.imread(cMapFileName)[DipEnh.astype(int),Strike.astype(int)]
Image.fromarray(outputEnh).save(outputEnhFileName)

######## COMPUTE SSE:

SSEgX=(FDx*FDx*khx*(khx+1)*N*N)/3.
SSEgY=(FDy*FDy*khy*(khy+1)*N*N)/3.

KerSum=np.ones([kx,ky])
array=array.astype(np.float64)
arraySq=array*array
arraySqSum=ndimage.convolve(arraySq,KerSum)
arraySum=ndimage.convolve(array,KerSum)
arraySumSq=arraySum*arraySum
#Image.fromarray(arraySqSum).save("arraySqSum.tif")
#Image.fromarray(arraySumSq/N).save("arraySumSq.tif")
#adjCoeff=9
SSE=arraySqSum-arraySumSq/N - (SSEgX + SSEgY)
Image.fromarray(SSE/N).save("cMSE.tif")

#Image.fromarray(arraySqSum-arraySumSq/N).save("SEEsums.tif")
#SEEmatch=SSEgX + SSEgY
#Image.fromarray(SEEmatch).save("SEEmatch.tif")
