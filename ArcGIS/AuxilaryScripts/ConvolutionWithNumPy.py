blockMemory = 0 ### Debug placeholder
import numpy as np
import arcpy
import os
#########                   ########
##        Class Definitions       ##
class recycledArrays: 
    def __init__(self, bH, bW): ##  ## bH = block height, bW=block width
        arcpy.AddMessage("recycledArrays bH=" + str(bH) + " , bW=" + str(bW)) #DEBUG
        self.arrSource    =  np.zeros((bH,     bW),      np.float64)
        self.inputMaskArr =  np.zeros((bH,     bW),      np.float64) # Holds 1.'s for NAN, else 0.
        self.inputNdArr   =  np.zeros((bH,     bW),      np.float64) # Input arr, 0.'s for NAN, else data
        def makeNewOut(dataType):
            return np.zeros((bH-2*rV,bW-2*rH), dataType)
        def makeLargeArrs(listArr):
            for Ras in listArr: setattr(self, Ras, makeNewOut(np.float64))
        def prepArrays(listArr, smName, smType):
            makeLargeArrs(listArr)
            setattr(self, smName, makeNewOut(smType))
        makeLargeArrs(["arrConvo","arrTemp","arrOut", "arrTempMask", "arrMaskOut"])
        self.arrConvoSm  = makeNewOut(np.float32)

class outRasterFile:   # datatype = "32_BIT_FLOAT", "16_BIT_UNSIGNED", "8_BIT_UNSIGNED", etc
    def __init__(self, address, seedName, datatype, castingType, noData):
        self._address = str(address)
        self._fileName = os.path.basename(address)
        arcpy.AddMessage("Raster address: " + address)  ### DEBUG
        self._directory = address[:len(address)-len(self._fileName)]
        self._castingType = castingType
        self._seedName = seedName
        self.tempRasterList = []
        self._noData = noData
        self._datatype = datatype
        self._numBands = 1
    def saveArr(self, arr, LLP):
        fileTemp = arcpy.CreateUniqueName(str(arcpy.env.workspace) + self._seedName)
        arcpy.AddMessage("fileTemp: " + fileTemp)
        arr[np.isnan(arr)] = float(self._noData)
        tempRaster = arcpy.NumPyArrayToRaster(arr, LLP, cW, cH, 0)
        tempRaster.save(fileTemp)
        del(arr)
        del(tempRaster)
        self.tempRasterList.append(fileTemp[:])
    def saveMosaic(self):
        arcpy.AddMessage("Raster directory: " + self._directory)
        arcpy.AddMessage("Raster filename: " + self._fileName)
        arcpy.AddMessage("tempRasterList: " + ';'.join(self.tempRasterList))
        arcpy.MosaicToNewRaster_management(';'.join(self.tempRasterList),
                                           self._directory, self._fileName,
                                           coordinate_system_for_the_raster = inputCoordinateSystem,
                                           pixel_type = self._datatype,
                                           number_of_bands = self._numBands)
        arcpy.AddMessage("1 " + self._noData +".00")
        arcpy.SetRasterProperties_management(self._address, nodata="1 " + self._noData +".00")
        arcpy.CalculateStatistics_management(self._address)
        for tempRasterFile in self.tempRasterList:
            if arcpy.Exists(tempRasterFile):
                arcpy.Delete_management(tempRasterFile)

class coordsOfBlock:
    def __init__(self, cW, cH, rH, rW, blockInStartX, blockInStartY, bW, bH):
        self._in = arcpy.Point(blockInStartX, blockInStartY)
        self._out = arcpy.Point(blockInStartX+cW*rH, blockInStartY+cH*rV)
        self._inWidth = bW
        self._inHeight = bH
        self._outWidth = bW - 2 * rH
        self._outHeight = bH - 2 * rV

class voidClass: pass

class refIterator:
    def __init__(self, start): self.position = start
    def increment(self):self.position +=1
    def incrReturn(self):
        self.increment()
        return self.position

#########                      ########
##        Function Definitions       ##
#                                     #
def convolution(objArr, convoMatrix, bH, bW): #output stored in objArr.arrTemp
    objArr.arrTemp*=0
    objArr.arrTempMask*=0
    objArr.arrOut*=0
    objArr.arrMaskOut*=0
    for px in convoMatrix: 
        fromX,toX, fromY,toY, mult = rH+px[0],bW-rH+px[0], rV+px[1],bH-rV+px[1], px[2] #0 = x coord, 1 = y coord, 2 = multiplier]
        np.multiply( objArr.inputNdArr[  fromY:toY, fromX:toX], mult,    objArr.arrTemp       )
        np.add(      objArr.arrTemp,                                   *[objArr.arrOut     ]*2)
        np.multiply( objArr.inputMaskArr[fromY:toY, fromX:toX], 1,    objArr.arrTempMask   )
        np.add(      objArr.arrTempMask,                               *[objArr.arrMaskOut ]*2)
    objArr.arrTemp = np.where(objArr.arrMaskOut>0, objArr.arrOut*np.nan, objArr.arrOut)

def prepareForConvo(inputArr, objArr):
    np.isnan(    inputArr, objArr.inputMaskArr) #put all nodata values into
    np.add(      0,    np.nan_to_num(inputArr), objArr.inputNdArr)

def dfWriteConvo(ref, objArr, bH, bW, dataPoints, coord):
    np.add(0, objArr.arrConvo, objArr.arrConvoSm, dtype = np.float32)
    objRasters[ref.Convo].saveArr( objArr.arrConvoSm, coord._out )

def processOneBlock(coord, objArr, inSourceRaster, elevationFactor, objRasters, ref, toRun):
    bH, bW = coord._inHeight, coord._inWidth
    np.multiply(1, arcpy.RasterToNumPyArray(inSourceRaster, coord._in, bW, bH, irNDV), objArr.arrSource)
    np.multiply(elevationFactor, np.where(objArr.arrSource==irNDV, objArr.arrSource*np.nan, objArr.arrSource), objArr.arrSource)
    arcpy.AddMessage(str(irNDV))
    prepareForConvo(objArr.arrSource, objArr)     # this fills arrays used for convolution, prepare null values as zeros
    convolution(objArr, mtrx, bH, bW)
    np.add(objArr.arrTemp, 0, objArr.arrConvo)
    dataPoints = [[x,y] for x in range(-rH,rH+1) for y in range(-rV,rV+1)] #list dataPoints contains kernel matrix for deviation
    for funct in toRun: funct(ref, objArr, bH, bW, dataPoints, coord) #see function prepareOutput()

def processTileset(blockCounter, processOneBlockParams, blockMemory,PSs, blockHeight, blockWidth):
    if (blockHeight < 2*rV ) | (blockWidth<2*rH):
        return
    objArr = recycledArrays(blockHeight, blockWidth)                    # different block dimensions on edge blocks
    coordLLPPs = [coordsOfBlock(cW, cH, rH, rV, *PS) for PS in PSs]
    for coordLLPP in coordLLPPs:
        stepMemory = 0
        processOneBlock(coordLLPP, objArr, *processOneBlockParams) #defined outside the processTileset,
        blockCounter[0]+=1                                              # as processTileset is called multiple times
        arcpy.AddMessage("Block number " + str(blockCounter[0]) + " processed, out of " + str(totalBlocks))
    arcpy.AddMessage("All blocks computed")                                                                       ###DEBUG
    del objArr

def incr(self, num):
    self.counter+=num
    return self.counter

#########                  ########
##        Input Collection       ##
outputIter      =         refIterator(-1)                                    # iterator object to remember our input
grabInput       = lambda: arcpy.GetParameterAsText(outputIter.incrReturn())  # anonymus function to collect input
inSourceRaster  =         arcpy.Raster(grabInput())   #input raster layer
elevationFactor =         float(       grabInput())   #vertical exaggeration ## PROPER NAME IS z-factor
rV              =         int(         grabInput())   #radius North-South
rH              =         int(         grabInput())   #radius East-West
convoKernelStr  =         str(         grabInput())   #obtain convolution matrix as string

propFloat       =  ["32_BIT_FLOAT",    np.float32, str(-2**20) ]
out = voidClass()
out.Convo = grabInput()
##        Global Variable Definitions       ##
inputCoordinateSystem = inSourceRaster.spatialReference
cW,cH = inSourceRaster.meanCellWidth, inSourceRaster.meanCellHeight      # cell size
kernelWidth, kernelHeight = rH*2+1, rV*2+1                               # kernel size
kernelArea = float(kernelHeight*kernelWidth)                             # kernel area
irNDV = inSourceRaster.noDataValue                                       # input NODATA value
convoKernelStrSeq = convoKernelStr.split(',')                   #Converting Kernel from string to a matrix of numbers
convoKernelFloatSeq = [float(u) for u in convoKernelStrSeq]
itr = voidClass()
itr.counter = -1
mtrx = [[-rH+incr(itr,1)%kernelWidth, -rV+itr.counter/kernelWidth, u] for u in convoKernelFloatSeq]

###       Computing block sizes
blockWidth, blockHeight = 2**11, 2**11                                   #  maximum block size
irWidth, irHeight = inSourceRaster.width, inSourceRaster.height          #  raster size
steps = irWidth / (blockWidth-2*rH-1)                                     # # find number of steps
rows = irHeight / (blockHeight-2*rV -1)                                   # #                rows
totalBlocks = (steps+1)*(rows+1)                                         #  Total number of blocks
edgeWidth  = irWidth  % (blockWidth -2*rH -1)                             # # Edge width
edgeHeight = irHeight % (blockHeight-2*rV -1)                             # #      height
startX,startY = inSourceRaster.extent.XMin,inSourceRaster.extent.YMin    # Lower Left corner of source image
shiftX,shiftY = cW*(blockWidth-2*rH),cH*(blockHeight-2*rV)               # Displacement for regular tiles
# #         PROGRAM STARTS HERE       # #
PSs = [[shiftX*step +startX, shiftY*row +startY, blockWidth, blockHeight] for step in range(steps) for row in range(rows)] # Regular tiles
REs = [[shiftX*steps+startX, shiftY*row +startY, edgeWidth,  blockHeight] for row in range(rows)]                          # Right edge tiles
TEs = [[shiftX*step +startX, shiftY*rows+startY, blockWidth, edgeHeight ] for step in range(steps)]                        # Top edge tiles
OB =  [[shiftX*steps+startX, shiftY*rows+startY, edgeWidth,  edgeHeight]]                                                  # Corner block
Tiles = [[PSs,blockHeight,blockWidth], [REs,blockHeight,edgeWidth], [TEs,edgeHeight,blockWidth], [OB,edgeHeight,edgeWidth]]# Combine into list
blockCounter  = [0]                                                       # Start counting blocks
ref, refIndexing = voidClass(), refIterator(-1)                        # Set up referencing for rasters
objRasters, toRun    =  [],[]                                           # ## # initialize lists for: rasters, function calls
ref.Convo = refIndexing.incrReturn()
objRasters.append(outRasterFile(out.Convo, "\\outConvo", *propFloat))
toRun.append(dfWriteConvo)
for TileSet  in Tiles: processTileset(blockCounter, [inSourceRaster, elevationFactor, objRasters, ref, toRun], blockMemory,*TileSet)
for oRas     in objRasters: oRas.saveMosaic()
