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
            return np.zeros((bH,bW), dataType)
        def makeLargeArrs(listArr):
            for Ras in listArr: setattr(self, Ras, makeNewOut(np.float64))
        def prepArrays(listArr, smName, smType):
            makeLargeArrs(listArr)
            setattr(self, smName, makeNewOut(smType))
        def makeAVA(out):
            self.arrAVA = np.zeros((3,bH,bW), np.uint8)
            self.arrLOW = np.zeros((bH,bW), np.int)
            prepArrays(["arrStrike"], "arrStrikeSm", np.uint16)
            prepArrays(["arrDipX","arrDipY","arrDip"], "arrDipSm", np.uint8)
        makeAVA(out)

class outRasterFile:   # datatype = "32_BIT_FLOAT", "16_BIT_UNSIGNED", "8_BIT_UNSIGNED", etc
    def __init__(self, address, seedName, datatype, castingType, noData):
        self._fileName = os.path.basename(address)
        arcpy.AddMessage("Raster address: " + address)
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
        tempRaster = arcpy.NumPyArrayToRaster(arr, LLP, cW, cH)
        tempRaster.save(fileTemp)
        arcpy.AddMessage("fileTemp: " + fileTemp)
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
        for tempRasterFile in self.tempRasterList:
            if arcpy.Exists(tempRasterFile):
                arcpy.Delete_management(tempRasterFile)

class coordsOfBlock:
    def __init__(self, cW, cH, blockInStartX, blockInStartY, bW, bH):
        self._in = arcpy.Point(blockInStartX, blockInStartY)
        self._out = arcpy.Point(blockInStartX, blockInStartY)
        self._inWidth   = bW  #AVA script and Convolution script are built from
        self._inHeight  = bH  #using legacy code of All-In-One First Derivative tool
        self._outWidth  = bW  #AVA edges do not shrink, so "in" and "out" dimensions 
        self._outHeight = bH  #are the same

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
def classAVAmap(dip, strike):# 0-360, 0-90
    high = dip / 90.0
    mid = strike / 60.0
    huePrime = mid
    mid %= 2.0
    mid  -= 1
    mid = abs(mid)
    mid *= high
    low = high / 2.0
    low = 0.5 - low
    mid +=  low
    high += low
    lowC = int(low*254)
    midC = int(mid*254)
    hghC = int(high*254)
    if huePrime<1:
        rgb = [lowC, midC, hghC]
    elif huePrime<2:
        rgb = [midC, lowC, hghC]
    elif huePrime<3:
        rgb = [hghC, lowC, midC]
    elif huePrime<4:
        rgb = [hghC, midC, lowC]
    elif huePrime<5:
        rgb = [midC, hghC, lowC]
    elif huePrime<6:
        rgb = [lowC, hghC, midC]
    else:
        rgb = [lowC, midC, hghC]
    return rgb

def colourMapClassic():
    ClassicCM = np.zeros([3,90,360],np.uint8)
    for dip in range(90):
        for strike in range(360):
            rgb = classAVAmap(dip, strike)
            for channel in range(3):
                ClassicCM[channel][dip][strike] = rgb[channel]
    return ClassicCM
    
def classicAVA(ref, objArr, bH, bW, dataPoints, coord, colourMap):
    refRas = "AVA"
    objRasters[getattr(ref,refRas)]._numBands = 3
    objArr.arrLOW *= 0
    redMap = np.ravel(colourMap[0])
    np.multiply(objArr.arrDipSm, 360, objArr.arrLOW)
    np.add(objArr.arrStrikeSm, *[objArr.arrLOW]*2)
    np.take(redMap, objArr.arrLOW, out = objArr.arrAVA[0], mode = 'wrap')
    grnMap = np.ravel(colourMap[1])
    np.take(grnMap, objArr.arrLOW, out = objArr.arrAVA[1], mode = 'wrap')
    bluMap = np.ravel(colourMap[2])
    np.take(bluMap, objArr.arrLOW, out = objArr.arrAVA[2], mode = 'wrap')
    objRasters[getattr(ref,refRas)].saveArr( getattr(objArr,str("arr"+refRas)), coord._out )

def processOneBlock(coord, objArr, inStrikeRaster, inDipRaster, objRasters, ref, toRun, colourMap):
    bH, bW = coord._inHeight, coord._inWidth
    np.multiply(1, arcpy.RasterToNumPyArray(inStrikeRaster, coord._in, bW, bH), objArr.arrStrikeSm)
    np.multiply(1, arcpy.RasterToNumPyArray(inDipRaster, coord._in, bW, bH), objArr.arrDipSm)
    dataPoints = 0
    classicAVA(ref, objArr, bH, bW, dataPoints, coord, colourMap)

def processTileset(blockCounter, processOneBlockParams, blockMemory,PSs, blockHeight, blockWidth):
    objArr = recycledArrays(blockHeight, blockWidth)                    # different block dimensions on edge blocks
    coordLLPPs = [coordsOfBlock(cW, cH, *PS) for PS in PSs]
    for coordLLPP in coordLLPPs:
        stepMemory = 0
        processOneBlock(coordLLPP, objArr, *processOneBlockParams) #defined outside the processTileset,
        blockCounter[0]+=1                                              # as processTileset is called multiple times
        arcpy.AddMessage("Block number " + str(blockCounter[0]) + " processed, out of " + str(totalBlocks))
    arcpy.AddMessage("All blocks computed")
    del objArr

def incr(self, num):
    self.counter+=num
    return self.counter

#########                  ########
##        Input Collection       ##
outputIter      =         refIterator(-1)                                    # iterator object to remember our input
grabInput       = lambda: arcpy.GetParameterAsText(outputIter.incrReturn())  # anonymus function to collect input
inStrikeRaster    =         arcpy.Raster(grabInput())
inDipRaster       =         arcpy.Raster(grabInput())
inColourMapRef  = grabInput()
if inColourMapRef <> '': #THIS GIVES AN OPTION TO USE DIFFERENT COLOUR MAP
    inColourMap     =         arcpy.Raster(inColourMapRef)
    colourMap = arcpy.RasterToNumPyArray(inColourMap)
else:
    colourMap = colourMapClassic()
    
propUByte       =  ["8_BIT_UNSIGNED",  np.uint8,   str(100)    ]
out = voidClass()
out.AVA = grabInput()

##        Global Variable Definitions       ##
inputCoordinateSystem = inStrikeRaster.spatialReference
cW,cH = inStrikeRaster.meanCellWidth, inStrikeRaster.meanCellHeight      # cell size
irNDV = inStrikeRaster.noDataValue

###       Computing block sizes
blockWidth, blockHeight = 2**11, 2**11                                   #  maximum block size
irWidth, irHeight = inStrikeRaster.width, inStrikeRaster.height          #  raster size
steps = irWidth / (blockWidth-1)                                     # # find number of steps
rows = irHeight / (blockHeight -1)                                   # #                rows
totalBlocks = (steps+1)*(rows+1)                                         #  Total number of blocks
edgeWidth  = irWidth  % (blockWidth -1)                             # # Edge width
edgeHeight = irHeight % (blockHeight-1)                             # #      height
startX,startY = inStrikeRaster.extent.XMin,inStrikeRaster.extent.YMin    # Lower Left corner of source image
shiftX,shiftY = cW*(blockWidth-1),cH*(blockHeight-1)               # Displacement for regular tiles

# #         PROGRAM STARTS HERE       # #
PSs = [[shiftX*step +startX, shiftY*row +startY, blockWidth, blockHeight] for step in range(steps) for row in range(rows)] # Regular tiles
REs = [[shiftX*steps+startX, shiftY*row +startY, edgeWidth,  blockHeight] for row in range(rows)]                          # Right edge tiles
TEs = [[shiftX*step +startX, shiftY*rows+startY, blockWidth, edgeHeight ] for step in range(steps)]                        # Top edge tiles
OB =  [[shiftX*steps+startX, shiftY*rows+startY, edgeWidth,  edgeHeight]]                                                  # Corner block
Tiles = [[PSs,blockHeight,blockWidth], [REs,blockHeight,edgeWidth], [TEs,edgeHeight,blockWidth], [OB,edgeHeight,edgeWidth]]# Combine into list
blockCounter  = [0]                                                       # Start counting blocks

ref, refIndexing = voidClass(), refIterator(-1)                        # Set up referencing for rasters
objRasters, toRun    =  [],[]                                           # ## # list: rasters, function calls

ref.AVA = refIndexing.incrReturn()
objRasters.append(outRasterFile(out.AVA, "\\outAVA", *propUByte))
toRun.append(classicAVA)

for TileSet  in Tiles: processTileset(blockCounter, [inStrikeRaster, inDipRaster, objRasters, ref, toRun, colourMap], blockMemory,*TileSet)
for oRas     in objRasters: oRas.saveMosaic()
