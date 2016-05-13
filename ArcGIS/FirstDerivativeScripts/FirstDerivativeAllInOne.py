import numpy as np
import arcpy
import os
blockMemory = 0.0 # placeholder for memory leak debugging
#########                   ########
##        Class Definitions       ##
class recycledArrays: 
    def __init__(self, bH, bW): ##  ## bH = block height, bW=block width
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
        def makeAVA(out):
            self.arrAVA = np.zeros((3,bH-2*rV,bW-2*rH), np.uint8)
            self.arrLOW = np.zeros((bH-2*rV,bW-2*rH), np.int)
            if not out.Strike : prepArrays(["arrStrike"], "arrStrikeSm", np.uint16)
            if not out.Dip    : prepArrays(["arrDipX","arrDipY","arrDip"], "arrDipSm", np.uint8)
        makeLargeArrs(["arrXgrad","arrYgrad","arrTemp","arrOut", "arrTempMask", "arrMaskOut"])
        boolDeviation =  out.MSE <> "" or out.MAE <> "" or out.MSSD <> "" or out.MD <> ""
        if out.Xgrad        : self.arrXgradSm  = makeNewOut(np.float32)
        if out.Ygrad        : self.arrYgradSm  = makeNewOut(np.float32)
        if out.Strike       : prepArrays(["arrStrike"], "arrStrikeSm", np.uint16)
        if out.Dip          : prepArrays(["arrDipX","arrDipY","arrDip"], "arrDipSm", np.uint8)
        if boolDeviation    : prepArrays(["arrDevPt","arrDev","arrDevXgrad","arrDevYgrad","arrDevTemp"], "arrDevSm", np.float32)
        if out.AVA          : makeAVA(out)

class outRasterFile:   # datatype = "32_BIT_FLOAT", "16_BIT_UNSIGNED", "8_BIT_UNSIGNED", etc
    def __init__(self, address, seedName, datatype, castingType, noData):
        self._fileName = os.path.basename(address)
        self._directory = address[:len(address)-len(self._fileName)]
        self._castingType = castingType
        self._seedName = seedName
        self.tempRasterList = []
        self._noData = noData
        self._datatype = datatype
        self._numBands = 1
    def saveArr(self, arr, LLP):
        fileTemp = arcpy.CreateUniqueName(str(arcpy.env.workspace) + self._seedName)
        tempRaster = arcpy.NumPyArrayToRaster(arr, LLP, cW, cH)
        tempRaster.save(fileTemp)
        arcpy.AddMessage("Workspace file: " + fileTemp)
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
    def __init__(self, cW, cH, rH, rW, blockInStartX, blockInStartY, bW, bH):
        self._in = arcpy.Point(blockInStartX, blockInStartY)
        self._out = arcpy.Point(blockInStartX+cW*rH, blockInStartY+cH*rV)
        self._inWidth = bW
        self._inHeight = bH
        self._outWidth = bW - 2 * rH
        self._outHeight = bH - 2 * rV

class fileAddr: pass

class references: pass

class refIterator:
    def __init__(self, start): self.position = start
    def increment(self):self.position +=1
    def incrReturn(self):
        self.increment()
        return self.position

#########                      ########
##        Function Definitions       ##
def convolution(objArr, convoMatrix, bH, bW): #output stored in objArr.arrTemp
    objArr.arrTemp*=0
    objArr.arrTempMask*=0
    objArr.arrOut*=0
    objArr.arrMaskOut*=0
    for px in convoMatrix:
        fromX,toX, fromY,toY, mult = rH+px[0],bW-rH+px[0], rV+px[1],bH-rV+px[1], px[2] #0 = x coord, 1 = y coord, 2 = multiplier]
        np.multiply( objArr.inputNdArr[  fromY:toY, fromX:toX], mult,    objArr.arrTemp       )#"mult" is the grid spacing multiplier
        np.add(      objArr.arrTemp,                                   *[objArr.arrOut     ]*2)
        np.multiply( objArr.inputMaskArr[fromY:toY, fromX:toX], 1,    objArr.arrTempMask   )
        np.add(      objArr.arrTempMask,                               *[objArr.arrMaskOut ]*2)
    objArr.arrTemp = np.where(objArr.arrMaskOut>0, objArr.arrOut*np.nan, objArr.arrOut)
def prepareForConvo(inputArr, objArr):
    np.isnan(    inputArr, objArr.inputMaskArr)
    np.add(      0,    np.nan_to_num(inputArr), objArr.inputNdArr)
def neighbourhoodAverage(inputArr, bH, bW, objArr): #Warning: This function does not ignore null values.
    objArr.arrTemp*=0
    for x in range(-rH,rH+1):
        for y in range(-rV,rV+1):
            np.add(objArr.arrTemp, inputArr[rV+y:bH-rV+y, rH+x:bW-rH+x],objArr.arrTemp)
    np.divide(   objArr.arrTemp,       kernelArea,          objArr.arrTemp    )
def pointDeviation(y,x,objArr, bH, bW): #Warning: This function does not ignore null values.
    np.add(      objArr.arrSource[rV+y:bH-rV+y, rH+x:bW-rH+x], 0, objArr.arrDevPt    )      # = Z(x,y)
    np.multiply( objArr.arrXgrad,      x*cW,                      objArr.arrDevXgrad )      # =  A*x
    np.subtract( objArr.arrDevPt,      objArr.arrDevXgrad,        objArr.arrDevPt    )      # =  Z(x,y) - A*x
    np.multiply( objArr.arrYgrad,      y*cH,                      objArr.arrDevYgrad )      # =  B*y
    np.subtract( objArr.arrDevPt,      objArr.arrDevYgrad,        objArr.arrDevPt    )      # =  Z(x,y) - A*x - B*y
    np.subtract( objArr.arrDevPt,      objArr.arrTemp,            objArr.arrDevPt    )      # =  Z(x,y) - A*x - B*y - D
dummyFunction = lambda a,b,c,d,e,f: True #Dummy is only used by AVA
def writeXgradF(  ref, objArr): np.add(0, objArr.arrXgrad, objArr.arrXgradSm, dtype = np.float32)
def writeYgradF(  ref, objArr): np.add(0, objArr.arrYgrad, objArr.arrYgradSm, dtype = np.float32)
def writeStrikeF( ref, objArr):
    np.arctan2(  objArr.arrYgrad,      objArr.arrXgrad,    objArr.arrStrike)                # [1.] Find arctan from gradients.
    objArr.arrStrike    *= (180/math.pi)                                                    # [2.] Convert to degrees.
    np.rint( objArr.arrStrike, objArr.arrStrike)# try using round to nearest instead of adding 0.5
    np.add(      objArr.arrStrike,    180,      objArr.arrStrikeSm,     dtype = np.uint16)  # [3.] Type conversion. [UInt 0 to 359] # Changed from 179 to 180
def writeDipF(    ref, objArr):
    np.square(   objArr.arrXgrad,                           objArr.arrDipX    )             # [1.] Find Gradient squared in X.
    np.square(   objArr.arrYgrad,                           objArr.arrDipY    )             # [2.] -//------------------ in Y.
    np.add(      objArr.arrDipX,       objArr.arrDipY,      objArr.arrDip     )             # [3.] Add [1.] and [2.] together.
    np.sqrt(     objArr.arrDip,                             objArr.arrDip     )             # [4.] Take square root of [3.]
    np.arctan(   objArr.arrDip,                             objArr.arrDip     )             # [5.] Take arctan of [4.]
    objArr.arrDip       *= (180/math.pi)                                                    #       convert to degrees.
    np.rint( objArr.arrDip, objArr.arrDip) # try using round to nearest instead of adding 0.5
    np.multiply( objArr.arrDip,        1,       objArr.arrDipSm,         dtype = np.uint8)  # [6.] Type conversion.
def dFuncReturn(procFunc, refRas):
    def tempFunc(ref, objArr, bH, bW, dataPoints, coord):
        procFunc(ref, objArr)
        objRasters[getattr(ref,refRas)].saveArr( getattr(objArr,str("arr"+refRas+"Sm")), coord._out )
    return tempFunc
def ErrorFuncMSE(  objArr ): #These functions calculate deviation in variety of ways:       # Mean Squared Error
    np.square(                                   *[objArr.arrDevPt  ]*2 )
def ErrorFuncMAE(  objArr ):                                                                # Mean Absolute Error
    np.absolute(                                 *[objArr.arrDevPt  ]*2 )
def ErrorFuncMSSD( objArr ):                                                                # Mean Signed Squared Deviation
    np.absolute(            objArr.arrDevPt,       objArr.arrDevTemp    )
    np.multiply(            objArr.arrDevTemp,   *[objArr.arrDevPt  ]*2 )
def ErrorFuncMD(   objArr ):                                                                # Maximum Absolute Deviation
    np.absolute(                                 *[objArr.arrDevPt  ]*2 )
    np.maximum(             objArr.arrDevPt,     *[objArr.arrDevTemp]*2 )
    np.add(0,               objArr.arrDevTemp,        objArr.arrDevPt   )
def computeErr(objArr, bH, bW, dataPoints, ErrorFunc): #Returns an Error function out of the ones above
    objArr.arrDevTemp*=0
    objArr.arrOut*=0
    neighbourhoodAverage(objArr.arrSource, bH, bW, objArr )        # => Computes average into objArr.arrTemp
    for dataPoint in dataPoints:
        pointDeviation(dataPoint[0], dataPoint[1], objArr, bH, bW) #compute point deviation => modifies objArr.arrPointDev
        ErrorFunc(objArr)
        np.add(objArr.arrDevPt, objArr.arrOut, objArr.arrOut)
    np.divide(objArr.arrOut, float(len(dataPoints)), objArr.arrDevSm, dtype = np.float32)#Calculates Error, casts into 32 float
def eFuncReturn(ProcessingFunc, refRas ):
    def tempFunction(ref, objArr, bH, bW, dataPoints, coord):
        computeErr(objArr, bH, bW, dataPoints, ProcessingFunc)
        objRasters[getattr(ref, refRas)].saveArr(   objArr.arrDevSm, coord._out)
    return tempFunction
def prepareOutput(propType,   oPut, outPlace, funCall,   toRun,ref,refIndexing,objRasters):
    setattr(ref, oPut, refIndexing.incrReturn()) # incrementing reference counter
    objRasters.append(outRasterFile(outPlace, "\\out" + oPut, *propType))
    toRun.append(funCall)
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
    
def classicAVA(ref, objArr, bH, bW, dataPoints, coord):
    if not out.Strike: writeStrikeF( "a", objArr)#arrays will be created but not saved, since we're calling FuncReturn
    if not out.Dip:    writeDipF(    "a", objArr)#skipping instead to the inner function for making arrays
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

def processOneBlock(coord, objArr, inSourceRaster, elevationFactor, objRasters, ref, toRun):
    bH, bW = coord._inHeight, coord._inWidth
    np.multiply(elevationFactor, arcpy.RasterToNumPyArray(inSourceRaster, coord._in, bW, bH), objArr.arrSource)
    prepareForConvo(objArr.arrSource, objArr)     # this fills arrays used for convolution, prepare null values as zeros
    convolution(objArr, [[x,y,x] for x in range(-rH,rH+1) for y in range(-rV,rV+1)], bH, bW) #Convolution gradient X:
    np.multiply(objArr.arrTemp, normalizationX, objArr.arrXgrad) #multiply by normalization factor computed at input collection
    convolution(objArr, [[x,y,y] for x in range(-rH,rH+1) for y in range(-rV,rV+1)], bH, bW) # -//-                Y:
    np.multiply(objArr.arrTemp, normalizationY, objArr.arrYgrad)
    dataPoints = [[x,y] for x in range(-rH,rH+1) for y in range(-rV,rV+1)] #list dataPoints contains kernel matrix for deviation
    for funct in toRun: funct(ref, objArr, bH, bW, dataPoints, coord) #see function prepareOutput()

def processTileset(blockCounter, processOneBlockParams, blockMemory,PSs, blockHeight, blockWidth):
    objArr = recycledArrays(blockHeight, blockWidth)                    # different block dimensions on edge blocks
    coordLLPPs = [coordsOfBlock(cW, cH, rH, rV, *PS) for PS in PSs]
    for coordLLPP in coordLLPPs:
        stepMemory = 0 #debug placeholder
        processOneBlock(coordLLPP, objArr, *processOneBlockParams) #defined outside the processTileset,
        blockCounter[0]+=1                                              # as processTileset is called multiple times
        arcpy.AddMessage("Block number " + str(blockCounter[0]) + " processed, out of " + str(totalBlocks))
    arcpy.AddMessage("All blocks computed")
    del objArr

#########                  ########
##        Input Collection       ##
outputIter      =         refIterator(-1)                                    # iterator object to remember our input
grabInput       = lambda: arcpy.GetParameterAsText(outputIter.incrReturn())  # anonymus function to collect input
inSourceRaster  =         arcpy.Raster(grabInput())   #input raster layer
elevationFactor =         float(       grabInput())   #vertical exaggeration ## PROPER NAME IS z-factor
rV              =         int(         grabInput())   #radius North-South
rH              =         int(         grabInput())   #radius East-West
propFloat       =  ["32_BIT_FLOAT",    np.float32, str(-2**20) ]
propUInt        =  ["16_BIT_UNSIGNED", np.uint16,  str(400)    ]
propUByte       =  ["8_BIT_UNSIGNED",  np.uint8,   str(100)    ]
outputs = [ ["Xgrad", propFloat, dFuncReturn(writeXgradF,   "Xgrad")],["Ygrad",  propFloat, dFuncReturn(writeYgradF,  "Ygrad") ],
            ["Dip",   propUByte, dFuncReturn(writeDipF,     "Dip")  ],["Strike", propUInt,  dFuncReturn(writeStrikeF, "Strike")],
            ["MSE",   propFloat, eFuncReturn(ErrorFuncMSE,  "MSE")  ],["MAE",    propFloat, eFuncReturn(ErrorFuncMAE, "MAE")   ],
            ["MSSD",  propFloat, eFuncReturn(ErrorFuncMSSD, "MSSD") ],["MD",     propFloat, eFuncReturn(ErrorFuncMD,  "MD")    ],
            ["AVA",   propUByte, classicAVA                         ]] # from 4 to 12
out = fileAddr()
for output in outputs:
    setattr(out, output[0], grabInput())
inColourMapRef    = grabInput()
if inColourMapRef <> '': #if user wants to use a custom colour map
    inColourMap     =         arcpy.Raster(inColourMapRef)
    colourMap = arcpy.RasterToNumPyArray(inColourMap)
else:
    colourMap = colourMapClassic() 

#########                             ########
##        Global Variable Definitions       ##
#                                            #
inputCoordinateSystem = inSourceRaster.spatialReference
cW,cH = inSourceRaster.meanCellWidth, inSourceRaster.meanCellHeight      # cell size
kernelWidth, kernelHeight = rH*2+1, rV*2+1                               # kernel size
kernelArea = float(kernelHeight*kernelWidth)                             # kernel area
irNDV = inSourceRaster.noDataValue                                       # input NODATA value - do i even need this??
###       Computing block sizes
blockWidth, blockHeight = 2**11, 2**11                                   #  maximum block size
irWidth, irHeight = inSourceRaster.width, inSourceRaster.height          #  raster size
steps = irWidth / (blockWidth-rH-1)                                     # # find number of steps
rows = irHeight / (blockHeight-rV -1)                                   # #                rows
totalBlocks = (steps+1)*(rows+1)                                         #  Total number of blocks
edgeWidth  = irWidth  % (blockWidth -rH -1)                             # # Edge width
edgeHeight = irHeight % (blockHeight-rV -1)                             # #      height
startX,startY = inSourceRaster.extent.XMin,inSourceRaster.extent.YMin    # Lower Left corner of source image
shiftX,shiftY = cW*(blockWidth-2*rH-1),cH*(blockHeight-2*rV-1)               # Displacement for regular tiles ## added mult by 2
###     Computing normalization factors
normalizationX = 3/(cW*rH*(rH+1)*kernelArea)
normalizationY = 3/(cH*rV*(rV+1)*kernelArea)
###  Environment setting, probably unnecessary, but just in case
arcpy.env.overwriteOutput = True              
#  ########                     ######  #
# #         PROGRAM STARTS HERE       # #
#                                       #
PSs = [[shiftX*step +startX, shiftY*row +startY, blockWidth, blockHeight] for step in range(steps) for row in range(rows)] # Regular tiles
REs = [[shiftX*steps+startX, shiftY*row +startY, edgeWidth,  blockHeight] for row in range(rows)]                          # Right edge tiles
TEs = [[shiftX*step +startX, shiftY*rows+startY, blockWidth, edgeHeight ] for step in range(steps)]                        # Top edge tiles
OB =  [[shiftX*steps+startX, shiftY*rows+startY, edgeWidth,  edgeHeight]]                                                  # Corner block
Tiles = [[PSs,blockHeight,blockWidth], [REs,blockHeight,edgeWidth], [TEs,edgeHeight,blockWidth], [OB,edgeHeight,edgeWidth]]# Combine into list
blockCounter  = [0]                                                       # Start counting blocks
ref, refIndexing = references(), refIterator(-1)                        # Set up referencing for rasters
objRasters, toRun    =  [],[]                                           # ## # list: rasters, function calls
prepOutArgs   = [toRun, ref, refIndexing, objRasters]
for output in outputs: #         prepareOutput(dataType, name,        address,                 function,  *args)
    if getattr(out, output[0]) : prepareOutput(output[1], output[0],  getattr(out, output[0]), output[2], *prepOutArgs )
for TileSet  in Tiles: processTileset(blockCounter, [inSourceRaster, elevationFactor, objRasters, ref, toRun], blockMemory,*TileSet)
for oRas     in objRasters: oRas.saveMosaic()
