import arcpy       
import math

def getCircKernYY(kernelRadius, spacing):
    GradNSraw     = generateNorthSouthWeigths(                 kernelRadius)
    return fixCircCurvKern(GradNSraw, kernelRadius, spacing)

def getCircKernXX(kernelRadius, spacing):
    GradEWraw     = generateEastWestWeigths(                   kernelRadius)
    return fixCircCurvKern(GradEWraw, kernelRadius, spacing)

def getCircKernXY(kernelRadius, spacing):
    GradEWraw     = generateEastWestWeigths(                   kernelRadius)
    GradNSraw     = generateNorthSouthWeigths(                 kernelRadius)
    GeoTorsion    = Hadamard(GradEWraw,       GradNSraw,       kernelRadius)
    CircleWeights = generateCircleWeigths(                     kernelRadius)
    GeoTorCirc    = Hadamard(GeoTorsion,      CircleWeights,   kernelRadius)
    GeoTorSq      = Hadamard(GeoTorsion,      GeoTorsion,      kernelRadius)
    GeoTorSqCirc  = Hadamard(GeoTorSq,  CircleWeights,    kernelRadius)
    torNorm       = sumArray(GeoTorSqCirc,                     kernelRadius)
    curNorm       = torNorm*spacing*spacing
    return normalizeKernel(GeoTorCirc, 1/curNorm, kernelRadius)

######### This function will normalize kernels in XX and YY
def fixCircCurvKern(GradRaw, kernelRadius, spacing):
    GradEWraw     = generateEastWestWeigths(                   kernelRadius)
    GradNSraw     = generateNorthSouthWeigths(                 kernelRadius)
    GeoTorsion    = Hadamard(GradEWraw,       GradNSraw,       kernelRadius)
    CircleWeights = generateCircleWeigths(                     kernelRadius)
    GeoTorCirc    = Hadamard(GeoTorsion,      CircleWeights,   kernelRadius)
    GeoTorSq      = Hadamard(GeoTorsion,      GeoTorsion,      kernelRadius)
    GeoTorSqCirc  = Hadamard(GeoTorSq,       CircleWeights,    kernelRadius)
    torNorm       = sumArray(GeoTorSqCirc,                     kernelRadius) #This gets torsion normalization XXYY
    #
    GradSq        = Hadamard(GradRaw,         GradRaw,       kernelRadius)  #find squared gradient
    GradQu        = Hadamard(GradSq,          GradSq,        kernelRadius)  #find quadrated gradient
    #
    GradSqCirc    = Hadamard(GradSq,          CircleWeights, kernelRadius)  #find circle of squared gradient
    GradQuCirc    = Hadamard(GradQu,          CircleWeights, kernelRadius)  #find circle of quadrated gradient
    #
    areaCirc      = sumArray(CircleWeights,                  kernelRadius)  #find area of the circle
    sumGradSqCirc = sumArray(GradSqCirc,                     kernelRadius)
    gradSqCircAvg = sumGradSqCirc/areaCirc                  #find average value
    #
    OutKernelRaw  = addToKernel(GradSq,      -gradSqCircAvg,     kernelRadius) #subtract
    OutKernelCirc = Hadamard(OutKernelRaw,    CircleWeights, kernelRadius)
    gradNormQu    = sumArray(GradQuCirc,                     kernelRadius)
    curNorm       = (gradNormQu - torNorm)*spacing*spacing # note that nominator is 2 orders below denominator
    return normalizeKernel(OutKernelCirc, 1/curNorm, kernelRadius)

def getCircKernEW(kernelRadius, spacing):
    GradEWraw     = generateEastWestWeigths(                   kernelRadius)
    CircleWeights = generateCircleWeigths(                     kernelRadius)
    GradEWcirc    = Hadamard(GradEWraw,         CircleWeights, kernelRadius)
    GradEWsq      = Hadamard(GradEWraw,         GradEWraw,     kernelRadius)
    GradEWsqCirc  = Hadamard(GradEWsq,          CircleWeights, kernelRadius)
    gradEWnorm    = sumArray(GradEWsqCirc,                     kernelRadius)
    return normalizeKernel(GradEWcirc, 1/(spacing*gradEWnorm), kernelRadius)

def getCircKernNS(kernelRadius, spacing):
    GradNSraw     = generateNorthSouthWeigths(                 kernelRadius)
    CircleWeights = generateCircleWeigths(                     kernelRadius)
    GradNScirc    = Hadamard(GradNSraw,         CircleWeights, kernelRadius)
    GradNSsq      = Hadamard(GradNSraw,         GradNSraw,     kernelRadius)
    GradNSsqCirc  = Hadamard(GradNSsq,          CircleWeights, kernelRadius)
    gradNSnorm    = sumArray(GradNSsqCirc,                     kernelRadius)
    return normalizeKernel(GradNScirc, 1/(spacing*gradNSnorm), kernelRadius)

def evaluateIntegral(r,j):
   return 0.5*((r**2)*math.asin(j/r)+j*math.sqrt(r**2-(j)**2))

def calculateOverlay(j,k,r):
   Pj         = j-0.5
   Pk         = k-0.5
   Tj         = j+0.5
   Tk         = k+0.5
   if (Pj**2+Pk**2) >= (r**2): return 0.0
   if (Tj**2+Tk**2) <= (r**2): return 1.0
   Pjp         = math.sqrt(r**2-(Tk)**2)
   Pkp         = math.sqrt(r**2-(Tj)**2)
   Gj           = (Pjp, Pj)[Pj >= Pjp]
   Gk           = math.sqrt(r**2-Gj**2)
   Qk           = (Pkp, Pk)[Pk >= Pkp]
   Qj           = math.sqrt(r**2-Qk**2)
   Ipj           = Gj-Pj
   Ipk           = Qk-Pk
   areaLeft       = Ipj
   areaBottom     = Ipk
   areaIntercept  = Ipj*Ipk
   areaUnderCurve = evaluateIntegral(r,Qj) - evaluateIntegral(r,Gj) - (Qj-Gj)*Qk
   overlayArea = areaUnderCurve+areaLeft+areaBottom-areaIntercept
   return overlayArea

def getOverlay(j,k,r):
   if (j == 0) | (k == 0): return 1.0
   return calculateOverlay(j,k,r)

def foldList(A,n):
    return [A[i:i+n] for i in range(0, len(A), n)]

def generateCircleWeigths( kernelRadius ): # as integer
    a = []
    for j in range(-kernelRadius, kernelRadius+1):
         for k in range(-kernelRadius, kernelRadius+1):
             a.append(getOverlay(abs(j),abs(k),kernelRadius+0.5))
    return foldList(a,int(2*kernelRadius+1))

def generateCOVARWeigths(A,B, kernelRadius ): # as integer
    a = []
    for j in range(-kernelRadius, kernelRadius+1):
         for k in range(-kernelRadius, kernelRadius+1):
             a.append(A[j,k]*B[j,k])
    return foldList(a,int(2*kernelRadius+1))

def generateEastWestWeigths( kernelRadius ): # as integer
    a = []
    for j in range(-kernelRadius, kernelRadius+1):
         for k in range(-kernelRadius, kernelRadius+1):
             a.append(float(k))
    return foldList(a,int(2*kernelRadius+1))

def generateNorthSouthWeigths( kernelRadius ): # as integer
    a = []
    for j in range(-kernelRadius, kernelRadius+1):
         for k in range(-kernelRadius, kernelRadius+1):
             a.append(-float(j))
    return foldList(a,int(2*kernelRadius+1))

def Hadamard(A, B, kernelRadius):
    a = []
    for j in range(0, 2*kernelRadius+1):
         for k in range(0, 2*kernelRadius+1):
             a.append(A[j][k]*B[j][k])
    return foldList(a,int(2*kernelRadius+1))

def sumArray(A, kernelRadius):
    a = 0
    for j in range(0, 2*kernelRadius+1):
         for k in range(0, 2*kernelRadius+1):
             a+=A[j][k]
    return a

def normalizeKernel(A, factor, kernelRadius):
    a = []
    for j in range(0, 2*kernelRadius+1):
         for k in range(0, 2*kernelRadius+1):
            a.append(A[j][k]*factor)
    return  foldList(a,int(2*kernelRadius+1))

def addToKernel(A, factor, kernelRadius):
    a = []
    for j in range(0, 2*kernelRadius+1):
         for k in range(0, 2*kernelRadius+1):
            a.append(float(A[j][k])+float(factor))
    return  foldList(a,int(2*kernelRadius+1))

## PROGRAM BEGINS HERE ##
#                       #
kernelRadius = int(arcpy.GetParameterAsText(0))
spacing = float(arcpy.GetParameterAsText(1))
#
myKernelXX = getCircKernXX(kernelRadius, spacing)
myStringXX = ",".join([','.join(d) for d in [[str(c) for c in b] for b in myKernelXX]])
arcpy.SetParameter(2,myStringXX)
#
myKernelYY = getCircKernYY(kernelRadius, spacing)
myStringYY = ",".join([','.join(d) for d in [[str(c) for c in b] for b in myKernelYY]])
arcpy.SetParameter(3,myStringYY)
#
myKernelXY = getCircKernXY(kernelRadius, spacing)
myStringXY = ",".join([','.join(d) for d in [[str(c) for c in b] for b in myKernelXY]])
arcpy.SetParameter(4,myStringXY)

