import arcpy       
import math

## ## ## Square Kernel
def getRegKernEW(kernelRadius, spacing):
    GradEWraw     = generateEastWestWeights(                   kernelRadius)
    GradEWsq      = Hadamard(GradEWraw,         GradEWraw,     kernelRadius)
    gradEWnorm    = sumArray(GradEWsq,                         kernelRadius)
    return normalizeKernel(GradEWraw, 1/(spacing*gradEWnorm),  kernelRadius)

def getRegKernNS(kernelRadius, spacing):
    GradNSraw     = generateNorthSouthWeights(                 kernelRadius)
    GradNSsq      = Hadamard(GradNSraw,         GradNSraw,     kernelRadius)
    gradNSnorm    = sumArray(GradNSsq,                         kernelRadius)
    return normalizeKernel(GradNSraw, 1/(spacing*gradNSnorm),  kernelRadius)

## ## ## Circular Kernel
def getCircKernEW(kernelRadius, spacing):
    GradEWraw     = generateEastWestWeights(                   kernelRadius)
    CircleWeights = generateCircleWeights(                     kernelRadius)
    GradEWcirc    = Hadamard(GradEWraw,         CircleWeights, kernelRadius)
    GradEWsq      = Hadamard(GradEWraw,         GradEWraw,     kernelRadius)
    GradEWsqCirc  = Hadamard(GradEWsq,          CircleWeights, kernelRadius)
    gradEWnorm    = sumArray(GradEWsqCirc,                     kernelRadius)
    return normalizeKernel(GradEWcirc, 1/(spacing*gradEWnorm), kernelRadius)

def getCircKernNS(kernelRadius, spacing):
    GradNSraw     = generateNorthSouthWeights(                 kernelRadius)
    CircleWeights = generateCircleWeights(                     kernelRadius)
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

def generateCircleWeights( kernelRadius ): # as integer
    a = []
    for j in range(-kernelRadius, kernelRadius+1):
         for k in range(-kernelRadius, kernelRadius+1):
             a.append(getOverlay(abs(j),abs(k),kernelRadius+0.5))
    return foldList(a,int(2*kernelRadius+1))

def generateEastWestWeights( kernelRadius ): # as integer
    a = []
    for j in range(-kernelRadius, kernelRadius+1):
         for k in range(-kernelRadius, kernelRadius+1):
             a.append(float(k))
    return foldList(a,int(2*kernelRadius+1))

def generateNorthSouthWeights( kernelRadius ): # as integer
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

kernelRadius = int(arcpy.GetParameterAsText(0))
spacing = float(arcpy.GetParameterAsText(1))
isCircularShape = bool(arcpy.GetParameter(2))
if isCircularShape:
    myKernelX = getCircKernEW(kernelRadius, spacing)
    myKernelY = getCircKernNS(kernelRadius, spacing)
else:
    myKernelX = getRegKernEW(kernelRadius, spacing)
    myKernelY = getRegKernNS(kernelRadius, spacing)

myStringX = ",".join([','.join(d) for d in [[str(c) for c in b] for b in myKernelX]])
arcpy.SetParameter(3,myStringX)
myStringY = ",".join([','.join(d) for d in [[str(c) for c in b] for b in myKernelY]])
arcpy.SetParameter(4,myStringY)   

