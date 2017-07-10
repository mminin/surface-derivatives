Repo for standalone scripts for making surface derivatives.
Should be easier to work with for batch operations, etc.
Users of ArcGIS may find this useful, as it does not require using Arc licence, 
nor is it bound by limitations of 32-bit address space, so you can take advantage of 64-bit processor,
which would make calculations a lot faster.
Also this code is a lot simpler and shorter than ArcGIS plugin version.

Perhaps the most importaint dependency here is Python-GDAL bindings in osgeo,
install thouse via http://trac.osgeo.org/osgeo4w/wiki

If working from Windows, you should use "OSGeo4W Shell" to run this script.

If default installation does not work, use advanced installation, and make sure to select python libraries (skipped by default), particularly python-gdal-dev and scipy. Scipy provides handles to ITK subroutines, which (i suspect) is what does the convolution.

To run this script, say from either Linux terminal or OSGeo4W shell:

    python FD.py yourDTMFile.tiff [kernelSizeX=3 kernelSizeY=3 gridSpacing=50]
For instance for file Mikhail_test.tif with 3x3 kernel and grid spacing of 50 m/px, say:

    python FD.py Mikhail_test.tif
Which is the same as saying

    python FD.py Mikhail_test.tif 3 3 50
For larger kernel you can say

    python FD.py Mikhial_test.tif 7 7 50
