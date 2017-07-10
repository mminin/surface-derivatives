Repo for standalone scripts for making surface derivatives.
Should be easier to work with for batch operations, etc.
Users of ArcGIS may find this useful, as it does not require using Arc licence, 
nor is it bound by limitations of 32-bit address space, so you can take advantage of 64-bit processor,
which would make calculations a lot faster.
Also this code is a lot simpler and shorter than ArcGIS plugin version.

Perhaps the most importaint dependency here is Python-GDAL bindings in osgeo,
install thouse via http://trac.osgeo.org/osgeo4w/wiki
Use advanced installation, and make sure to select python libraries (skipped by default), particularly python-gdal-dev and scipy.

Scipy provides handles to ITK subroutines, which (i suspect) is what does the convolution.
