#!/bin/bash

file=$1
filename=${file%.*f}
rm -f ${filename}_500.tiff
rm -f ${filename}_84.tiff
gdalwarp -tr 500 500 $1 ${filename}_500.tiff
gdalwarp -t_srs EPSG:4326 -dstalpha -dstnodata "0 0 0" ${filename}_500.tiff ${filename}_84.tiff
rm -f ${filename}_500.tiff
mv ${filename}_84.tiff $1
