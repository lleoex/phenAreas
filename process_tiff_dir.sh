#!/bin/bash

tif_root=/home/gonchukov-lv/src/phenomenaAreas/tiffs
script_root=/home/gonchukov-lv/src/phenomenaAreas

dir=$1

if [[ -d $dir ]]
then
    for file in $dir/*.tiff
    do
	if [[ -f $file ]]
	then
    	    echo `date` $file
    	    $script_root/prepare_for_tiling.sh $file > /dev/null
    	    gdal2tiles.py --zoom=5-10 -r near --processes=8 $file
    	    rm -f $file
    	fi
    done
fi
                                            