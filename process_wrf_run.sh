#!/bin/bash



# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/gonchukov-lv/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/gonchukov-lv/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/gonchukov-lv/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/gonchukov-lv/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate tiles                                

dt_str=$1

script_root=/home/gonchukov-lv/src/phenomenaAreas
tif_root=/home/gonchukov-lv/src/phenomenaAreas/tiffs
result_dir=/mnt/home/export/wmstiles
http_result_dir="http://95.173.159.38/primugms/~lgonchukov/wmstiles"

for dom in 1 2 3
do
    echo dom $dom
    python $script_root/main.py $dom $dt_str > $script_root/d0${dom}.out
done

## make tiles from rgba tifs
echo $tif_root
for dom in wrf1 wrf3
do
    for dir in $tif_root/$dom/*
    do
	if [[ -d $dir ]]
	then
	    cd $dir/$dt_str
	    $script_root/process_tiff_dir.sh $dir/$dt_str > $dir.log &
	fi
    done
    echo "waiting while tiling"
    wait
done
                    
##publish tiles to 10.11.25.67 via ssh
rm -f $result_dir/last.list
for dom in wrf1 wrf3
do
    for dir in $tif_root/$dom/*
    do
    if [[ -d $dir ]]
    then
        cd $dir
	echo $dir
        basedir=`basename $dir`
        rm -f $dir/$dt_str/*.tiff
        cp $script_root/cbar_$basedir.png $dir/$dt_str/
        tar -czf $dt_str.tar.gz ./$dt_str
        echo $dt_str.tar
        tag_dir=$result_dir/$dom/$basedir
	if [[ ! -d $tag_dir ]]
        then
            mkdir $tag_dir
        fi
        cp $dt_str.tar.gz $tag_dir
        echo "$http_result_dir/$dom/$basedir/$dt_str.tar.gz" >> $result_dir/last.list
        echo " out=${dom}Z${basedir}Z${dt_str}.tar.gz" >> $result_dir/last.list
#        echo ssh -t dataloader@10.11.25.67 "'cd /var/data/wmstiles/wrf5/$basedir ; tar -xf /var/data/wmstiles/wrf5/$basedir/$dt.tar'"
#        ssh -t dataloader@10.11.25.67 "cd /var/data/wmstiles/wrf5/$basedir ; tar -xf /var/data/wmstiles/wrf5/$basedir/$dt.tar"
#        rm $dt.tar
#        rm /mnt/tiles_data/wmstiles/wrf5/$basedir/$dt.tar
#        rm -rf $dir/$dt
#        rm -f /mnt/tiles_data/wmstiles/wrf5/$basedir/list.txt
#        #for d in `ls -d /mnt/tiles_data/wmstiles/wrf5/$basedir/????-??-??-??/`; do echo `basename $d` >> /mnt/tiles_data/wmstiles/wrf5/$basedir/list.txt; done
#        echo `ls /mnt/tiles_data/wmstiles/wrf5/$basedir/` > /mnt/tiles_data/wmstiles/wrf5/$basedir/list.txt
    fi
done                                                                                                                                            
done
