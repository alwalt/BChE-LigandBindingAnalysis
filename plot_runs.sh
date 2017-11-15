for i in proj*/; do
    cp make_rms_rg_xvgs_index_mod.py "$i"
    cp auto-plot-script.py "$i"
done

dir1=$(pwd)
for folder in ~/Thesis/kmeans/noncholine/*; do 
  [ -d "$folder" ] && cd "$folder";
	for i in 00*/; do  
		tail -1 "$i"/rmsd.xvg;
		done
done
cd $dir1

#dir2=$(pwd)
#for folder in ~/Thesis/kmeans/noncholine/*; do
#  [ -d "$folder" ] && cd "$folder" && ./auto-plot-script.py
#done
#cd $dir2

#cp **/*.png ~/Thesis/kmeans/noncholine/images/

