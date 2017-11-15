for i in proj*/; do
    cp -rf make_rms_rg_xvgs_index_mod.py "$i"
    #cp auto-plot-script.py "$i"
done

dir1=$(pwd)

for folder in $dir1/*; do 
  [ -d "$folder" ] && cd "$folder" && ./make_rms_rg_xvgs_index_mod.py
done
cd $dir1


