dir2=$(pwd)

for folder in $dir2/*; do
  [ -d "$folder" ] && cd "$folder" && ./auto-plot-script.py
done
cd $dir2

find . -name \*.png -exec cp {} $dir2/plots \;
