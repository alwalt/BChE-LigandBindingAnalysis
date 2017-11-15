for i in /home/server/data/kmeans/choline_set_1/proj*/; do
    cp Summarize_kmeans_results.pl "$i"
done

dir1=$(pwd)
for folder in /home/server/data/kmeans/choline_set_1/*; do 
  [ -d "$folder" ] && cd "$folder" && perl Summarize_kmeans_results.pl
done
cd $dir1


