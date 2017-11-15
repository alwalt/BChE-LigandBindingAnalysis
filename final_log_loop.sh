#dont forget to copy prepvectorskmeans into folders

dir1=$(pwd)

PROJ  = $1
TRIAL = $2

for PROJ in proj"$PROJ"; do
	cd proj"$PROJ";
	cp kmeans_trials/trial."$TRIAL".kmeans.100.txt ../p"$PROJ"_centers.txt;	
	perl Kmeans_clustering_BChE.pl -data p"$PROJ".txt -k 20 -name p"$PROJ"_final_fit -nume 17 -clu p"$PROJ"_centers.txt -hold;
	awk '{ if ($17<=.4) print $0, "       0"; else print $0, "       1"}' "$PROJ".log > "$PROJ"_io.log;
	cd proj82"$i";
	perl summarize_kmeans_results.pl 100 p82"$i"_kmeans_results.txt;
	cd $dir1
done



