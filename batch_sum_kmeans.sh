#dont forget to copy prepvectorskmeans into folders

dir1=$(pwd)

for i in {02..20}; do
	cp -f summarize_kmeans_results.pl proj82"$i";
	cd proj82"$i";
	perl summarize_kmeans_results.pl 100 p82"$i"_kmeans_results.txt;
	cd $dir1
done



