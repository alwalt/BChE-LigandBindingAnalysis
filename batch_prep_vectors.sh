#dont forget to copy prepvectorskmeans into folders

dir1=$(pwd)

for i in {02..20}; do
	cd proj82"$i";
	perl prep_vectors_for_Kmeans.pl p82"$i"_vector.txt 5 6 7 11 12 13 20 > p82"$i".txt;
	awk '($24 >= 1000) {print $0}' p82"$i".txt > p82"$i"eq.txt; 
	#python Multiple_kmeans.py p82"$i"eq.txt;
	cd $dir1
done



