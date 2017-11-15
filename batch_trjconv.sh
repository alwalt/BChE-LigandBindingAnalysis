#dont forget to copy prepvectorskmeans into folders

dir1=$(pwd)
rmsf='/home/server/Desktop/rmsf_analysis'

for i in {06..20}; do
	cd proj82"$i";
	echo "1" | trjcat -f */*xtc -o $rmsf/proj82"$i"_cat -n proj82"$i".ndx -b 10000 -cat -nosort;
	echo "1" | g_rmsf -f $rmsf/proj82"$i"_cat.xtc -s proj82"$i".gro -n proj82"$i".ndx -o $rmsf/proj82"$i" -res;
	cd $dir1;
done



