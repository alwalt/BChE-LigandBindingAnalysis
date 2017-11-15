#!/usr/bin/perl


########## global variables ####################
$usage="\nUsage: \.\/vectors_bche_gro453\.pl \[Project\] \[\# of Clones\] \n\n";
$proj     = $ARGV[0] || die "$usage\n";
$maxclone = $ARGV[1] || die "$usage\n";
$homedir  = `pwd`; chomp $homedir;

############ iterate through max run & max clone ##########################

$currentclone = 1;
while($currentclone <= $maxclone){
	$workdir = "$homedir"."/proj$proj"."/00$currentclone"."/";
	chdir $workdir;
	$tprfile = "$homedir"."/proj$proj"."/00$currentclone"."/proj$proj.tpr"; 
	$ndxfile = "$homedir"."/proj$proj"."/00$currentclone"."/proj$proj.ndx";
        $outxtc  = "$homedir"."/proj$proj"."/00$currentclone"."/proj$proj".".xtc";
	$newframe = "$currentclone"."frame";	

	system("echo 24 | trjconv -f $outxtc -s $tprfile -sep -n $ndxfile -o $newframe.pdb");

 	$currentclone++;
}
close(OUT);

