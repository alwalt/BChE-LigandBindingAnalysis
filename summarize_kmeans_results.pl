#! /usr/bin/perl
$startcount = 0;


########## global variables ####################
$usage="\nUsage: \.\/summarize_kmeans_results.pl  \[Number of trials] [Output filename]\n\n";
$numtrials  = $ARGV[0] || die "$usage\n";
$outfile    = $ARGV[1] || die "$usage\n";

$numtrials = $numtrials;

for($i=1;$i<=$numtrials;$i++){
	$clustercount{$i} = 0;
}

for($i=1;$i<=$numtrials;$i++){

  $numclusters = 0;
  $infile = "trial."."$i".".kmeans.100.txt";

  `head -120 $infile > temp.txt`;

  open (INP, "<temp.txt") or die "Can't open $infile\n";
  while ($line = <INP>){
        chomp ($line);
        foreach($line) { s/^\s+//;s/\s+$//; s/\s+/ /g; }
        @info = split(/ /,$line);

	if(@info[0] eq "Clusters"){ $startcount = 0; }

	if($startcount == 1){ $numclusters++; }

	if(@info[1] eq "group"){ $startcount = 1; }

  }
  $clustercount{$i} = $numclusters;   
}


close(INP); 
`rm temp.txt`; 
`touch temp.txt`;

open(TMP,">temp.txt");
for($i=1;$i<=$numtrials;$i++){
  	if($clustercount{$i} > 0){ 
		printf TMP "%3d clusters in trial %3d\n",$clustercount{$i},$i;	
	}	
}
close(TMP);

`sort -nr temp.txt > $outfile`;
`rm temp.txt`;

