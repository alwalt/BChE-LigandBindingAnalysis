#! /usr/bin/perl
# to cluster in 15D, remove these fields: 5,6,7 (COPpro); 11,12,13 (COPphe); 20 (POP angle)
# use spaces between chosen columns after flag
# also move the proj/run/clone/time fields to end of lines
@tocut = '';

########## global variables ####################
$usage="\nUsage: \.\/delete_desired_fields.pl  [Input File\]  [cols to remove]\n
\tBe sure to start your counting at 1, not at 0 (and use space delimeter for cols to remove)\n\n";
$inpfile   = $ARGV[0] || die "$usage\n";
$numinput  = $#ARGV; 

for($i=1;$i<=$numinput;$i++){
  chomp $ARGV[$i];
  $delme = $ARGV[$i] - 1;
  $j = $i - 1;
  @tocut[$j] = $delme;
}
$numcuts = scalar @tocut;

open (INP, "<$inpfile") or die "Can't open $inpfile\n";
while ($line = <INP>){
  chomp ($line);
  for($line) { s/^\s+//;s/\s+$//; s/\s+/ /g; }
  @lines = split(/ /,$line);
  $numcols = scalar @lines;

  for($i=0;$i<$numcols;$i++){

    $total = 0;
    # if($i is not found in @tocut){
    for($j=0;$j<$numcuts;$j++){
      if($i == @tocut[$j]){ $total++; } 
    }    

    if($total == 0){ 
	if($i > 3){ printf STDOUT "%10.3f",@lines[$i]; }
    }
  }
  printf STDOUT "%6d",@lines[0];
  printf STDOUT "%6d",@lines[1];
  printf STDOUT "%6d",@lines[2];
  printf STDOUT "%6d",@lines[3];
  print STDOUT "\n";
}
close(INP);

