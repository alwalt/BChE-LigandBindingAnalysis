#!/usr/bin/perl
@input = ();

########## global variables ####################
$usage="\nUsage: \.\/vector_file_checkup\.pl \[file name\]\n\n";
$proj = $ARGV[0] || die "$usage\n";

open(KEY,"<$proj");
$iter = 0;
@input = '';
@inputold = '';

while(defined($origline=<KEY>)) {
  $line = $origline;
  chomp $line;
  for($line) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
  @input = split(/ /,$line);
  $arraysize = scalar(@input);
  $pron = @input[0];
  $runn = @input[1];
  $clon = @input[2];
  $timn = @input[3];
  $c7n  = @input[7];
  $c13n = @input[13];
  $c15n = @input[15];

  # look for identical times #  
  if($timn == $timo){ print STDOUT "$pron $runn $clon $timn redundant time\n"; }

  # compare columns 7, 13 & 15 to determine if identical values were used #
  if(($c7n == $c7o)&&($c13n == $c13o)&&($c15n == $c15o)){  print STDOUT "$pron $runn $clon $timn identical values\n";  }

  # look for time jumps #
  $diff = $timn - $timo;
  if($diff > 100.1){
	if($clon == $cloo){	
		print STDOUT "$pron $runn $clon $timn timejump\n";
	}
  }

  $proo = @input[0];
  $runo = @input[1];
  $cloo = @input[2];
  $timo = @input[3];
  $c7o  = @input[7];
  $c13o = @input[13];
  $c15o = @input[15];

  @inputold = @input;
  $iter++;
}
close(KEY);
