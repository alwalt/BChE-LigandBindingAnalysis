#! /usr/bin/perl


$projectname = $ARGV[0];
$max_run = $ARGV[1];
$thorough = $ARGV[2];
$script_home = $ENV{'PWD'};
$results =$ENV{'PWD'}.'/'.$projectname.'/';

%ligand = (1 => "l_tyrosine", 2 => "d_tyrosine", 3 => "l_tryptophan", 4 => "d_tryptophan");

$NoLigand = scalar(keys %ligand);

for($j=1;$j<=$NoLigand;$j++)
{
   chdir $script_home;
   $command= "perl Multiple_ICM_Dockings.pl $projectname $ligand{$j} $max_run $thorough";
   system($command);
}

open(W,'>',$results."Total_Results.log") || die "Please give me output filename $!";

for($f=1;$f<=$NoLigand;$f++)
{  

   open(Readlog,'<', "$results$ligand{$f}".'/'.$projectname.'_'.$ligand{$f}.'.log') || die $!;
   while($line=<Readlog>)
   {
      chomp($line);
      foreach($line) { s/^\s+//;s/\s+$//; s/\s+/ /g; }
      my @temp=split(' ',$line);
      $ICMscore = $temp[1];
      close(Readlog)||die $!;
   }
   printf W "%-15s\t%15s\n",$ligand{$f},$ICMscore;
   
}
close(W)||die $!;


