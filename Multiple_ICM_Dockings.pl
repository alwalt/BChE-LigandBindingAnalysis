#! /usr/bin/perl
#use Cwd qw(abs_path);

#ICM set up project with name "Trail" then set up Receptor. Create ligand in mol file. Then set-up receptors file and ligand files in the icm_home.
#Run this script to begin the docking process. 
#The script with input information: ligandname 	number_of_simulation

$project=$ARGV[0];
$ligand=$ARGV[1];
$max_run = $ARGV[2];
$thorough=$ARGV[3];

#Start time
$start = gmtime(time());
##########################

#$script_home = abs_path($0); #with this we can have fullpathway
$script_home = $ENV{'PWD'};
$data_home = $script_home.'/'.$project.'/'.$ligand.'/';

#$data_home ="/home/walter/Desktop/ICMdockings/$project/$ligand/"; #Adjust this path to match your output files folder
$icm_home = "/home/server/icm-3.7-2b/"; #******************************************
#$icm_home = $ENV{'ICMHOME'}.'/'; # $ICMHOME environment varible stores icm pathway.

print $data_home."\n";
system("mkdir -p -v $data_home");

for($i=1; $i<=$max_run; $i++) {                      
   chdir $data_home;
   $obfile="rm $ligand"."_dock$i.ob";
   if(-e $ligand."_dock$i.ob") { system($obfile); }
   chdir $icm_home;
   system("./icm64 _dockScan $project input=$ligand.mol -s confs=50 thorough=$thorough outdir=$data_home");    
   chdir $data_home; 
   $obfile="mv $project"."_$ligand"."1.ob $ligand"."_dock$i.ob";
   system($obfile);
   $date = `date`;
   print "Docking $i complete ... $date\n\n";
}            

#Create ICM script
open(ICM,'>',$icm_home."$project") || die "Please give me output filename $!"; #adjust the ICMscript 
print ICM "#!$icm_home"."icm -s\n";
print ICM "for i=1, $max_run\n"; 
print ICM "s_obname= \"$data_home$ligand"."_dock\"+i+\".ob\";\n";
print ICM "s_sdfname= \"$data_home$ligand"."_dock\"+i+\".sdf\";\n"; 
#print ICM 'read  stack s_cnfname'."\n";
print ICM 'read object s_obname'."\n";
print ICM 'load stack a_'."\n";
print ICM 'write Energy(stack) s_sdfname'."\n"; 
print ICM "endfor\n";
print ICM "quit\n";
close(ICM)||die $!;
#Ending creating ICM script

#running the command to generate sdf files that were used to created a file log file
chdir $icm_home;
system("./icm64 -s $project"); 
#end running the command to generate sdf files

#create the final log file
chdir $data_home;

open(W,'>',"temp.log") || die "Please give me output filename $!";

for($index=1;$index<=$max_run;$index++)
{
   open(Read,'<',$ligand."_dock$index.sdf")||die $!;
   while($line=<Read>)
   {
      chomp($line);
      foreach($line) { s/^\s+//;s/\s+$//; s/\s+/ /g; }
      my @temp=split(' ',$line);
      $ICMscore = $temp[0];
      close(Read)||die $!;
   }
   print W $ligand."_dock$index\t$ICMscore\n";
   
}
close(W)||die $!;
$outputName=$project."_".$ligand.".log";
`sort -n -k2 temp.log > $outputName`;
`rm temp.log`;
`rm \*.sdf`;
chdir $icm_home;
system("rm $project"); # Remove ICM script file.

$end = gmtime(time());
print "Start Process at $start\n";
print "End Process at   $end\n\n\n";
