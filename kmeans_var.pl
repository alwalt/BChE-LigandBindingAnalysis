#!/usr/bin/perl
use Math::Trig;
use Math::Complex;
# Updated by Eric & Sam on 7/15/11 to account for redundant lines & line jumps
# modified code is followed by original code 


########## global variables ####################
$usage="\nUsage: \.\/vectors_bche_gro453\.pl \[Project\] \[\Inh Letter Code\] \[\# of Clones\] \n\n";
$proj     = $ARGV[0] || die "$usage\n";
$code     = $ARGV[1] || die "$usage\n";
$maxclone = $ARGV[2] || die "$usage\n";
$homedir  = `pwd`; chomp $homedir;
$outfile  = "$homedir/p$proj"."_vector.txt";
chomp $code;
open (OUT, ">$outfile") or die "Can't open $outfile\n";
# size of HC tails, used to define atoms in each chemical group for vector analysis #
#$HCsize = substr($proj,-1);
$HCsize = 4;
#$keyfile = "/home/server/FAHData/BCHE_project1/vector_key.txt";
$keyfile = "/home/server/Thesis/Analysis/noncholine";
open(KEY,">$keyfile");
printf KEY "%-27s %-27s %-27s %-27s %-27s %-27s %-27s %-27s","COPpro","Dvector","PheCOP","Axis","Normal","Theta Phe-O-P","OCH1","OCH2";
close(KEY);
print $code;
$maxclone = $maxclone + 1;

############ iterate through max run & max clone ##########################

$currentclone = 1;
while($currentclone < $maxclone){

	$oldtime = 0;	 
	# define the work directory and go there #
	#$workdir = "/home/server/FAHData/BCHE_project1/PROJ$proj"."/RUN$currentrun"."/CLONE$currentclone"."/";
	$workdir  = "/home/server/data/kmeans/choline_set_1/proj$proj"."/00$currentclone"."/";
	# opendir(DIR,$workdir); 
	chdir $workdir;
	#$tprfile = "/home/server/FAHData/BCHE_project1/PROJ$proj"."/RUN$currentrun"."/CLONE$currentclone"."/frame0.tpr"; 
	$tprfile  = "/home/server/data/kmeans/choline_set_1/proj$proj"."/00$currentclone"."/proj$proj.tpr"; 
        
        #rmsd & gyrate file (added by Walter)
	$rmsdfile = "/home/server/data/kmeans/choline_set_1/proj$proj"."/00$currentclone"."/rmsd.xvg";
	$gyrfile  = "/home/server/data/kmeans/choline_set_1/proj$proj"."/00$currentclone"."/gyrate.xvg";
	$coulfile = "/home/server/data/kmeans/choline_set_1/proj$proj"."/00$currentclone"."/coul.xvg";
        $ljfile = "/home/server/data/kmeans/choline_set_1/proj$proj"."/00$currentclone"."/lj.xvg";


	# This block converts the original data into new concatenated xtc files #
	# This was only needed the first time around, and is thus commented out #
	#$numxtcfiles = `ls \*xtc | wc | awk '{print \$1}'`;
	#@xtcfiles = ();
        #if($numxtcfiles > 0){
	#   for($i=0;$i<$numxtcfiles;$i++){
	#     $xtc = "frame$i".".xtc";	
	#    if($i==0){
	#      push (@xtcfiles, $xtc);
	#    }else{ 
	#      $newxtc = "new$i".".xtc";
	#      $starttime = $i * 1000;	
	#      `trjconv -f $xtc -t0 $starttime -o $newxtc >& /dev/null`;
        #      push (@xtcfiles, $newxtc); 
	#    }  	
        #   }	
           $outxtc  = "proj$proj".".xtc";
	#   `trjcat -o $outxtc -f @xtcfiles >& /dev/null`;


	# make pdb's at all times and run vector calc on them #
	# this system call that makes the pdb files results in extra/missing/misnumbered pdb files
	# for instance, 1771-01-26 has an extra frame_61.pdb that says it's time is only 6000 ps
	# duplicating frames @ certain times (i.e. frame 40 in 1771-01-26) 
	
	system("echo 0 | trjconv -f $outxtc -s $tprfile -sep -o frame.pdb");
	#system("echo 0 | -f proj8202.xtc -s proj8202.tpr -sep -o frame.pdb"); 
	my @pdbfiles = `ls frame\*pdb`;
	$count = scalar(@pdbfiles);	   
	$numframes = $count;
	for($k=1;$k<$numframes;$k++){
	 	 $currentpdb = "frame"."$k".".pdb";
		 $temp = `head -2 $currentpdb | tail -1`;
 	  	 chomp $temp;
  	  	 for($temp) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
		 @templine = split(/ /,$temp);
		 $currenttime = $templine[11];
		 $timetest = $currenttime - $oldtime;
		   print STDOUT "$currentrun  $currentclone  $k  $currentpdb  $currenttime  $oldtime  $timetest\n";
	         $oldtime = $currenttime;	



      if($timetest > 0){  
	open(PDB,"<$currentpdb") || die "\n\tError reading from $currentpdb\n\n";
	if($HCsize==1){
		@PheCOP=(1,4,7,10,13,15);
		@axis = (2,10);
		@norm=(10,4);
		@theta=(1,2,3);
		@dcopi=(2);
		@OCH1=(2,17);
		@OCH2=(2,21);
	}
	if($HCsize==4){
		@PheCOP=(8318,8319,8311,8312,8314,8316);
		@axis = (8321,8312);
		@norm=(8312,8319);
		@theta=(8322,8321,8318);
		@dcopi=(8321);
		@OCH1=(8321,8357);
		@OCH2=(8321,8332);
	}
	if($HCsize==5){
		@PheCOP=(1,27,30,33,36,38);
		@axis = (2,33);
		@norm=(33,27);
		@theta=(1,2,3);
		@dcopi=(2);
		@OCH1=(2,19);
		@OCH2=(2,23);
	}


	#defining variables
	$protx = 0;
	$proty = 0;
	$protz = 0;
	$phenx = 0;
	$pheny = 0;
	$phenz = 0;

	########  store .pdb info for analysis and rewriting   ############
	$i = 0;
	while(defined($origline=<PDB>)) {
	  $line = $origline;
 	  chomp $line;
  	  for($line) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
  	  @input = split(/ /,$line);


  	#### only look at inhibitor coordinates ####
  	  if(@input[0] eq 'ATOM'){
	  	    $atomnum    = @input[1];
   		    $atomname   = @input[2];
   		    $x[$atomnum]= @input[6];
   		    $y[$atomnum]= @input[7];
   		    $z[$atomnum]= @input[8];
      		    if((@input[0] eq 'ATOM')&&(@input[3] ne 'SOL')&&(@input[3] ne '$code')&&(@input[3] ne 'UNK')&&(@input[3] ne 'CL')&&(@input[2] ne 'NA')){ 
			$protx = $protx + $x[$atomnum];
			$proty = $proty + $y[$atomnum];
			$protz = $protz + $z[$atomnum];
			$i++; 
 	    	    }
	  }
	}
	close(PDB);

	$protx = $protx/$i;
	$proty = $proty/$i;
	$protz = $protz/$i;

	#get vector from center of inhibitor to COP
	$dx = $x[@dcopi[0]]-$protx;
	$dy = $y[@dcopi[0]]-$proty;
	$dz = $z[@dcopi[0]]-$protz;
	

	#get phenyl center of position
	$phenx = 0; 
	$pheny = 0; 
	$phenz = 0; 
	foreach $num (@PheCOP){
		$phenx = $phenx + $x[$num];
		$pheny = $pheny + $y[$num];
		$phenz = $phenz + $z[$num];
	}
	$phenx = $phenx/6;
	$pheny = $pheny/6;
	$phenz = $phenz/6;
	

	#find axis vector for phenyl ring
	$axisX = $x[@axis[0]]-$x[@axis[1]];
	$axisY = $y[@axis[0]]-$y[@axis[1]];
	$axisZ = $z[@axis[0]]-$z[@axis[1]];
	
	#find normal vector to ring
	#find vectors from center of ring to C on ring
	$x1 = $x[@norm[0]] - $phenx;
	$y1 = $y[@norm[0]]-$pheny;
	$z1 = $z[@norm[0]]-$phenz;
	$x2 = $x[@norm[1]]-$phenx;
	$y2 = $y[@norm[1]]-$pheny;
	$z2 = $z[@norm[1]]-$phenz;
	$xnorm = ($y1*$z2)-($y2*$z1);
	$ynorm = ($x1*$z2)-($x2*$z1);
	$znorm = ($x1*$y2)-($x2*$y1);


	#find angle theta between phenyl-O-P
	$X1 = $x[@theta[0]]-$x[@theta[1]];
	$Y1 = $y[@theta[0]]-$y[@theta[1]];
	$Z1 = $z[@theta[0]]-$z[@theta[1]];
	$X2 = $x[@theta[2]]-$x[@theta[1]];
	$Y2 = $y[@theta[2]]-$y[@theta[1]];
	$Z2 = $z[@theta[2]]-$z[@theta[1]];
	$dot = ($X1*$X2)+($Y1*$Y2)+($Z1*$Z2);
	$sum1 = ($X1**2)+($Y1**2)+($Z1**2);
	$sum2 = ($X2**2)+($Y2**2)+($Z2**2);
	$mag1 = sqrt($sum1);
	$mag2 = sqrt($sum2);
	$product = $mag1*$mag2;
	$div = $dot/$product;
	$theta = acos($div);
	$theta = rad2deg($theta);
	

	#CH vectors
	$xCH1 = $x[@OCH1[1]]-$x[@OCH1[0]];
	$yCH1 = $y[@OCH1[1]]-$y[@OCH1[0]];
	$zCH1 = $z[@OCH1[1]]-$z[@OCH1[0]];
	$xCH2 = $x[@OCH2[1]]-$x[@OCH2[0]];
	$yCH2 = $y[@OCH2[1]]-$y[@OCH2[0]];
	$zCH2 = $z[@OCH2[1]]-$z[@OCH2[0]];


	open(RMSD,"<$rmsdfile") || die "\n\tError reading from $rmsdfile\n\n";
	open(GYR,"<$gyrfile") || die "\n\tError reading from $gyrfile\n\n";
	open(COUL,"<$coulfile") || die "\n\tError reading from $gyrfile\n\n";
        open(LJ,"<$ljfile") || die "\n\tError reading from $gyrfile\n\n";



	my %tablet;
	while(defined($oline=<RMSD>)) {
        	$myline = $oline;
        	chomp $myline;
        	for ($myline) {  s/^\s+//; s/\s+$//; s/\s+/ /g }
        	@rmsd = split(/ /,$myline);
        	$time = @rmsd[0] + 0; #add + 0
        	$tablet{$time} = @rmsd[1];
        }


        my %tablet2;
        while(defined($oline2=<GYR>)) {
                $myline2 = $oline2;
                chomp $myline2;
                for ($myline2) {  s/^\s+//; s/\s+$//; s/\s+/ /g }
                @gyr = split(/ /,$myline2);
                $time2 = @gyr[0] + 0; #add + 0
                $tablet2{$time2} = @gyr[1];
        }


        my %tablet3;
        while(defined($oline3=<COUL>)) {
                $myline3 = $oline3;
                chomp $myline3;
                for ($myline3) {  s/^\s+//; s/\s+$//; s/\s+/ /g }
                @coul = split(/ /,$myline3);
                $time3 = @coul[0] + 0; #add + 0
                $tablet3{$time3} = @coul[1];
        }

        my %tablet4;
        while(defined($oline4=<LJ>)) {
                $myline4 = $oline4;
                chomp $myline4;
                for ($myline4) {  s/^\s+//; s/\s+$//; s/\s+/ /g }
                @lj = split(/ /,$myline4);
                $time4 = @lj[0] + 0; #add + 0
                $tablet4{$time4} = @lj[1];
        }



	$qr = $oldtime + 0;
	printf OUT "%-10s %-10d %-10d %9.2f",$proj,$currentrun,$currentclone,$currenttime;
	printf OUT "%9.3f %9.3f %9.3f ",$protx,$proty,$protz;
	printf OUT "%9.3f %9.3f %9.3f ",$dx,$dy,$dz;
	printf OUT "%9.3f %9.3f %9.3f ",$phenx,$pheny,$phenz;
	printf OUT "%9.3f %9.3f %9.3f ",$axisX,$axisY,$axisZ;
	printf OUT "%9.3f %9.3f %9.3f ",$xnorm,$ynorm,$znorm;
	printf OUT "%9.3f ",$theta;
	printf OUT "%9.3f %9.3f %9.3f ",$xCH1,$yCH1,$zCH1;
	printf OUT "%9.3f %9.3f %9.3f ",$xCH2,$yCH2,$zCH2;
	printf OUT "%9.3f %9.3f ",$tablet{$qr}, $tablet2{$qr};
	printf OUT "%9.3f %9.3f ",$tablet3{$qr}, $tablet4{$qr};
	printf OUT "%9.3f ",$tablet3{$qr}+$tablet4{$qr};
	print OUT "\n";


  
	############### remove all excess files ##############
    }
   }	
    system("rm *pdb");
    $currentclone++;
}
close(OUT);

