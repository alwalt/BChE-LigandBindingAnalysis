### Subroutines and Modules ###
use Cwd;
use File::Copy;
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

########## global variables ####################
$usage="\nUsage: \.\/kmeans_structural_analysis\.pl \[Project\#\] \[Inhibitor 3-Letter Code\] \n\nKey_file.txt and final_trial.txt should be in current directory \n";
$proj     = $ARGV[0] || die "$usage\n";
$code     = $ARGV[1] || die "$usage\n";
chomp $proj;
chomp $code;
$dir = getcwd;
$totalsim = 0;
$outfile  = "$dir/$proj"."_AA.txt";

###amino acid key###
@amino = (ASN68,ILE69,ASP70,GLN71,PHE73,PRO74,GLY75,PHE76,SER79,MET81,TRP82,ASN83,TYR114,GLY115,GLY116,GLY117,PHE118,GLN119,THR120,TYR128,GLU197,SER198,ALA199,TRP231,ALA277,PRO285,LEU286,SER287,VAL288,GLU325,ALA328,PHE329,TYR332,ASN397,PHE398,TRP430,MET437,HID438,GLY439,TYR440,ILE442);


#$filename = 'test.txt';

#open ($fn, '>', $filename) or die "Yo, can't open the file '$filename' $!";
#print $fn "Cluster\n";
#print $fn "Group\n";
#foreach (@amino){
#	print $fn "$_\n";
#}

#### Read in Inhibitor Key File ####
$keyfile = "noncholine_key_file.txt";
#chomp $keyfile;
open(KI,"<$keyfile")|| die "\n\tError reading from $keyfile\n\n";
while (<KI>){
    if (/Proj $proj/i) {
        $spool = 1;
        next;
    }
    elsif ($_ =~ /Proj\s\d/ && $_ !~ /Proj $proj/) {
        $spool = 0;
    }
    elsif ($spool) {
                undef @singroup;
                $groupline = $_;
                for($groupline) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
                @logline = split(/ /,$groupline);
                $group = "@logline[0]";
                chomp $group;
                $atoms = "@logline[1..$#logline]";
                @singroup = split(/ /,$atoms);
                if ($group eq "PheCOP"){
                        push (@PheCOP,@singroup);
                        undef @singroup;
                }
                elsif ($group eq "OCH1"){
                        push (@OCH1, @singroup);
                        undef @singroup;
                }
                elsif ($group eq "OCH2"){
                        push (@OCH2, @singroup);
                        undef @singroup;
                }
                elsif ($group eq "Phos"){
                        push (@Phos, @singroup);
                        undef @singroup;
                }
    }
}


### Find Cluster Size and Number ###
$log = "proj$proj"."/$proj"."_trial.txt";
chomp $log;
open(LOG,"<$log")|| die "\n\tError reading from $log\n\n";
while (<LOG>) {
    if (/#  group  pop	Centers/i) {
        $spool = 1;
        next;
    }
    elsif (/Clusters \/ Data points \/ P\/R\/C\/T/i) {
        $spool = 0;
	}
    elsif ($spool) {
		$var = $_;
		for($var) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
		@logline = split(/ /,$var);
		$energy{$logline[0]} = $logline[1];
		$totalsim += @logline[1];
    }
}


### Rename/Organize by Cluster Size ###
$k=0;
foreach my $clusternum (sort { $energy{$b} <=> $energy{$a} or $a cmp $b } keys %energy){
	$clusternew{$clusternum} = $k;
	#printf"%-8s %s%s\n", $clusternew{$clusternum};
	$k++;
}


### Make cluster directories ###
mkdir "$dir"."/proj$proj"."/clusters";
#mkdir clusters;
$clusterfolder = "$dir"."/proj$proj"."/clusters";
#print $clusterfolder;

for($j=0; $j<$k; $j++){
	mkdir "$clusterfolder"."/$j";
}


### Look for a pdb's associated cluster number and move ###
seek (LOG, 0, 0);
while (<LOG>) {
    if (/Clusters \/ Data points \/ P\/R\/C\/T/i) {
		$spool = 1;

		next;
	}
    elsif ($spool) {
		$var2 = $_;
		for($var2) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
		@logline2 = split(/ /,$var2);
		$frame = @logline2[24]/100;
		$clone = @logline2[23];
		$oldclusternum = @logline2[0];
		$curlocation = "$dir"."/proj$proj"."/00$clone"."/$clone"."frame$frame".".pdb";
		$newlocation = "$clusterfolder"."/$clusternew{$oldclusternum}";
		#move($curlocation, $newlocation) or die "Error when moving clone:$clone frame$frame.pdb: $!";
    }
}
close(LOG);


### Move into cluster folders and find distances ###
for($m=0; $m<$k; $m++){
	chdir "$clusterfolder"."/$m";
	$cluster = $m;
	my @pdbfiles = `ls *pdb`;
	$count = scalar(@pdbfiles);
	$clustertotal = $count;
	
	### Find distances for each pdb ###
	foreach $pdb(@pdbfiles){
		chomp $pdb;
		$temp = `head -2 $pdb | tail -1`;
		chomp $temp;
		for($temp) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
		@templine = split(/ /,$temp);
		$currenttime = $templine[11];


		### Open/Read pdb line and save protein coordinates ###
		open(PDB,"<$pdb") || die "\n\tError reading from $pdb\n\n";
		$i = 0;
		while(defined($origline=<PDB>)) {
			$line = $origline;
			chomp $line;
			for($line) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
			@input = split(/ /,$line);

			###Save protein coordinates###
			if(@input[0] eq 'ATOM'){
				$atomnum    = @input[1];
				$atomname   = @input[2];
				$resname[$atomnum] = @input[3];
				$resnum[$atomnum]  = @input[5];
				$x[$atomnum]= @input[6];
				$y[$atomnum]= @input[7];
				$z[$atomnum]= @input[8];
				if((@input[0] eq 'ATOM')&&(@input[3] ne 'SOL')&&(@input[3] ne '$code')&&(@input[3] ne 'UNK')&&(@input[3] ne 'CL')&&(@input[2] ne 'NA')){ 
					$protx[$atomnum] = $x[$atomnum];
					$proty[$atomnum] = $y[$atomnum];
					$protz[$atomnum] = $z[$atomnum];
					$i++;
				}
			}
		}
		close(PDB);

		### Find the inhibitor group's atom distance from protein ###

		foreach $num (@PheCOP){
			for ( $i=1; $i<8310; $i++){
				$X = $protx[$i]-$x[$num];
				$Y = $proty[$i]-$y[$num];
				$Z = $protz[$i]-$z[$num];
				$D2 = $X**2 + $Y**2 + $Z**2;
				$D = sqrt($D2);
				if ($D < 4){
					@pheorg=();
					$phe = "$resname[$i]$resnum[$i]";
					for($phe) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
					push @pheline, $phe;
					@pheorg = uniq(@pheline);
				}
			}
		}
		#maybe print pheorg?###
		foreach $num (@Phos){
			#print $num;
			for ( $i=1; $i<8310; $i++){
				$X = $protx[$i]-$x[$num];
				$Y = $proty[$i]-$y[$num];
				$Z = $protz[$i]-$z[$num];
				$D2 = $X**2 + $Y**2 + $Z**2;
				$D = sqrt($D2);
				if ($D < 4){
					@phoorg=();
					$pho = "$resname[$i]$resnum[$i]";
					for($pho) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
					push @pholine, $pho;
					@phoorg = uniq(@pholine);
				}
			}
		}
	
		foreach $num (@OCH1){
			#print $num;
			for ( $i=1; $i<8310; $i++){
				$X = $protx[$i]-$x[$num];
				$Y = $proty[$i]-$y[$num];
				$Z = $protz[$i]-$z[$num];
				$D2 = $X**2 + $Y**2 + $Z**2;
				$D = sqrt($D2);
				if ($D < 3){
					@oc1org=();
					$oc1 = "$resname[$i]$resnum[$i]";
					for($oc1) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
					push @oc1line, $oc1;
					@oc1org = uniq(@oc1line);
				}
			}
		}

		foreach $num (@OCH2){
			for ( $i=1; $i<8310; $i++){
				$X = $protx[$i]-$x[$num];
				$Y = $proty[$i]-$y[$num];
				$Z = $protz[$i]-$z[$num];
				$D2 = $X**2 + $Y**2 + $Z**2;
				$D = sqrt($D2);
				if ($D < 3){
					@oc2org=();
					$oc2 = "$resname[$i]$resnum[$i]";
					for($oc2) {  s/^\s+//; s/\s+$//; s/\s+/ /g; }
					push @oc2line, $oc2;
					@oc2org = uniq(@oc2line);
				}
			}
		}

		push(@phecluster,@pheorg);
		@pheorg=();
		undef @pheline;

		push(@phocluster,@phoorg);
		@phoorg=();
		undef @pholine;

		push(@oc1cluster,@oc1org);
		@oc1org=();
		undef @oclline;

		push(@oc2cluster,@oc2org);
		@oc2org=();
		undef @oc2line;
	}

	####Finding total interaction time###
	
        $phecount{$_}++ foreach @phecluster;
        $phocount{$_}++ foreach @phocluster;
        $oc1count{$_}++ foreach @oc1cluster;
        $oc2count{$_}++ foreach @oc2cluster;
	
	#my %grades;
	foreach $ami(@amino){
		$cphe = $phecount{$ami};
		$cpho = $phocount{$ami};
		$coc1 = $oc1count{$ami};
		$coc2 = $oc2count{$ami};
		if ($cphe == 0 && $cpho == 0 && $coc1 == 0 && $coc2 == 0){
			$grades{$cluster}{$ami} = 0;
		}
		elsif ($cphe >= $cpho && $cphe >= $coc1 && $cphe >= $coc2){
			$grades{$cluster}{$ami} = "Phe";
		}
		elsif ($cpho >= $coc1 && $cpho >= $coc2){
			$grades{$cluster}{$ami} = "Phos";
		}
		elsif ($coc1 >= $coc2){
			$grades{$cluster}{$ami} = "Alk-L";
		}
		else{
			$grades{$cluster}{$ami} = "Alk-R";
		}
	}


        @phecluster = ();
	@phocluster = ();
	@oc1cluster = ();
	@oc2cluster = ();
	undef %phecount; 
	undef %phocount;
	undef %oc1count;
	undef %oc2count;
	

	

}

open (OUT, ">$outfile") or die "Can't open $outfile\n";
printf OUT "%-8s", "Clus";
 
foreach $amis(@amino){
	printf OUT "%-8s", $amis;
}

printf OUT "\n";

for($z=0; $z<$k; $z++){
	printf OUT "%-8s", $z;
	foreach $amis(@amino){
		printf OUT "%-8s", $grades{$z}{$amis};
	}
	printf OUT "\n";
}


close OUT;
