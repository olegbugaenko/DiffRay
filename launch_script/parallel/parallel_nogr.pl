#!/usr/bin/perl -w
$iage = @ARGV[0];
$minSector = @ARGV[1];
$maxSector = @ARGV[2];
$n = 1;
$isector = $minSector;
mkdir "../../result/ngc/nograins_parallel/$iage";
while($isector <= $maxSector)
{
	$nProcesses = 0;
	for($iprocess = 0; $iprocess < $n; $iprocess++)
	{
		if($rSector > $maxSector) 
		{
			break;
		}
		$rSector = $isector + $iprocess;
		print "Start process $iprocess ($rSector)\n";
		my $pid = fork;
		if(not defined $pid)
		{
			print "Unable to fork.\n";
			next;
		}
		if($pid)
		{
			$nProcesses++;
			printf("Entering $pid for $rSector ($nProcesses)\n");
			open(COMMANDS, "> commands$iage.$rSector.ini")||
				die "CANNOT OPEN COMMANDS FILE\n";

			$theta = 3.14*($rSector-0.05)/(2*$maxSector);
			print COMMANDS "print directory ../../result/ngc/nograins_parallel/$iage\n";
			print COMMANDS "print loglevel debug\n";
			print COMMANDS "input directory ../../data/ngc1569_nograin\n";
			print COMMANDS "input database ../../data/database\n";
			print COMMANDS "object geometry cones\n";
			print COMMANDS "integration maxiterations 1\n";
			print COMMANDS "integration precision 0.01\n";
			print COMMANDS "integration predictive on\n";
			print COMMANDS "calculation opacity off\n";
			print COMMANDS "calculation grains off\n";
			print COMMANDS "calculation continuum on\n";
			print COMMANDS "calculation bands off\n";
			# print COMMANDS "input bands ../../data/bands\n";
			print COMMANDS "distance file Emis_Lines_SectorNo${rSector}_Age$iage.00Myr.dat 2\n";
			print COMMANDS "object age $iage\n";
			print COMMANDS "apperture full\n";
			print COMMANDS "object phi 0.000\n";
			print COMMANDS "object theta ${theta}\n";
			#print COMMANDS "calc abundancescalc linesprint bands\n";
			print COMMANDS "print fluxes Sector${rSector}_Age$iage.00Myr_Fluxes.dat\n";
			close (COMMANDS);
			
		}
		else
		{
			system "./DiffRay3D.exe commands$iage.$rSector.ini";
			print "Process ($$) exiting\n";
			exit;
		}
	}
	for (1 .. $nProcesses) {
	   my $pid = wait();
	   print "Parent saw $pid exiting\n";
	   
	}
	print "Finished iter $nProcesses\n";
	$isector += $nProcesses;
}