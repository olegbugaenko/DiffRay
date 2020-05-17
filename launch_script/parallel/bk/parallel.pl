#!/usr/bin/perl -w
$iage = @ARGV[0];
$minSector = @ARGV[1];
$maxSector = @ARGV[2];
$n = 4;
$isector = $minSector;
mkdir "../../result/ngc/ism0p05/parallel/$iage";
while($isector <= $maxSector)
{
	$nProcesses = 0;
	for($iprocess = 0; $iprocess < 4; $iprocess++)
	{
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

			$phi = 3.14*($rSector-0.5)/(2*$maxSector);
			print COMMANDS "print directory ../../result/ngc/ism0p05/parallel/$iage\n";
			print COMMANDS "input directory ../../data/ngc1569_ism0p05_new\n";
			print COMMANDS "input database ../../data/database\n";
			print COMMANDS "input bands ../../data/bands\n";
			print COMMANDS "distance file Emis_Lines_SectorNo${rSector}_Age$iage.00Myr.dat 2\n";
			print COMMANDS "object age $iage\n";
			print COMMANDS "integration precision 0.005\n";
			print COMMANDS "integration maxiterations 12\n";
			print COMMANDS "apperture full\n";
			print COMMANDS "object phi ${phi}\n";
			print COMMANDS "object theta 0.000\n";
			print COMMANDS "calc abundancescalc linesprint bands\n";
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