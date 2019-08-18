#include "test.h"
#include "reader.h"
#include "lines.h"
#include "geometry.h"
#include "app.h";
#include "continuum.h";
#include "geom3d.h";
#include "solver.h";
#include "diffray.h";
#include "abund.h";
#include "grain_temp.h";
#include "physics.h"

int main(int argc, char **argv)
{

	if(argc > 1) {
		sprintf(App::commands_path, argv[1]);
	}
	else
	{
		sprintf(App::commands_path, "commands.ini");
	}

	printf("Reading commands from %s\n\n", App::commands_path);

	App::readCommands();

	//printf("!!!!!!!!%s\n", App::output_dir);

	printf("output: %s\ninput: %s\nage: %.2lf; app-width: %lfx%lf\n", App::output_dir, App::input_dir, App::age, App::phi_width, App::theta_width);

	
	//read geometry

	char path[500];

	

	for(int i=0;i<App::nSectors;i++)
	{
		sprintf(path,"%s/Emis_Cont_SectorNo%d_Age%.2lfMyr.dat",App::input_dir,i+1,App::age);
		CGeometry::readRadiuses(path,i);
	}

	for(int i=0;i<CGeometry::nSectors;i++)
	{
		printf("nSectors[%d] = %d (%le - %le)\n", i, CGeometry::nRadiuses[i],CGeometry::outer_radius[i][0],CGeometry::outer_radius[i][CGeometry::nRadiuses[i]]);
	}
	
	//read lines
	CLine::readen = 0;
	printf("READING LINES\n");

	double **resarr;
	char ovrflname[500];
	char flname[500];
	char opflname[500];
	char trflname[500];
	char abflname[500];
	char inflname[500];
	char intempfile[500];

	printf("dbentry\n");
	CLine::readDatabase();

	//CLine::init();

	//double ***CLine::emits;

	if(App::CalcLines)
	{
		for(int i=0;i<CGeometry::nSectors;i++)
		{
			sprintf(flname, "%s/Emis_Lines_SectorNo%d_Age%.2lfMyr.dat",App::input_dir,i+1,App::age);
			printf("Lines %d\n",i+1);
			CLine::readLines(flname,i);
			printf("\nReaden\n");
			//printf("Emts %d = %lf\n", i, CLine::emits[i][113][0]);
			//printf("Prev %d = %lf\n", i, CLine::emits[0][113][0]);
			
		}
		//printf("READEN: %le - %le\n",CLine::emits[1][113][0],CLine::emits[1][83][24]);
	}
	

	sprintf(flname, "%s/Continuum%d_SB99_Age%.2lfMyr.dat",App::input_dir,1,App::age);
	printf("Reading mesh\n");
	CContinuum::readMesh(flname);
	printf("Readen mesh\n");
	if(CContinuum::cellCount<10)
	{
		sprintf(flname, "%s/Continuum%d_Last_Age%.2lfMyr.dat",App::input_dir,1,App::age);
		CContinuum::readMesh(flname);

	}
	printf("Readen mesh 2\n");
	printf("Mesh count: %d\n", CContinuum::cellCount);
	printf("%le - %le - %le\n", CContinuum::anu[0], CContinuum::anu[1089], CContinuum::anu[CContinuum::cellCount-1]);


	for(int i=0;i<CGeometry::nSectors;i++)
	{
		sprintf(ovrflname, "%s/Overview%d_Age%.2lfMyr.dat",App::input_dir,i+1,App::age);
		sprintf(flname, "%s/Emis_Cont_SectorNo%d_Age%.2lfMyr.dat",App::input_dir,i+1,App::age);
		sprintf(opflname, "%s/Opac_Cont_SectorNo%d_Age%.2lfMyr.dat",App::input_dir,i+1,App::age);
		sprintf(trflname, "%s/Emis_Trans_SectorNo%d_Age%.2lfMyr.dat",App::input_dir,i+1,App::age);
		sprintf(abflname, "%s/Abund%d_Age%.2lfMyr.dat",App::input_dir,i+1,App::age);
		sprintf(inflname, "%s/Continuum%d_SB99_Age%.2lfMyr.dat",App::input_dir,i+1,App::age);
		sprintf(intempfile, "%s/Grain_Temp%d_Age%.2lfMyr.dat",App::input_dir,i+1,App::age);

		Physics::readPhysics(ovrflname, i);


		if(App::CalcCont)
		{
			printf("Reading cont %d\n", i);
			CContinuum::readEmissivity(flname,i);

		}
		
		if(App::CalcOpac)
		{
			printf("reading opacity...\n");
			CContinuum::readOpacity(opflname,i);
		}
		
		if(App::CalcLines && App::CalcCont)
		{
			printf("reading trans...\n");
			CContinuum::readTransitions(trflname,i);
			printf("transitions_readen!");
		}

		if(App::CalcAbund)
		{
			printf("reading abunds...\n");
			Abund::readAbunds(abflname,i);
		}


		if(App::CalcCont)
		{
			CContinuum::readInwardFluxes(inflname);
		}


		if(App::CalcGrainTemp) 
		{
			printf("\nreading temps..., %s\n", intempfile);
			GrainTemp::readTemps(intempfile, i);
		}

		
		/*for(int i=0;i<CContinuum::nbands;i++)
		{
			printf("%s %lf - %lf\n", CContinuum::bands[i].label, CContinuum::bands[i].leftCont, CContinuum::bands[i].rightCont);
		}*/
		
		
	}

	

	if(App::CalcIRBands && CContinuum::nbands <= 0)
	{
		CContinuum::readBands();
		CContinuum::assignBands();
	}

	FILE *TEST_CONT;

	TEST_CONT = fopen("test_cont_read.dat","w+");

	for(int il=0;il<CGeometry::nRadiuses[0];il++)
	{
		for(int i=0;i<CContinuum::cellCount;i++)
		{
			if(fabs(CContinuum::anu[i]-1.0)<0.005)
			{
				fprintf(TEST_CONT, "%le -> %le\n", CContinuum::anu[i], CContinuum::emits[0][il][i]);
			}
		}
	}

	fclose(TEST_CONT);
	
	if(App::AppMode == 0)
	{

		App::init();

		//printvec(obsray.angle);
		printf("Calculate rays");
	
		//exit(1);
		CDiffRay::runDiffRay();
	}
	else
	{
		for(int i=0;i<App::nApps;i++)
		{
			App::initAperture(i);
			CDiffRay::runDiffRay();	
		}	
		
	}


	//printvec(App::rayToObj.start);
	
	//CSolver::getPoints();

	/*vector3 *iPoints;
	iPoints = new vector3[2];
	iPoints = intersectSphere(App::rayToObj,6.95e+21);

	for (int i = 0; i < 2; ++i)
	{
		printvec(iPoints[i]);
	}*/
}
