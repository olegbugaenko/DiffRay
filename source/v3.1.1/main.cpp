#include "test.h"
#include "reader.h"
#include "lines.h"
#include "geometry.h"
#include "app.h";
#include "continuum.h";
#include "geom3d.h";
#include "solver.h";
#include "diffray.h" 

int main()
{
	/*FILE *F;
	F = fopen("data/readme.txt","r");
	char *mychar;
	mychar = new char[500];
	while(CReader::readLine(F,mychar,CReader::callable) != 0)
	{
		printf("\n ---- \n%s\n",mychar);
	};
	fclose(F);
	int **myarr;
	myarr = CBasics::int2Arr(myarr,8,10,1);
	myarr[4][1] = -8;
	printf("%d\n",myarr[4][1]);
	printf("%d\n",myarr[1][1]);

	for(int i=0;i<8;i++)
	{
		printf("|");
		for(int j=0; j<10; j++)
		{
			printf("%d |",myarr[i][j]);
		}
		printf("\n");
	}
	
	int inval = 34;
	int res = CTest::test(inval);
	printf("test: %d\n",res);*/

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
	char flname[500];
	char opflname[500];
	char trflname[500];

	

	printf("dbentry\n");
	CLine::readDatabase();

	//CLine::init();

	//double ***CLine::emits;

	if(App::CalcLines)
	{
		for(int i=0;i<CGeometry::nSectors;i++)
		{
			sprintf(flname, "%s/Emis_Lines_SectorNo%d_Age%.2lfMyr.dat",App::input_dir,i+1,App::age);
			CLine::readLines(flname,i);
			//printf("Emts %d = %lf\n", i, CLine::emits[i][113][0]);
			//printf("Prev %d = %lf\n", i, CLine::emits[0][113][0]);
			
		}
		printf("READEN: %le - %le\n",CLine::emits[1][113][0],CLine::emits[1][83][24]);
	}
	

	sprintf(flname, "%s/Continuum%d_SB99_Age%.2lfMyr.dat",App::input_dir,1,App::age);
	CContinuum::readMesh(flname);

	if(CContinuum::cellCount<10)
	{
		sprintf(flname, "%s/Continuum%d_Last_Age%.2lfMyr.dat",App::input_dir,1,App::age);
		CContinuum::readMesh(flname);

	}

	printf("Mesh count: %d\n", CContinuum::cellCount);
	printf("%le - %le - %le\n", CContinuum::anu[0], CContinuum::anu[1089], CContinuum::anu[CContinuum::cellCount-1]);


	for(int i=0;i<CGeometry::nSectors;i++)
	{
		sprintf(flname, "%s/Emis_Cont_SectorNo%d_Age%.2lfMyr.dat",App::input_dir,i+1,App::age);
		sprintf(opflname, "%s/Opac_Cont_SectorNo%d_Age%.2lfMyr.dat",App::input_dir,i+1,App::age);
		sprintf(trflname, "%s/Emis_Trans_SectorNo%d_Age%.2lfMyr.dat",App::input_dir,i+1,App::age);

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
		
		if(App::CalcLines)
		{
			printf("reading trans...\n");
			CContinuum::readTransitions(trflname,i);
		}
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

	
	App::init();

	//printvec(obsray.angle);
	printf("Calculate rays");
	
	//exit(1);
	CDiffRay::runDiffRay();
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
