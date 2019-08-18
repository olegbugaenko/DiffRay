#include "test.h"
#include "reader.h"
#include "lines.h"
#include "geometry.h"
#include "app.h";
#include "continuum.h";
#include "geom3d.h";

int main()
{
	FILE *F;
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
	printf("test: %d\n",res);

	F = fopen("data/commands.ini","r");
	char **commands;
	
	commands = CReader::read_commands(F);
	while(commands != NULL)
	{
		for(int icm=0;icm<10;icm++)
		{
			printf("'%s' ",commands[icm]);
		}
		printf("\n");
		commands = CReader::read_commands(F);
	};
	
	fclose(F);

	//read geometry

	char path[500];

	sprintf(path,"data/test/Emis_Cont_SectorNo1_Age60.00Myr.dat");

	for(int i=0;i<App::nSectors;i++)
	{
		CGeometry::readRadiuses(path,i);
	}

	for(int i=0;i<CGeometry::nSectors;i++)
	{
		printf("nSectors[%d] = %d (%le - %le)\n", i, CGeometry::nRadiuses[i],CGeometry::outer_radius[i][0],CGeometry::outer_radius[i][CGeometry::nRadiuses[i]]);
	}

	//read lines
	printf("READING LINES\n");

	double **resarr;
	char flname[500];
	char opflname[500];
	char trflname[500];

	sprintf(flname, "data/test/Emis_Lines_SectorNo1_Age60.00Myr.dat");

	//CLine::init();

	//double ***CLine::emits;

	for(int i=0;i<CGeometry::nSectors;i++)
	{
		CLine::readLines(flname,i);
		//printf("Emts %d = %lf\n", i, CLine::emits[i][113][0]);
		//printf("Prev %d = %lf\n", i, CLine::emits[0][113][0]);
		
	}
	printf("READEN: %le - %le\n",CLine::emits[1][113][0],CLine::emits[1][83][24]);

	sprintf(flname, "data/test/Continuum1_SB99_Age60.00Myr.dat");
	CContinuum::readMesh(flname);

	printf("Mesh count: %d\n", CContinuum::cellCount);
	printf("%le - %le - %le\n", CContinuum::anu[0], CContinuum::anu[1089], CContinuum::anu[CContinuum::cellCount-1]);

	/*sprintf(flname, "data/test/Emis_Cont_SectorNo1_Age60.00Myr.dat");
	sprintf(opflname, "data/test/Opac_Cont_SectorNo1_Age60.00Myr.dat");
	sprintf(trflname, "data/test/Emis_Trans_SectorNo1_Age60.00Myr.dat");
	for(int i=0;i<CGeometry::nSectors;i++)
	{
		printf("Reading cont %d\n", i);
		CContinuum::readEmissivity(flname,i);
		printf("reading opacity...\n");
		CContinuum::readOpacity(opflname,i);
		printf("reading trans...\n");
		CContinuum::readTransitions(trflname,i);
	}
	//printf("Mesh count: %d\n", CContinuum::cellCount);
	printf("%le - %le - %le\n", CContinuum::emits[0][3][5], CContinuum::emits[0][5][3], CContinuum::emits[0][3][CContinuum::cellCount-1]);
	printf("Opacs: %le - %le - %le\n", CContinuum::opacs[0][3][5], CContinuum::opacs[0][5][3], CContinuum::opacs[0][3][CContinuum::cellCount-1]);
	printf("Trans: %le - %le - %le\n", CContinuum::trans[0][3][5], CContinuum::trans[0][5][3], CContinuum::trans[0][3][CContinuum::cellCount-1]);
*/
	vector3 start = vector3(3.0,0.0,0.0);
	vector3 angle = vector3(0.0,7*M_PI/4,M_PI/4);

	TRay radray = createRay(start,angle);

	/*vector3* intersections;
	intersections = new vector3[2];

	intersections = intersectSphere(radray,25.0);

	for(int isv=0;isv<2;isv++)
	{
		printf("VEC%d: ",isv);
		printvec(intersections[isv]);
	}*/

	vector3 x = vector3(-1.0,3.0,5.0);
	vector3 y = vector3(-3.0,2.0,5.0);
	vector3 z = vector3(3.0,0.0,2.0);

	TPlane result = getFromPoints(x,y,z);
	printplane(result);
	printf("------------PLANE---------------\n");
	vector3 axis = vector3(1.0,0.0,0.0);
	axis = rotate(axis,vector3(0.0,M_PI/4.0,M_PI/4.0));
	printvec(axis);
	TPlane sector = getFromPhi(axis,M_PI);
	printplane(sector);
}