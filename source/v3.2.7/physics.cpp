#include "basics.h"
#include "reader.h"
#include "physics.h"
#include "geometry.h"
#include "string.h"


int Physics::nsectors = 0;
int Physics::iSector = 0;
int Physics::nrows = 0;
int Physics::readen = 0;
double **Physics::Te;
double **Physics::Ne;

int Physics::initSectorials(int nsectors)
{
	//assuming we have only 1 sector (SS)
	Physics::Te = new double*[nsectors];
	Physics::Ne = new double*[nsectors];
	Physics::nsectors = nsectors;
}

int Physics::readPhysics(char fname[500],int iSector)
{

	if(Physics::nsectors<=0)
	{
		Physics::initSectorials(CGeometry::nSectors);
	}
	FILE *FREAD;

	FREAD = fopen(fname,"r");

	if(!FREAD || FREAD == NULL)
	{
		printf("%s\n", fname);
		CBasics::throwError("ERROR OPENING FILE WITH NULL POINTER!");
	}

	Physics::iSector = iSector;
	printf("Reading physics\n");
	FREAD = fopen(fname,"r");
	CReader::read_file_array(FREAD,2550,1,0,Physics::getRadii,Physics::putPhysics,Physics::putPhysicsVal);
	
	fclose(FREAD);

	return 1;
}

int Physics::putPhysics(int nlines)
{
	
	int is = Physics::iSector;
	{
		Physics::Te[is] = new double [Physics::nrows];
		Physics::Ne[is] = new double [Physics::nrows];

		for(int ir=0;ir<Physics::nrows;ir++)
		{
			Physics::Te[is][ir] = 0.0;
			Physics::Ne[is][ir] = 0.0;
		}
	}
}

int Physics::getRadii(int nrows)
{
	Physics::nrows = nrows;
}

int Physics::putPhysicsVal(int raw, int col, double val)
{
	
	if(col == 1) // not radii
	{
		Physics::Te[Physics::iSector][raw] = pow(10.0,val);
	}
	else
	if(col == 4) 
	{
		Physics::Ne[Physics::iSector][raw] = pow(10.0,val);
	}
	//printf("%d %d %d - %lf\n",CLine::iSector,raw,col,CLine::emits[CLine::iSector][raw][col-1] );
	
	return 1;
}
