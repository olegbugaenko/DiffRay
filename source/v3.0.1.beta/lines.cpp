#include "basics.h"
#include "reader.h"
#include "lines.h"
#include "geometry.h"


int CLine::linesCount;
int CLine::nrows;
double ***CLine::emits;
int CLine::nsectors;
int CLine::iSector;


int CLine::initSectorials(int nsectors)
{
	//assuming we have only 1 sector (SS)
	CLine::emits = new double**[nsectors];
	CLine::nsectors = nsectors;
}

int CLine::readLines(char fname[500],int iSector)
{

	if(CLine::nsectors<=0)
	{
		CLine::initSectorials(CGeometry::nSectors);
	}
	FILE *FREAD;

	FREAD = fopen(fname,"r");

	if(!FREAD || FREAD == NULL)
	{
		printf("%s\n", fname);
		CBasics::throwError("ERROR OPENING FILE WITH NULL POINTER!");
	}

	CLine::iSector = iSector;

	FREAD = fopen(fname,"r");
	CReader::read_file_array(FREAD,2550,2,0,CLine::getRadii,CLine::putLine,CLine::putLineVal);
	fclose(FREAD);

	return 1;
}

int CLine::putLine(int nlines)
{
	CLine::linesCount = nlines-1;

	int is = CLine::iSector;
	{
		CLine::emits[is] = new double *[CLine::nrows];

		for(int ir=0;ir<CLine::nrows;ir++)
		{
			CLine::emits[is][ir] = new double[CLine::linesCount];

			for(int il=0;il<CLine::linesCount;il++)
			{
				CLine::emits[is][ir][il] = 0.0;
			}
		}
	}
}

int CLine::getRadii(int nrows)
{
	CLine::nrows = nrows;
}

int CLine::putLineVal(int raw, int col, double val)
{
	
	if(col>0) // not radii
	{
		CLine::emits[CLine::iSector][raw][col-1] = pow(10.0,val);
	}
	//printf("%d %d %d - %lf\n",CLine::iSector,raw,col,CLine::emits[CLine::iSector][raw][col-1] );
	
	return 1;
}