#include "basics.h"
#include "reader.h"
#include "continuum.h"
#include "geometry.h"


int CContinuum::cellCount = 0;
int CContinuum::nrows;
double ***CContinuum::emits;
int CContinuum::nsectors = 0;
int CContinuum::iSector;
double ***CContinuum::opacs;
double ***CContinuum::trans;
double *CContinuum::anu;
double CContinuum::endSkip;
int CContinuum::skipInterval;


int CContinuum::initSectorials(int nsectors)
{
	//assuming we have only 1 sector (SS)
	printf("initSectorials\n");
	CContinuum::emits = new double**[nsectors];
	CContinuum::opacs = new double**[nsectors];
	CContinuum::trans = new double**[nsectors];
	CContinuum::nsectors = nsectors;

}

int CContinuum::readMesh(char fname[500])
{
	if(CContinuum::nsectors<=0)
	{
		CContinuum::initSectorials(CGeometry::nSectors);
	}
	FILE *FREAD;

	FREAD = fopen(fname,"r");

	if(!FREAD || FREAD == NULL)
	{
		printf("%s\n", fname);
		return 0;
		//CBasics::throwError("ERROR OPENING FILE WITH NULL POINTER!");
	}

	CReader::read_file_array(FREAD,255,1,0,CContinuum::setMeshSize,0,CContinuum::putMeshEnergy);
	fclose(FREAD);

	return 1;
}

int CContinuum::setMeshSize(int nrows)
{
	if(CContinuum::cellCount<=0)
	{
		CContinuum::cellCount = nrows;
		CContinuum::anu = new double[nrows];

		for(int i=0;i<CGeometry::nSectors;++i)
		{
			int nrads = CGeometry::nRadiuses[i]+1;
			CContinuum::emits[i] = new double*[nrads];
			CContinuum::opacs[i] = new double*[nrads];
			CContinuum::trans[i] = new double*[nrads];
			for(int ir=0;ir<CGeometry::nRadiuses[i];++ir)
			{
				CContinuum::emits[i][ir] = new double[CContinuum::cellCount];
				CContinuum::opacs[i][ir] = new double[CContinuum::cellCount];
				CContinuum::trans[i][ir] = new double[CContinuum::cellCount];
				for (int ic = 0; ic < CContinuum::cellCount; ++ic)
				{
					CContinuum::emits[i][ir][ic] = 0.00;
					CContinuum::opacs[i][ir][ic] = 0.00;
					CContinuum::trans[i][ir][ic] = 0.00;
				}
			}
		}
	}
}

int CContinuum::putMeshEnergy(int raw, int col, double val)
{
	
	if(col<=0) // not radii
	{
		CContinuum::anu[raw] = val;
	}
	//printf("%d %d %d - %lf\n",CLine::iSector,raw,col,CLine::emits[CLine::iSector][raw][col-1] );
	
	return -1;
}

int CContinuum::readEmissivity(char fname[500],int iSector)
{

	printf("FUNTION readEmissivity\n");
	if(CContinuum::nsectors<=0)
	{
		CContinuum::initSectorials(CGeometry::nSectors);
	}
	FILE *FREAD;

	FREAD = fopen(fname,"r");

	if(!FREAD || FREAD == NULL)
	{
		printf("%s\n", fname);
		printf("WARNING: FILE NOT FOUND");
		return 0;
	}

	CContinuum::iSector = iSector;
	printf("Reader\n");
	CReader::read_file_array(FREAD,255*1024,1,0,CContinuum::getRadii,0,CContinuum::putContEmissivity);
	printf("Readen\n");
	fclose(FREAD);

	return 1;
}

int CContinuum::getRadii(int nrows)
{
	CContinuum::nrows = nrows;

	if(nrows != CGeometry::nRadiuses[CContinuum::iSector] && nrows+1 != CGeometry::nRadiuses[CContinuum::iSector])
	{
		printf("Radiuses[%d] = %d, but %d continuums found!",CContinuum::iSector,CGeometry::nRadiuses[CContinuum::iSector],nrows);
		CBasics::throwError("ERROR! Wrong continuum data!");
	}
}

int CContinuum::putContEmissivity(int raw, int col, double val)
{
	
	if(col>0 && col-1<CContinuum::cellCount) // not radii
	{
		CContinuum::emits[CContinuum::iSector][raw][col-1] = val;
	}
	//printf("%d %d %d - %lf\n",CLine::iSector,raw,col,CLine::emits[CLine::iSector][raw][col-1] );
	
	return 1;
}

int CContinuum::readOpacity(char fname[500],int iSector)
{

	printf("CNS: %d\n", CContinuum::nsectors);
	if(CContinuum::nsectors<=0)
	{
		CContinuum::initSectorials(CGeometry::nSectors);
	}
	FILE *FREAD;

	FREAD = fopen(fname,"r");

	if(!FREAD || FREAD == NULL)
	{
		printf("%s\n", fname);
		printf("WARNING: FILE NOT FOUND");
		return 0;
	}

	CContinuum::iSector = iSector;

	CReader::read_file_array(FREAD,255*1024,1,0,CContinuum::getRadii,0,CContinuum::putContOpacity);
	fclose(FREAD);

	return 1;
}

int CContinuum::putContOpacity(int raw, int col, double val)
{

	if(col>0 && col-1<CContinuum::cellCount) // not radii
	{
		

		CContinuum::opacs[CContinuum::iSector][raw][col-1] = val;
	}
	
	return 1;
}

int CContinuum::readTransitions(char fname[500],int iSector)
{

	printf("CNS: %d\n", CContinuum::nsectors);
	if(CContinuum::nsectors<=0)
	{
		CContinuum::initSectorials(CGeometry::nSectors);
	}
	FILE *FREAD;

	FREAD = fopen(fname,"r");

	if(!FREAD || FREAD == NULL)
	{
		printf("%s\n", fname);
		printf("WARNING: FILE NOT FOUND");
		return 0;
	}

	CContinuum::iSector = iSector;

	CReader::read_file_array(FREAD,255*1024,1,0,CContinuum::getRadii,0,CContinuum::putTransitions);
	fclose(FREAD);

	return 1;
}

int CContinuum::putTransitions(int raw, int col, double val)
{

	if(col>0 && col-1<CContinuum::cellCount) // not radii
	{
		

		CContinuum::trans[CContinuum::iSector][raw][col-1] = val;
	}
	
	return 1;
}

/*int CContinuum::readLines(char fname[500],int iSector)
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
}*/