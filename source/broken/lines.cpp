#include "basics.h"
#include "reader.h"
#include "lines.h"
#include "geometry.h"
#include "string.h"


int CLine::linesCount;
int CLine::nrows;
double ***CLine::emits;
int CLine::nsectors;
int CLine::iSector;

char CLine::linesCapDB[100000][14];;
double *CLine::linesEn;
int CLine::ndatabase;
int CLine::readen;
int CLine::lineIds[100];

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
	printf("Reading lines\n");
	FREAD = fopen(fname,"r");
	CReader::read_file_array(FREAD,2550,2,0,CLine::getRadii,CLine::putLine,CLine::putLineVal, CLine::getLineList);
	printf("here");

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

int CLine::getLineList(int irow,char str[1000])
{
	//printf("%d; %d\n", irow, CLine::readen);
	int lineIds[100];
	int linePos[100];

	int idSorted[100];
	printf("-----\n");
	if(irow==1 && CLine::readen==0)
	{
		CLine::readen = 1;
		char * pch;
		int il = 0;
		for(int i=0;i<CLine::ndatabase;i++)
		{
			
			pch=strcasestr(str,CLine::linesCapDB[i]);
			
			int pos = pch - str + 1;

			
			if(pos>4)
			{
				printf("%s  ->  %d\n", CLine::linesCapDB[i], pos);
				bool found = false;
				for(int j=0;j<il;j++)
				{
					if(linePos[j] == pos)
					{
						found = true;
					}
				}
				if(!found)
				{	
					linePos[il] = pos;
					lineIds[il] = i;
					il++;
				}
			}

			
			
		}

		for(int i=0;i<il;i++)
		{
			for(int j=0;j<il-i-1;j++)
			{
				if(linePos[j]>linePos[j+1])
				{
					int pos_next = linePos[j+1];
					int id_next  = lineIds[j+1];
					linePos[j+1] = linePos[j];
					lineIds[j+1] = lineIds[j];
					linePos[j] = pos_next;
					lineIds[j] = id_next;
				}
			}
		}

		for(int i=0;i<il;i++)
		{
			CLine::lineIds[i] = lineIds[i];
			printf("%d --- %d (%d)\n", lineIds[i],linePos[i],il);
		}

		printf("TTTTT\n");

	}

	return 1;
}

int CLine::readDatabase()
{
	FILE *LDB;

	LDB = fopen("data/database/line_ch.data","r");
	if(LDB == NULL)
	{
		printf("cant open\n");
		exit(1);
	}

	FILE *LDB2;
	
	LDB2 = fopen("data/database/line_en.data","r");
	if(LDB2 == NULL)
	{
		printf("cant open\n");
		exit(1);
	}
	
	int nCount = CReader::calc_lines_number(LDB,15);
	int nC2    = CReader::calc_lines_number(LDB2,15);

	CLine::ndatabase = nCount;

	CLine::linesEn    = new double[nCount+1];

	rewind(LDB);
	rewind(LDB2);

	char line[20];

	int ic=0;

	printf("read database...");
	while(!feof(LDB))
	{
		fgets(line,20,LDB);
		
		double el = 0;

		if(!feof(LDB2))
		{
			fscanf(LDB2,"%le",&el);	
		}
		
		if(el>1.e-35)
		{
			el = 912.0/(el);	
		}
		else
		{
			el = 1.e+8;
		}

		

		if(ic<nCount)
		{
			CLine::linesEn[ic] = el;
			//CLine::linesCapDB[ic] = new char[14];
			if(ic<=10)
			{
				printf("%s\n", line);
			}

			char subbuff[12];
			memcpy( subbuff, &line[1], 11 );
			subbuff[11] = '\0';

			sprintf(CLine::linesCapDB[ic],"%s",subbuff);

			if(ic<=10)
			{
				printf("TO: %s\n", CLine::linesCapDB[ic]);
			}
		}
		ic++;
	}
	printf("%d/%d\n",ic,nCount);
	//exit(1);
	fclose(LDB);
	fclose(LDB2);
	printf("Start\n");

	FILE *ft;

	ft = fopen("fop.dat","w+");
	for(int i=0;i<nCount;i++)
	{
		fprintf(ft,"%s - %le\n", CLine::linesCapDB[i], CLine::linesEn[i]);
	}
	fclose(ft);
	//exit(1);
}