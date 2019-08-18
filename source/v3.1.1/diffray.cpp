#include "basics.h"
#include "integration.h"
#include "app.h"
#include "continuum.h"
#include "lines.h"
#include "iterator.h"
#include "diffray.h"

int CDiffRay::runDiffRay()
{
	int mode = 1;
	CIteration::initIteration(mode);

	double dltLines = 1.0;
	double dltCont  = 1.0;

	while((dltCont>0.02 || dltLines>0.02) && CIteration::nIteration < 10)
	{
		CIteration::doIteration();

		if(App::CalcLines)
		{
			dltLines = CIteration::getdLine();
		}
		else
		{
			dltLines = 0.01;
		}

		if(App::CalcCont)
		{
			dltCont  = CIteration::getdCont();	
		}
		else
		{
			dltCont = 0.01;
		}

		printf("DELTA[%d]: %le;%le\n", CIteration::nIteration, dltLines, dltCont);

		if(App::CalcLines)
		{
			char fname[255];
			sprintf(fname,"Iter_%d_lines_%d.dat",CIteration::nIteration,mode);
			FILE *FP;
			FP = fopen(fname,"w+");
			
			for(int il=0;il<CLine::linesCount;il++)
			{
				fprintf(FP, "%s\t", CLine::linesCapDB[CLine::lineIds[il]]);
			}
			
			fprintf(FP, "\n");

			for(int il=0;il<CLine::linesCount;il++)
			{
				fprintf(FP, "%lf\t", log10(CIntegration::lumfactor*CIntegration::flux_lines[il]*pow(App::distance,2.0)));
			}
			
			fclose(FP);
		}

		if(App::CalcCont)
		{
			char fname[255];
			sprintf(fname,"%s/contarr.txt",App::output_dir);
			FILE *FP;

			printf("%s\n", fname);
			FP = fopen(fname,"w+");

			fprintf(FP, "En, Ryd\t\tDiffuse\t\tAttenuated direct\n");
			
			for(int il=0;il<CContinuum::cellCount;il++)
			{
				fprintf(FP, "%le\t\t%le\t\t%le\n", CContinuum::anu[il], CIntegration::flux_continua[il],0.0);
			}

			fclose(FP);	
		}

		printf("Dist: %le\n", App::distance);

		CIteration::nIteration++;	
	}
}