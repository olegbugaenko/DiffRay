#include "basics.h"
#include "integration.h"
#include "app.h"
#include "continuum.h"
#include "lines.h"
#include "iterator.h"
#include "diffray.h"
#include "geometry.h"
#include "abund.h"
#include "grain_temp.h"

int compareRadii (const void * pa, const void * pb)
{
  const double *a = (double *) pa;
  const double *b = (double *) pb;
  if(a[0] > b[0])
  	return -1;
  
  if(a[0] < b[0])
  	return 1;
  
  return 0;
}

int CDiffRay::runDiffRay()
{
	int mode = App::integration_mode;
	CIteration::initIteration(mode);

	double dltLines = 1.0;
	double dltCont  = 1.0;
	App::isStatMode = false;



	while((dltCont > App::precision || dltLines > App::precision) && CIteration::nIteration < App::maxiterations)
	{
		Abund::refreshStatistics();
		CLine::refreshStatistics();
		GrainTemp::refreshStatistics();
		
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

		if(App::AppMode != 0)
		{
			CIntegration::lumfactor *= (1/(4*M_PI));
		}


		if(App::CalcCont)
		{
			char fname[255];
			sprintf(fname,"%s/contarr.txt",App::output_dir);
			FILE *FP;

			printf("Printing to %s\n", fname);
			FP = fopen(fname,"w+");

			fprintf(FP, "En, Ryd\t\tDiffuse\t\tAttenuated direct\t\tOutLine\n");
			
			for(int il=0;il<CContinuum::cellCount;il++)
			{
				fprintf(FP, "%le\t\t%le\t\t%le\t\t%le\n", CContinuum::anu[il], CIntegration::flux_continua[il],CIntegration::flux_intrinsic[il]*App::cov_fac,CIntegration::flux_transitions[il]);
			}

			fclose(FP);	

			FILE *FP2;

			char fname2[255];
			sprintf(fname2,"%s/spectra.dat",App::output_dir);
			
			printf("%s\n", fname2);
			FP2 = fopen(fname2,"w+");

			fprintf(FP, "En, Ryd\t\tAttenuated direct\t\tDiffuse\t\tOutLine\n");
			
			for(int il=1;il<CContinuum::cellCount;il++)
			{
				//double cn = 1.8e-11*pow(CContinuum::anu[il],2.0)/(CContinuum::anu[il]- CContinuum::anu[il-1]);
				//double cn = pow(App::distance/CGeometry::outer_radius[0][CGeometry::nRadiuses[0]],2.0);//*CIntegration::lumfactor;
				double cn = 1.0;
				fprintf(FP2, "%lf\t\t%le\t\t%le\t\t%le\n", 912.0/CContinuum::anu[il],CIntegration::flux_intrinsic[il]*cn, CIntegration::flux_continua[il]*cn,CIntegration::flux_transitions[il]*cn*App::cov_fac);
				//printf("Printing\n");
			}

			fclose(FP2);	
		}

		if(App::CalcLines)
		{
			char fname[255];
			sprintf(fname,"%s/linesarr.txt",App::output_dir);
			FILE *FP;

			printf("%s\n", fname);
			FP = fopen(fname,"w+");

			fprintf(FP, "Line\t\tFlux at observer\t\tLuminosity\n");
			
			for(int il=0;il<CLine::linesCount;il++)
			{
				fprintf(FP, "%s\t\t%le\t\t%lf\t\t%lf\n", CLine::linesCapDB[CLine::lineIds[il]], CIntegration::lumfactor*CIntegration::flux_lines[il],log10(1.e-50 + 4*M_PI*CIntegration::lumfactor*CIntegration::flux_lines[il]*pow(App::distance,2.0)),CIntegration::flux_lines[il]/(CIntegration::flux_lines[0]));
			}

			fclose(FP);	
		}

		if(App::CalcLines)
		{
			char fname[255];
			sprintf(fname,"%s/emergent.txt",App::output_dir);
			FILE *FP;

			printf("%s\n", fname);
			FP = fopen(fname,"w+");

			fprintf(FP, "Line\t\tFlux at observer\t\tLuminosity\n");
			
			for(int il=0;il<CLine::linesCount;il++)
			{
				fprintf(FP, "%s\t\t%le\t\t%lf\t\t%lf\n", 
					CLine::linesCapDB[CLine::lineIds[il]], 
					CIntegration::lumfactor*CLine::cumulative[il],
					log10(1.e-50 + 4*M_PI*CIntegration::lumfactor*CLine::cumulative[il]),
					CLine::cumulative[il]/(CLine::cumulative[0]+1.e-50)
				);
			}

			fclose(FP);	
		}

		if(App::CalcLines)
		{
			char fname[255];
			sprintf(fname,"%s/lifr.txt",App::output_dir);
			FILE *FP;

			printf("%s\n", fname);
			FP = fopen(fname,"w+");

			fprintf(FP, "#");
			
			for(int il=0;il<CLine::linesCount;il++)
			{
				fprintf(FP, "%s\t", 
					CLine::linesCapDB[CLine::lineIds[il]]
				);
			}

			fprintf(FP, "\n");
			
			for(int il=0;il<CLine::linesCount;il++)
			{
				fprintf(FP, "%le\t", 
					CLine::cumulative[il]
				);
			}

			fclose(FP);	
		}

		if(App::CalcLines)
		{
			char fname[255];
			sprintf(fname,"%s/lifr2.txt",App::output_dir);
			FILE *FP;

			printf("%s\n", fname);
			FP = fopen(fname,"w+");

			fprintf(FP, "#");
			
			for(int il=0;il<CLine::linesCount;il++)
			{
				fprintf(FP, "%s\t", 
					CLine::linesCapDB[CLine::lineIds[il]]
				);
			}

			fprintf(FP, "\n");
			
			for(int il=0;il<CLine::linesCount;il++)
			{
				fprintf(FP, "%le\t", 
					CIntegration::flux_lines[il]
				);
			}

			fclose(FP);	
		}

			
		if(App::CalcLines)
		{
			char fname[255];
			sprintf(fname,"%s/abmass.txt",App::output_dir);
			FILE *FP;

			printf("%s\n", fname);
			FP = fopen(fname,"w+");

			for(int i=0;i<Abund::nElements;i++){
				fprintf(FP, "%s\t\t",Abund::elements[i]);
			}

			fprintf(FP, "\n");

			
			for(int i=0;i<Abund::nElements;i++){
				fprintf(FP, "%lf\t\t",log10(Abund::abmass[i]/Abund::mass));
			}

			fclose(FP);	

			sprintf(fname,"%s/abemso.txt",App::output_dir);
			
			printf("%s\n", fname);
			FP = fopen(fname,"w+");

			for(int i=0;i<Abund::nElements;i++){
				fprintf(FP, "%s\t\t",Abund::elements[i]);
			}

			fprintf(FP, "\n");

			
			for(int i=0;i<Abund::nElements;i++){
				fprintf(FP, "%lf\t\t",log10(Abund::abemso[i]/Abund::emso));
			}

			fclose(FP);	
		}
		

		if(CContinuum::nbands > 0)
		{

			char fname[255];
			sprintf(fname,"%s/bands.txt",App::output_dir);
			FILE *FP;
			
			//printf("%s\n", fname);
			char strnm1[255];
			char strnm2[255];
			sprintf(strnm1,"mips24_app%d.dat",App::iApp);
			sprintf(strnm2,"mips70_app%d.dat",App::iApp);
			
			FP = fopen(fname,"w+");
			fprintf(FP, "Band\t\tFlux at observer\t\tFlux, Yan\t\tLuminosity\n");

			for(int i = 0; i < CContinuum::nbands; i++)
			{
				IRBand cBand = CContinuum::bands[i];

				double bandInten = 0.;
				double Yan = 0.0;

				int iT = cBand.nAnus-1;
				double trans = 1.0;

				for(int iCl = cBand.iNuLeft; iCl < cBand.iNuRight; iCl++)
				{
					double cEn = CContinuum::anu[iCl];
					//trans = 1.0;
					while(iT >= 0 && cEn > 912./cBand.ennu[iT])
					{
						trans = cBand.transmitions[iT];
						iT--;
					}
					
					//printf("Trans[%d]: %le; %le > %le\n", iT, trans, cEn, 912.0/cBand.ennu[iT]);
					bandInten += trans*((CIntegration::flux_continua[iCl] + CIntegration::flux_transitions[iCl])*App::cov_fac + CIntegration::flux_intrinsic[iCl])*2.179874099E-11*CContinuum::anu[iCl];
					Yan += trans*((CIntegration::flux_continua[iCl] + CIntegration::flux_transitions[iCl])*App::cov_fac + CIntegration::flux_intrinsic[iCl])*2.179874099E-11/(1.e-23*3.e+8/(912*1.e-10));
					
				}


				bandInten *= CIntegration::lumfactor;
				Yan *= CIntegration::lumfactor;

				fprintf(FP, "%s\t\t%le\t\t%le\t\t%lf\t\t%lf\n", cBand.label, bandInten, Yan, log10(1.e-50 + 4*M_PI*bandInten*pow(App::distance,2.0)),bandInten/(CIntegration::flux_lines[0]*CIntegration::lumfactor));
			}

			fclose(FP);	
		}

		CIteration::nIteration++;	
	}

	
	printf("LCP: %d\n",CLine::iStat);


	CIteration::obtainStatistics();


	/* Halt! We should sort statistics first! 
	   We can use qsort for this
	*/
	printf("OS\n");
	qsort (Abund::statistics, 500000, sizeof Abund::statistics[0], compareRadii);
	printf("SORTEDA\n");
	qsort (CLine::statistics, 500000, sizeof CLine::statistics[0], compareRadii);
	printf("SORTEDL\n");
	qsort (GrainTemp::statistics, 500000, sizeof GrainTemp::statistics[0], compareRadii);
	printf("SORTEDL\n");
	int i = 0;

	FILE *AS;

	char fnamedas[255];
	sprintf(fnamedas,"%s/gt_distribution_saved.txt",App::output_dir);

	AS = fopen(fnamedas, "w+");

	for(int i=0;i<500000;i++) 
	{
		if(GrainTemp::statistics[i][0] > 1.e-36) {
			fprintf(AS, "%le\t%le\n", GrainTemp::statistics[i][0], GrainTemp::statistics[i][1]);
		}		
	}

	fclose(AS);


	/*** Smoothing and averaging ***/
	for(int ir=0; ir<1; ir++)
	{
		for(int is=1;is<CLine::iStat-1; is++)
		{
			if(Abund::statistics[is][1] < (Abund::statistics[is-1][1] + Abund::statistics[is+1][1])/3)
			{
				for(int k=1;k<=Abund::nElements;k++)
				{
					Abund::statistics[is][k] = (Abund::statistics[is-1][k] + Abund::statistics[is][k] + Abund::statistics[is+1][k])/3;
				}
			}
		}
	}

	for(int ir=0; ir<1; ir++)
	{
		for(int is=1;is<CLine::iStat-1; is++)
		{
			if(CLine::statistics[is][1] < (CLine::statistics[is-1][1] + CLine::statistics[is+1][1])/3)
			{
				for(int k=1;k<=CLine::linesCount;k++)
				{
					CLine::statistics[is][k] = (CLine::statistics[is-1][k] + CLine::statistics[is][k] + CLine::statistics[is+1][k])/3;
				}
			}
		}
	}

	for(int ir=0; ir<1; ir++)
	{
		for(int is=1;is<CLine::iStat-1; is++)
		{
			if(GrainTemp::statistics[is][1] < (GrainTemp::statistics[is-1][1] + GrainTemp::statistics[is+1][1])/3)
			{
				for(int k=1;k<=GrainTemp::nBins+1;k++)
				{
					GrainTemp::statistics[is][k] = (GrainTemp::statistics[is-1][k] + GrainTemp::statistics[is][k] + GrainTemp::statistics[is+1][k])/3;
				}
			}
		}
	}

	int nLines = CLine::iStat;
	int nAbunds = CLine::iStat; 
	int nTemps = CLine::iStat; 

	int ik = 0;
	int ic = 0;
	int ikprev = 0;
	double av[100];
	for(int k=0;k<100;k++)
	{
		av[k] = 0.;
	}

	while(ic < CLine::iStat)
	{
		double dR = 0.01;

		for(int k=0;k<100;k++)
		{
			av[k] = 1.e-42;
		}
		int nElems = 0;

		while(dR < 1.e+19 && ik<CLine::iStat)
		{
			ik++;
			dR = fabs(Abund::statistics[ik][0] - Abund::statistics[ikprev][0]);
			
			for(int k=1;k<Abund::nElements;k++)
			{
				av[k] += Abund::statistics[ik][k];
			}
			nElems++;
			
		}

		ikprev = ik;

		//printf("av: %le; n=%d; LC=%d\n",av[1], nElems, CLine::iStat);

		for(int k=1;k<Abund::nElements;k++)
		{
			Abund::statistics[ic][k] = av[k]/nElems;
		}

		
		

		if(ik >= CLine::iStat)
		{
			Abund::statistics[ic][0] = Abund::statistics[CLine::iStat-1][0];
			nAbunds = ic+1;
			break;
		}
		else
		{
			Abund::statistics[ic][0] = Abund::statistics[ik][0];

		}

		ic++;

	}
	
	printf("Abunds selected: %d\n",nAbunds);

	char fnameda[255];
	sprintf(fnameda,"%s/abund_distribution.txt",App::output_dir);
	
	FILE *AbundF;

	AbundF = fopen(fnameda, "w+");
	while(i<nAbunds && i<500000)
	{
		if(Abund::statistics[i][1] < 1.e-36 || CLine::statistics[i][0] >= 1.e+25)
		{
			i++;
			continue;
		}
		int j = 0;
		while(j <= Abund::nElements && j < 100)
		{
			fprintf(AbundF, "%le\t", Abund::statistics[i][j]);
			j++;
		}
		fprintf(AbundF, "\n");
		i++;
	}
	fclose(AbundF);


	ik = 0;
	ic = 0;
	ikprev = 0;
	
	for(int k=0;k<100;k++)
	{
		av[k] = 0.;
	}

	while(ic < CLine::iStat)
	{
		double dR = 0.01;

		for(int k=0;k<100;k++)
		{
			av[k] = 1.e-42;
		}
		int nElems = 0;

		while(dR < 1.e+19 && ik<CLine::iStat)
		{
			ik++;
			dR = fabs(CLine::statistics[ik][0] - CLine::statistics[ikprev][0]);
			
			for(int k=1;k<CLine::linesCount;k++)
			{
				av[k] += CLine::statistics[ik][k];
			}
			nElems++;
			
		}

		ikprev = ik;

		//printf("av: %le; n=%d; LC=%d\n",av[1], nElems, CLine::iStat);

		for(int k=1;k<CLine::linesCount;k++)
		{
			CLine::statistics[ic][k] = av[k]/nElems;
		}

		
		

		if(ik >= CLine::iStat)
		{
			CLine::statistics[ic][0] = CLine::statistics[CLine::iStat-1][0];
			nLines = ic+1;
			break;
		}
		else
		{
			CLine::statistics[ic][0] = CLine::statistics[ik][0];

		}

		ic++;

	}
	
	printf("Lines selected: %d\n",nLines);


	i = 0;

	char fnameel[255];
	sprintf(fnameel,"%s/lineems_distribution.txt",App::output_dir);
	
	FILE *LineEmF;
	printf("LC: %d\n",CLine::iStat);
	LineEmF = fopen(fnameel, "w+");
	while(i<nLines && i<500000)
	{
		if(CLine::statistics[i][0] > 1.e-35 && CLine::statistics[i][0] < 1.e+25)
		{
			int j = 0;
			while(j <= CLine::linesCount && j < 100)
			{
				fprintf(LineEmF, "%le\t", CLine::statistics[i][j]);
				j++;
			}
			fprintf(LineEmF, "\n");
		}
		
		i++;
	}
	fclose(LineEmF);





	ik = 0;
	ic = 0;
	ikprev = 0;
	
	for(int k=0;k<100;k++)
	{
		av[k] = 0.;
	}

	while(ic < CLine::iStat)
	{
		double dR = 0.01;

		for(int k=0;k<100;k++)
		{
			av[k] = 1.e-42;
		}
		int nElems = 0;

		while(dR < 1.e+19 && ik<CLine::iStat)
		{
			ik++;
			dR = fabs(GrainTemp::statistics[ik][0] - GrainTemp::statistics[ikprev][0]);
			
			for(int k=1;k<GrainTemp::nBins+1;k++)
			{
				av[k] += GrainTemp::statistics[ik][k];
			}
			nElems++;
			
		}

		ikprev = ik;

		//printf("av: %le; n=%d; LC=%d\n",av[1], nElems, CLine::iStat);

		for(int k=1;k<GrainTemp::nBins+1;k++)
		{
			GrainTemp::statistics[ic][k] = av[k]/nElems;
		}

		
		

		if(ik >= CLine::iStat)
		{
			GrainTemp::statistics[ic][0] = GrainTemp::statistics[CLine::iStat-1][0];
			nTemps = ic+1;
			break;
		}
		else
		{
			GrainTemp::statistics[ic][0] = GrainTemp::statistics[ik][0];

		}

		ic++;

	}
	
	printf("Temps selected: %d\n",nTemps);


	i = 0;

	char fnametmp[255];
	sprintf(fnametmp,"%s/graintemp_distribution.txt",App::output_dir);
	
	FILE *GrainTempF;
	printf("TC: %d\n",CLine::iStat);
	GrainTempF = fopen(fnametmp, "w+");
	while(i<nTemps && i<500000)
	{
		if(GrainTemp::statistics[i][0] > 1.e-35 && GrainTemp::statistics[i][0] < 1.e+25)
		{
			int j = 0;
			while(j <= GrainTemp::nBins+1 && j < 100)
			{
				fprintf(GrainTempF, "%le\t", GrainTemp::statistics[i][j]);
				j++;
			}
			fprintf(GrainTempF, "\n");
		}
		
		i++;
	}
	fclose(GrainTempF);

}