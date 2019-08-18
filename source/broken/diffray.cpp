#include "basics.h"
#include "integration.h"
#include "app.h"
#include "continuum.h"
#include "lines.h"
#include "iterator.h"
#include "diffray.h"
#include "geometry.h"
#include "abund.h"

int CDiffRay::runDiffRay()
{
	int mode = 0;
	CIteration::initIteration(mode);

	double dltLines = 1.0;
	double dltCont  = 1.0;
	App::isStatMode = false;

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

		if(App::AppMode != 0)
		{
			CIntegration::lumfactor *= (1/(4*M_PI));
		}

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
			FILE *FDMP;
			FILE *FDMP2;

			//printf("%s\n", fname);
			char strnm1[255];
			char strnm2[255];
			sprintf(strnm1,"mips24_app%d.dat",App::iApp);
			sprintf(strnm2,"mips70_app%d.dat",App::iApp);
			
			FP = fopen(fname,"w+");
			FDMP = fopen(strnm1,"w+");
			fprintf(FP, "Band\t\tFlux at observer\t\tFlux, Yan\t\tLuminosity\n");

			FDMP2 = fopen(strnm2,"w+");
			//fprintf(FP, "Band\t\tFlux at observer\t\tLuminosity\n");

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
					if(i == 1)
					{
						double factor = pow(App::distance/CGeometry::outer_radius[0][CGeometry::nRadiuses[0]],2.0)*CIntegration::lumfactor;
						fprintf(FDMP,"%le\t%le\t%le\t%le\t%le\t%le\t%le\n",CContinuum::anu[iCl],CIntegration::flux_intrinsic[iCl]*factor, factor*CIntegration::flux_continua[iCl]*App::cov_fac, factor*CIntegration::flux_transitions[iCl], bandInten, trans,bandInten*factor/2.179874099E-11);
					}
					else
					if(i == 11)
					{
						double factor = pow(App::distance/CGeometry::outer_radius[0][CGeometry::nRadiuses[0]],2.0)*CIntegration::lumfactor;
						fprintf(FDMP2,"%le\t%le\t%le\t%le\t%le\t%le\t%le\n",CContinuum::anu[iCl],CIntegration::flux_intrinsic[iCl]*factor, factor*CIntegration::flux_continua[iCl]*App::cov_fac, factor*CIntegration::flux_transitions[iCl], bandInten, trans,bandInten*factor/2.179874099E-11);
					}
				}

				if(i == 1)
				{
					fclose(FDMP);
				}
				if(i == 11)
				{
					fclose(FDMP2);
				}

				bandInten *= CIntegration::lumfactor;
				Yan *= CIntegration::lumfactor;

				fprintf(FP, "%s\t\t%le\t\t%le\t\t%lf\t\t%lf\n", cBand.label, bandInten, Yan, log10(1.e-50 + 4*M_PI*bandInten*pow(App::distance,2.0)),bandInten/(CIntegration::flux_lines[0]*CIntegration::lumfactor));
			}

			fclose(FP);	
		}

		printf("Dist: %le\n", App::distance);

		CIteration::nIteration++;	
	}


	CIteration::obtainStatistics();

}