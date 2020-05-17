#include "basics.h"
#include "rt.h"
#include "geometry.h"
#include "continuum.h"
#include "lines.h"
#include "integration.h"
#include "iterator.h"

double *CIteration::line_fluxes;
double *CIteration::line_lumunosity;
double *CIteration::continuum_fluxes;
double *CIteration::totl_lines;
int CIteration::nIteration;
int CIteration::mode;
double *CIteration::line_fluxes_prev;
double *CIteration::continuum_fluxes_prev;

int CIteration::obtainStatistics()
{
	CIntegration::statisticRay();
}

int CIteration::initIteration(int mode)
{
	printf("initIteration\n");
	CIteration::line_fluxes = new double[CLine::linesCount];
	CIteration::line_lumunosity = new double[CLine::linesCount];
	CIteration::continuum_fluxes = new double[CContinuum::cellCount];
	CIteration::totl_lines = new double[CLine::linesCount];
	CIteration::line_fluxes_prev = new double[CLine::linesCount];
	CIteration::continuum_fluxes_prev = new double[CContinuum::cellCount];

	CIteration::nIteration = 1;
	CIteration::mode = mode;
	for(int i=0;i<CLine::linesCount;i++)
	{
		CIteration::line_fluxes_prev[i] = 0.0;
	}

	for(int i=0;i<CContinuum::cellCount;i++)
	{
		CIteration::continuum_fluxes_prev[i] = 0.0;
	}

	CIntegration::flux_continua = new double[CContinuum::cellCount];
	CIntegration::flux_lines = new double[CLine::linesCount];
	CIntegration::flux_transitions = new double[CContinuum::cellCount];
	CIntegration::flux_intrinsic = new double[CContinuum::cellCount];
	
	CRT::ray_cont = new double[CContinuum::cellCount];
	CRT::ray_lines = new double[CLine::linesCount];
	CRT::ray_trans = new double[CContinuum::cellCount];
	CRT::ray_intr = new double[CContinuum::cellCount];
}

int CIteration::doIteration()
{
	printf("doIteration");
	if(CIteration::nIteration>1)
	{
		for(int i=0;i<CLine::linesCount;i++)
		{
			CIteration::line_fluxes_prev[i] = CIteration::line_fluxes[i];
		}

		for(int i=0;i<CContinuum::cellCount;i++)
		{
			CIteration::continuum_fluxes_prev[i] = CIteration::continuum_fluxes[i];
		}
	}

	CIntegration::InitIntegration(CIteration::mode,6+4*nIteration);
	CIntegration::doCalc();

	for(int i=0;i<CLine::linesCount;i++)
	{
		CIteration::line_fluxes[i] = CIntegration::flux_lines[i];
	}

	for(int i=0;i<CContinuum::cellCount;i++)
	{
		CIteration::continuum_fluxes[i] = CIntegration::flux_continua[i];
	}

	CIntegration::finish();
	//exit(1);
}

double CIteration::getdCont()
{
	double dlt = 0.0;
	if(CIteration::nIteration < 2)
	{
		return 1.0;
	}

	for(int i=0;i<CContinuum::cellCount;i++)
	{
	
		double locDlt = pow((CIteration::continuum_fluxes[i] - CIteration::continuum_fluxes_prev[i])/(CIteration::continuum_fluxes_prev[i]+1.e-35),2.0);
		//if(locDlt > 1.e+50)
		//{
		//	printf("v=%le; next: %le; prev: %le\n",CContinuum::anu[i],CIteration::continuum_fluxes[i],CIteration::continuum_fluxes_prev[i]);
		//}
		
		dlt += locDlt;
	}

	dlt = dlt/CContinuum::cellCount;

	return dlt;
}

double CIteration::getdLine()
{
	double dlt = 0.0;
	if(CIteration::nIteration < 2)
	{
		return 1.0;
	}

	for(int i=0;i<CLine::linesCount;i++)
	{
		dlt += pow((CIteration::line_fluxes[i] - CIteration::line_fluxes_prev[i])/(CIteration::line_fluxes_prev[i]+1.e-80),2.0);
	}

	dlt = dlt/CLine::linesCount;

	return dlt;
}