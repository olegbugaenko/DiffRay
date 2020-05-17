#pragma once
#include "../basics.h"
#include "integrationRay.h"
#include "../continuum.h"
#include "../lines.h"
#include "../rt.h"
#include "../app.h"
#include "../solver.h"
#include "../const.h"

class CMatrix{
    public:
	static long double phi, theta, dphi, dtheta, dS;
	static long double calcPhi, calcTheta, calcDphi, calcDtheta, calcDS;
	static CIntegrationRay *rays[9];
	static double *ray_intr, *ray_cont, *ray_lines;
	static long double FFac;
	static double Scenter, Vcenter;

	static int setup(long double cphi, long double ctheta, long double cdphi, long double cdtheta)
	{
		printf("Init: %Le;%Le; %Le X %Le\n", cphi, ctheta, cdphi, cdtheta);
		CMatrix::phi = cphi;
		CMatrix::theta = ctheta;
		CMatrix::dphi = cdphi;
		CMatrix::dtheta = cdtheta;
		CMatrix::dS = dphi*dtheta;
		CMatrix::calcPhi = cphi;
		CMatrix::calcTheta = ctheta;
		CMatrix::calcDphi = cdphi;
		CMatrix::calcDtheta = cdtheta;
		CMatrix::calcDS = dphi*dtheta;
		CMatrix::Scenter = 0.;
		CMatrix::Vcenter = 0.;
		for(int ic = 0;ic < CContinuum::cellCount; ic++)
		{
			CMatrix::ray_intr[ic] = 0.;
			CMatrix::ray_cont[ic] = 0.;
		}
		CMatrix::initRays();
		return 1;
	}

	static int initRays()
	{
		
		for(int i = 0; i < 9; i++)
		{
			CMatrix::rays[i] = new CIntegrationRay(
				CMatrix::calcPhi + ((i % 3)-1)*CMatrix::calcDphi/3,
				CMatrix::calcTheta + (floor(i/3)-1)*CMatrix::calcDtheta/3,
				CMatrix::calcDphi/3,
				CMatrix::calcDtheta/3
			);
			CMatrix::rays[i]->initMesh(CContinuum::cellCount);
			CMatrix::rays[i]->initLines(CLine::linesCount);
		}
		return 1;
	}

	static int run(bool bWithStat, bool doIntrIter)
	{
		CMatrix::FFac = 0.;
		for(int i=0;i<9; i++)
		{
			App::rayIntegration = App::rayToObj;

			App::rayIntegration.angle.y = CMatrix::rays[i]->phi;
			App::rayIntegration.angle.z = CMatrix::rays[i]->theta;
			CRT::dS = CMatrix::rays[i]->dS;
			if(i == 4 && bWithStat)
			    App::punchStatistics = true;
			else
			    App::punchStatistics = false;
			CRT::calc_ray(1);
			for(int ic = 0; ic < CContinuum::cellCount; ic++)
			{
				if(doIntrIter)
					CMatrix::rays[i]->ray_intr[ic] = CRT::ray_intr[ic];
				
				CMatrix::rays[i]->ray_cont[ic] = CRT::ray_cont[ic];
			}
			for(int il = 0; il < CLine::linesCount; il++)
			{
				CMatrix::rays[i]->ray_line[il] = CRT::ray_lines[il];
			}
			if(CRT::Hcenter > 1.e-50)
			{
				CMatrix::Vcenter += CRT::Hcenter;
				CMatrix::Scenter += CRT::dS*App::distance*App::distance;
				if(App::maxPhiCen < CMatrix::rays[i]->phi)
					App::maxPhiCen = CMatrix::rays[i]->phi;
				if(App::minPhiCen > CMatrix::rays[i]->phi)
					App::minPhiCen = CMatrix::rays[i]->phi;
				if(App::maxThetaCen < CMatrix::rays[i]->theta)
					App::maxThetaCen = CMatrix::rays[i]->theta;
				if(App::minThetaCen > CMatrix::rays[i]->theta)
					App::minThetaCen = CMatrix::rays[i]->theta;
				
			} else {
				/*if(fabs(CMatrix::rays[i]->theta) < 1.3e-6 && fabs(CMatrix::rays[i]->phi) < 1.3e-6)
				{
					double dist = sqrt(pow(CMatrix::rays[i]->theta,2.0) + pow(CMatrix::rays[i]->phi,2.0));
					if(dist < 1.6e-6) //-6
					{
						printf("Angles %Le %Le dont intersect center!!!\n",CMatrix::rays[i]->phi,CMatrix::rays[i]->theta);
						printf("Ray: %Le %Le %Le\n",
							App::rayIntegration.angle.x,
							App::rayIntegration.angle.y,
							App::rayIntegration.angle.z);
						printf("DIST: %le\n",dist);
						char fname[255];
						FILE *FP;
						sprintf(fname, "%s/intersections_debug.dat",App::output_dir);
						FP = fopen(fname,"w+");
						for(i=0;i<CSolver::npoints;i++)
						{
							fprintf(FP,"Sector: %d/%d (%d - %lf); Layer: %d/%d; Point: %Le;%Le;%Le; R=%le\n", CSolver::points[i].sec_p, CSolver::points[i].sec_n,CSolver::points[i].is,CSolver::points[i].phi, CSolver::points[i].lay_p, CSolver::points[i].lay_n, CSolver::points[i].point.x, CSolver::points[i].point.y, CSolver::points[i].point.z,CSolver::points[i].R);
						}
						fclose(FP);

		sprintf(fname, "%s/paths.dat",App::output_dir);
		FP = fopen(fname,"w+");
		for(i=0;i<CSolver::npoints;i++)
		{
			fprintf(FP,"Sector: %d; Layer: %d; P=%le\n", CSolver::paths[i].sector, CSolver::paths[i].layer,CSolver::paths[i].path);
			//printvec(CSolver::points[i].point);
		}
		fclose(FP);


						exit(1);
				
					}
				}*/
			}
		
		}
		CMatrix::FFac = CMatrix::getRaysFi();
		return 1;
	}

	static void addFluxes(double &nPhots);

	static double getRaysFi();

	static bool getDelta(double nPhots)
	{
		printf("getDelta: %Le of %le\n",CMatrix::FFac, nPhots);
		for(int i = 0; i < 9; i++)
		{
			double deltaHere = fabs(1.e-50 + CMatrix::rays[i]->Fi - CMatrix::FFac)/(CMatrix::FFac + 1.e-45);
			//if more than 15% fluctuations within matrix found - there should be smthg interesting
			if(!CMatrix::isErrorAcceptable(nPhots, deltaHere)) {
				printf("%d [%Le;%Le]Delta here %le; %Le - %Le\n", i, CMatrix::rays[i]->phi, CMatrix::rays[i]->theta, deltaHere, CMatrix::rays[i]->Fi, CMatrix::FFac);
				return false;
			}
		}
		return true;
	}

	static bool freeMem()
	{
		for(int i = 0; i < 9; i++)
		{
			delete CMatrix::rays[i]->ray_intr;
			delete CMatrix::rays[i]->ray_line;
			delete CMatrix::rays[i]->ray_cont;
			delete CMatrix::rays[i];
		}
		return true;
	}

	static bool isErrorAcceptable(double nPhots, double deltaHere)
	{
		//Nothing to do here, fluxes are close to zero so far
		if(CMatrix::FFac * deltaHere <= SMALL_DELTA)
		    return true;


		if(deltaHere < 0.05)
		    return true;
		//Too great deviations. Even if fluxes are small comparing to 
		//total amount of photons - we should never accept deviations of order
		//since they can mean powerfull sources, weaken by small resolution
		if(deltaHere > 1.e+1)
		    return false;

		//get nPhots per radian
		//long double nPhotsPerRadian = nPhots / (App::phi_width * App::theta_width);
		if(CMatrix::dS/(App::phi_width * App::theta_width) < 9.e-5)
		{
			printf("TooSmall\n");
			return true;
		}
		long double deviation = CMatrix::FFac * deltaHere;

		//long double deviationPerRadian = deviation / CMatrix::dS / pow(1 + App::phi_width * App::theta_width / CMatrix::dS, 0.25);
		//if(deviation < 5.e-4 * nPhots)
		//	return true;
		//now, if deviationPerRadian is comparable to average photons
		//density - recalculate
		printf("%Le > 5.e-4 * %le?\n",deviation, nPhots);
		/*if(deviationPerRadian >= 0.1 * nPhotsPerRadian)
		{
			return false;
		}*/
		if(deviation > 5.e-4*(nPhots + CMatrix::FFac))
			return false;
		
		if(CMatrix::dS/(App::phi_width * App::theta_width) > 0.01)
		{
			printf("TooLarge\n");
			return false;
		}
		
		return true;
	}
};