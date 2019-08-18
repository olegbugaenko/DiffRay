#include "basics.h"
#include "geom3d.h"
#include "solver.h"
#include "rt.h"
#include "geometry.h"
#include "continuum.h"
#include "lines.h"
#include "integration.h"
#include "app.h"
#include "math.h"
#include "abund.h"

int CIntegration::mode = 1;
double *CIntegration::flux_lines;
double *CIntegration::flux_continua;
double *CIntegration::flux_transitions;
double *CIntegration::flux_intrinsic;
int CIntegration::nPhi = 1;
int CIntegration::nTheta = 1;
double CIntegration::lumfactor;
long long CIntegration::memPerGeom = 0;
long long CIntegration::memPerRT = 0;
int CIntegration::nRay = 0;

int CIntegration::statisticRay()
{
	App::isStatMode = true;
	App::rayIntegration = App::rayToObj;

	App::rayIntegration.angle.y += App::AppDPhi;
	App::rayIntegration.angle.z += App::AppDTheta;

	for(int i=0;i<Abund::nElements;i++)
	{
		Abund::abmass[i] = 0.;
		Abund::abemso[i] = 0.;
	}
	Abund::mass = 0.;
	Abund::emso = 0.;

	CRT::calc_ray(1);
	App::isStatMode = false;
}

int CIntegration::InitIntegration(int mode,int nPrecision)
{
	printf("Init Integration: %d,%d\n",mode,nPrecision);
	CIntegration::mode = mode;

	for(int i = 0; i < CContinuum::nbands; i++)
	{
		CContinuum::bands[i].totalLum = 0.;
	}

	for(int ic=0;ic<CContinuum::cellCount; ic++)
	{
		CIntegration::flux_continua[ic] = 0.0;
		CIntegration::flux_transitions[ic] = 0.0;
		CIntegration::flux_intrinsic[ic] = 0.0;
	}

	for(int il=0;il<CLine::linesCount;il++)
	{
		CIntegration::flux_lines[il] = 0.0;
		CLine::cumulative[il] = 0.0;
	
	}

	if(mode>0) //not outward only
	{
		CIntegration::nPhi   = nPrecision*2+1;
		CIntegration::nTheta = nPrecision+1;
		//printf("PREC: %le/%d\n",M_PI,nPrecision);
		//printf("nPhi: %d\n",CIntegration::nPhi);
	}
	else
	{
		CIntegration::nPhi = 1;
		CIntegration::nTheta = 1;
	}

	for(int i=0;i<Abund::nElements;i++)
	{
		Abund::abmass[i] = 0.;
		Abund::abemso[i] = 0.;
	}
	Abund::mass = 0.;
	Abund::emso = 0.;
	CIntegration::nRay = 0;

	
}

int CIntegration::doCalc()
{
	//printf("do Calc\n");
	struct sysinfo memInfo;

	App::rayIntegration = App::rayToObj;
	double lumfactor = 1.0;
	FILE *FPDUMP;
	FPDUMP = fopen("dump.dat","a+");
	sysinfo(&memInfo);
	long long free_before = memInfo.freeram*memInfo.mem_unit;
				
	long long sum = 0;
	CIntegration::memPerGeom = 0;
	CIntegration::memPerRT = 0;

	if(CIntegration::mode == 0)
	{
		CIntegration::nRay++;
		CRT::dS = 4*M_PI/3;
		CRT::calc_ray(0);

		for(int i=0;i<CLine::linesCount;i++)
		{
			CIntegration::flux_lines[i] += CRT::ray_lines[i];
		}



		for(int i=0;i<CContinuum::cellCount;i++)
		{
			CIntegration::flux_continua[i] += CRT::ray_cont[i];
			CIntegration::flux_transitions[i] += CRT::ray_trans[i]*CRT::dS;
			CIntegration::flux_intrinsic[i] += CRT::ray_intr[i];
		}

		
		
	}
	else
	{
		//printf("nPhi: %d\n", CIntegration::nPhi);
		double dphi = App::phi_width/CIntegration::nPhi;
		double dtheta = App::theta_width/CIntegration::nTheta;
		printf("nPhi: %d; nTheta: %d\n",CIntegration::nPhi,CIntegration::nTheta);
		double phi = App::AppDPhi - (App::phi_width-dphi)/2;
		double theta = App::AppDTheta - (App::theta_width-dtheta)/2;
		
		double phi_max = App::AppDPhi + (App::phi_width-dphi)/2;
		double theta_max = App::AppDTheta + (App::theta_width-dtheta)/2;

		//printf("phi: %le - %le; dphi: %le\n", phi, phi_max, dphi);
		//printf("theta: %le - %le; dtheta: %le\n", theta, theta_max, dtheta);
		
		
		//FILE *TEST_CONT;

		//TEST_CONT = fopen("last_iter.dat","w+");
		//FILE *FPD;
		//FPD = fopen("A_Test_Int.dat","w+");
				

		double STOT = 0.0;
		while(phi<phi_max+dphi/5.0)
		{

			theta = App::AppDTheta - (App::theta_width-dtheta)/2;

			while(theta<theta_max+dtheta/5.0)
			{
				//phi = 3.179-M_PI;
				//theta = 0.0;
				//App::setRay(phi,theta);
				CIntegration::nRay++;
				App::rayIntegration = App::rayToObj;

				App::rayIntegration.angle.y += phi;
				App::rayIntegration.angle.z += theta;
				
				//App::rayIntegration.angle.y = 0.5;
				//App::rayIntegration.angle.z = 0.5;
				

				CRT::dS = dphi*dtheta*cos(theta);
				STOT += CRT::dS;
				CRT::calc_ray(1);
				for(int i=0;i<CLine::linesCount;i++)
				{
					CIntegration::flux_lines[i] += CRT::ray_lines[i]*CRT::dS;
				}



				for(int i=0;i<CContinuum::cellCount;i++)
				{
					CIntegration::flux_continua[i] += CRT::ray_cont[i]*CRT::dS;
					CIntegration::flux_transitions[i] += CRT::ray_trans[i]*CRT::dS;
					CIntegration::flux_intrinsic[i] += CRT::ray_intr[i]*CRT::dS;
					if(fabs(CContinuum::anu[i]-1.0)<0.005)
					{
						//fprintf(TEST_CONT, "%le %le (+%le from %le;%le; dS=%le (%le;%le)) [%le;%le;%le]\n", CContinuum::anu[i], CIntegration::flux_continua[i],CRT::ray_cont[i],phi,theta,CRT::dS,dphi,dtheta,App::rayIntegration.angle.x,App::rayIntegration.angle.y,App::rayIntegration.angle.z);
					}
				}
				//fprintf(FPD,"%le\t%le\t%le\t%le\t%le\t%le\n",CContinuum::anu[2240],phi,theta,CIntegration::flux_continua[2240],CRT::ray_cont[2240],CRT::dS);
				theta += dtheta;
				//exit(1);
			}
			phi += dphi;
		}

		//fclose(TEST_CONT);
		//printf("nPrecision: %d",CIntegration::nPhi);
		printf("S2TOT: %le; phi: %le-%le (%le); theta: %le-%le\n", STOT,App::phi,phi_max,dphi,App::theta,theta_max);
		//fclose(FPD);
		
		lumfactor = 1.0/(STOT);
		//exit(1);
		/*for(int i=0; i<CContinuum::cellCount;i++)
		{
			CIntegration::flux_continua[i] *= (1/(4*M_PI));
		}*/
		
	

	}
	
	
	
	fprintf(FPDUMP,"%d:geom=%ld;rt=%ld\n",CIntegration::nPhi,CIntegration::memPerGeom,CIntegration::memPerRT);
	fclose(FPDUMP);

	CIntegration::lumfactor = lumfactor;
	
}

int CIntegration::finish()
{
	//delete[] CIntegration::flux_continua;
	//delete[] CIntegration::flux_lines;
}