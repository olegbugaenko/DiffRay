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
#include "grain_temp.h"
#include "integration/angular.h"

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
bool CIntegration::isCentralRay = false;

bool CIntegration::iterateOverSource(double phi, double theta, double dphi, double dtheta)
{
	if(App::isIrregularSource)
		return true;

	//if source within app position, and fit it comletely
	double objectWidth = 2 * CGeometry::inRadius[0];
	if(dphi * App::distance < objectWidth)
		return false;
	if(dtheta * App::distance < objectWidth)
		return false;

	while(phi > M_PI)
		phi -= 2*M_PI;

	if((phi - dphi/2.) * App::distance < objectWidth / 2.0 &&
	    (phi + dphi/2.) * App::distance > objectWidth / 2.0 &&
	    (theta - dtheta/2.) * App::distance < objectWidth / 2.0 &&
	    (theta + dtheta/2.) * App::distance > objectWidth / 2.0)
	{
		return false;
	}

	return true;
	
}

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
	for(int i=0;i<100;i++)
	{
		GrainTemp::massWeighted[i] = 0.0;
	}
	GrainTemp::massTotal = 0.0;
	Abund::mass = 0.;
	Abund::emso = 0.;
	CIntegration::nRay = 0;

}

int CIntegration::doCalc()
{
	//printf("do Calc\n");
	
	App::rayIntegration = App::rayToObj;
	double lumfactor = 1.0;
				
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

		double STOT = 0.0;
		int iskip_phi = 0;
		int iskip_theta = 0;
		int skip_param = 5;
		bool doIterationsOverSource = CIntegration::iterateOverSource(App::AppDPhi, App::AppDTheta, App::phi_width, App::theta_width);
		if(!doIterationsOverSource)
		{
			//init intrinsic from one central ray
			App::rayIntegration.angle.y += App::AppDPhi;
			App::rayIntegration.angle.z += App::AppDTheta;
			CRT::obtainIntrinsic();
			/*for(int ic = 0; ic < CContinuum::cellCount; ic++)
			{
				if(ic < 840)
					printf("F=%le\n",CRT::ray_intr[ic]);
				CIntegration::flux_intrinsic[ic] = CRT::ray_intr[ic];
			}*/
			//exit(1);
		}
		App::rayIntegration = App::rayToObj;
		CBasics::startClock();
		Angular* stepIterator = new Angular(App::rayIntegration.angle.y + App::AppDPhi, App::rayIntegration.angle.z + App::AppDTheta, App::phi_width, App::theta_width, 5, doIterationsOverSource);
		stepIterator->iterate(doIterationsOverSource);
		CBasics::endClock();
		printf("STEP_S_TOT: %le; Scen: %le; Vcen: %le; [%le]\n",stepIterator->STot, stepIterator->Scenter, stepIterator->Vcenter, stepIterator->SMcenter);
		printf("phi: %le - %le; theta: %le - %le;\n",App::minPhiCen, App::maxPhiCen, App::minThetaCen, App::maxThetaCen);
		/*
		while(phi<phi_max+dphi/5.0)
		{
			
			theta = App::AppDTheta - (App::theta_width-dtheta)/2;
			iskip_theta = 0;
			while(theta<theta_max+dtheta/5.0)
			{
				if(fabs(App::AppDPhi - phi) < dphi && fabs(App::AppDTheta - theta) < dtheta)
					CIntegration::isCentralRay = true;
				else
					CIntegration::isCentralRay = false;
				//Each 4'th ray should trigger statistics
				//At least, central and eadge angles should be included

				//phi = 3.179-M_PI;
				//theta = 0.0;
				//App::setRay(phi,theta);

				App::punchStatistics = false;
				if(iskip_phi == 0 || iskip_phi == CIntegration::nPhi - 1
					|| iskip_theta == 0 || iskip_theta == CIntegration::nTheta - 1 
					|| iskip_phi == (int)(CIntegration::nPhi/2) 
					|| iskip_theta ==(int)(CIntegration::nTheta/2))
				{
					if(skip_param == 5)
					{
						App::punchStatistics = true;
						skip_param = 0;
					}
					skip_param++;
					
				}

				CIntegration::nRay++;
				App::rayIntegration = App::rayToObj;

				App::rayIntegration.angle.y += phi;
				App::rayIntegration.angle.z += theta;
				
				//App::rayIntegration.angle.y = 0.5;
				//App::rayIntegration.angle.z = 0.5;
				

				CRT::dS = dphi*dtheta*cos(theta);
				STOT += CRT::dS;
				//CBasics::startClock();
				CRT::calc_ray(1);
				//CBasics::endClock();
				for(int i=0;i<CLine::linesCount;i++)
				{
					CIntegration::flux_lines[i] += CRT::ray_lines[i]*CRT::dS;
				}

				//printf("LINES BUG: %le <<<< %le?\n",CIntegration::flux_lines[2],CIntegration::flux_lines[1]);

				for(int i=0;i<CContinuum::cellCount;i++)
				{
					CIntegration::flux_continua[i] += CRT::ray_cont[i];//*CRT::dS;
					CIntegration::flux_transitions[i] += CRT::ray_trans[i]*CRT::dS;
					CIntegration::flux_intrinsic[i] += CRT::ray_intr[i];//*CRT::dS;
					if(fabs(CContinuum::anu[i]-1.0)<0.005)
					{
						//fprintf(TEST_CONT, "%le %le (+%le from %le;%le; dS=%le (%le;%le)) [%le;%le;%le]\n", CContinuum::anu[i], CIntegration::flux_continua[i],CRT::ray_cont[i],phi,theta,CRT::dS,dphi,dtheta,App::rayIntegration.angle.x,App::rayIntegration.angle.y,App::rayIntegration.angle.z);
					}
				}
				//printf("%le\t%le\t%le\t%le\t%le\t%le\n",CContinuum::anu[2240],phi,theta,CIntegration::flux_continua[843],CRT::ray_cont[843],CRT::dS);
				theta += dtheta;
				//exit(1);
				iskip_theta++;
				
			}
			phi += dphi;
			iskip_phi++;
			
		}*/

		
		lumfactor = 1.0/(STOT);
		
	

	}
	
	
	

	CIntegration::lumfactor = lumfactor;
	printf("INTEGRATION_ENDED\n");
}

int CIntegration::finish()
{
	//delete[] CIntegration::flux_continua;
	//delete[] CIntegration::flux_lines;
}