#include "basics.h"
#include "geom3d.h"
#include "solver.h"
#include "rt.h"
#include "geometry.h"
#include "continuum.h"
#include "lines.h"
#include "app.h"
#include "integration.h"
#include "abund.h"

double *CRT::ray_cont;
double *CRT::ray_lines;
double *CRT::ray_trans;
double *CRT::ray_intr;
double CRT::dS;

int CRT::get_mesh_ind(double en)
{
	int im=0;
	while(im<CContinuum::cellCount && en>CContinuum::anu[im])
	{
		im++;
	}

	return im;
}

int CRT::calc_ray(int mode)
{

	//struct sysinfo memInfo;
	//sysinfo(&memInfo);
	//printf("get points\n");
	//long long free_before = memInfo.freeram*memInfo.mem_unit;
	printf("Mode: %d\n", mode);
	CSolver::getPoints();
	//sysinfo(&memInfo);
	//long long free_after = memInfo.freeram*memInfo.mem_unit;
	//CIntegration::memPerGeom += free_before-free_after;
	//printf("got\n");
	CRT::calc_cont(mode);
	//sysinfo(&memInfo);
	//long long free_after2 = memInfo.freeram*memInfo.mem_unit;
	//CIntegration::memPerRT += free_after-free_after2;
	
}

int CRT::calc_cont(int mode)
{
	//printf("Calc cont\n");
	
	//FILE *RT;

	//RT = fopen("ray.info","w+");
	for(int ic=0;ic<CContinuum::cellCount; ic++)
	{
		CRT::ray_cont[ic] = 0.0;
		CRT::ray_trans[ic] = 0.0;
		CRT::ray_intr[ic] = 0.0;
	}

	for(int il=0;il<CLine::linesCount;il++)
	{
		CRT::ray_lines[il] = 0.0;
	}
	int nshells = CSolver::npoints-1;
	//printf("Reaady calc cont: %ld",CSolver::npoints);
	double path_total = 0.;
	double path_prev = 0.0;

	bool iCenterPassed = false;
	double dRFacTotl = 1.0;

	for(int ish=nshells;ish>=0;--ish)
	{
		int ilayer = CSolver::paths[ish].layer;
		int sector = CSolver::paths[ish].sector;
		//printf("%d - %d; %d; \n", ish, ilayer,sector);
		//getting continua emiss
		if(ilayer == -1)
		{
			iCenterPassed = true;

			if(mode == 0)
			{
				for(int ic = 0; ic<CContinuum::cellCount; ic++)
				{
					CRT::ray_intr[ic] = CContinuum::in_fluxes[ic]*App::cov_fac;

				}
				
				path_prev = 0.;
				path_total = 0.;
			}

			printf("Passed center\n");

			for(int ic = 0; ic<CContinuum::cellCount; ic++)
			{
				CRT::ray_intr[ic] = CContinuum::in_fluxes[ic]*App::cov_fac;

			}
		}

		if(mode==0 && !iCenterPassed)
		{
			continue;
		}

		double dRFac = 1.0;

		path_total += CSolver::paths[ish].path;

		if(path_prev>0.)
		{
			dRFac = pow(path_prev/path_total,2.0);
			
		}
		
		dRFacTotl *= dRFac;
		double opfac = 1.0;
			
		//fprintf(RT,"Path: %le; dR: %le; dRTot: %le\n",path_total, dRFac, dRFacTotl);
		//fprintf(RT,"%d CS: sector %d layer %d; path = %le\n",ish, sector, ilayer, CSolver::paths[ish].path);
		//fflush(RT);
		//continua
		for(int ic=0;ic<CContinuum::cellCount; ic++)
		{

			//getting Emissivity
			double locEmiss = 0.0;
			double transEmiss = 0.0;
			double emct = 0.0;
			if(ilayer>=0)
			{
				emct = CContinuum::emits[sector][ilayer][ic];
				locEmiss = (CContinuum::emits[sector][ilayer][ic]/(4.0*M_PI))*CSolver::paths[ish].path;
				transEmiss = (CContinuum::trans[sector][ilayer][ic]/(4.0*M_PI))*CSolver::paths[ish].path;
				
			}
			
			CRT::ray_cont[ic]  *= dRFac;
			CRT::ray_trans[ic] *= dRFac;
			CRT::ray_intr[ic] *= dRFac;

			opfac = 1.0;

			double tau = 0.;

			if(ilayer>=0)
			{
				tau = CSolver::paths[ish].path*CContinuum::opacs[sector][ilayer][ic];
				opfac = exp(-tau);
				
			}

			double dtau = 1.0;

			if(tau<1.e-4)
			{
				dtau = 1.0;
			}
			else
			{
				dtau = (1.0-exp(-tau))/tau;
			}
			//if(ic == 2240)
			//{
			//	fprintf(RT,"CONT: POEC = %le;(em= %le); opfc = %le; TOTL_B = %le\n",locEmiss, emct, opfac,CRT::ray_cont[ic]);
			//}

			if(mode == 0)
			{
				locEmiss *= 2*M_PI;
			}

			CRT::ray_cont[ic] *= opfac;
			CRT::ray_cont[ic] += locEmiss*dtau;

			CRT::ray_trans[ic] *= opfac;
			CRT::ray_trans[ic] += transEmiss*dtau;

			CRT::ray_intr[ic] *= opfac;
			//if(ic == 2240)
			//{
			//	fprintf(RT,"CONT: POEC = %le; opfc = %le; TOTL_B = %le\n",locEmiss, opfac,CRT::ray_cont[ic]);
			//	fflush(RT);
			//}
			
			//CContinuum::ray_cont[ic] = 0.0;
		}

		//printf("CONT_LOC: %le; dRFac: %le (%le/%le)^2\n",CRT::ray_cont[0],dRFac,path_prev,path_total);
	
		//lines
		for(int il=0;il<CLine::linesCount;il++)
		{
			double locEmiss = 0.0;
			if(ilayer>=0)
			{
				locEmiss = (CLine::emits[sector][ilayer][il]/(4*M_PI))*CSolver::paths[ish].path;
			}
			CRT::ray_lines[il]  *= dRFac;
			
			int lId = CLine::lineIds[il];

			double en = CLine::linesEn[lId];

			if(en<=0)
			{
				en = 1.0;
			}

			opfac = 1.0;
			int icell = CRT::get_mesh_ind(en);
			//printf("sec: %d; lay: %d; icell: %d\n", sector,ilayer,icell);
			//printf("opacs: %le\n", CContinuum::opacs[sector][ilayer][icell]);
			if(ilayer>=0)
			{
				opfac = exp(-CSolver::paths[ish].path*CContinuum::opacs[sector][ilayer][icell]);
			}
			CRT::ray_lines[il]  *= opfac;

			if(mode == 0)
			{
				locEmiss *= 4*M_PI*App::cov_fac;
			}

			CRT::ray_lines[il] += locEmiss;
			
			/*if(il==4)
			{
				printf("dR: %le\t; RL: %le; opfac: %le\n", dRFac, CRT::ray_lines[il],opfac);
			}*/
		}

		/*FILE *FP2;
		char fname2[255];
		sprintf(fname2,"%s/zonal.dat",App::output_dir);
		
		FP2 = fopen(fname2,"a+");*/
		
		double cn = 1.0;

		double ems = 0.;

		if(sector > -1 && ilayer > -1)
		{
			ems = CContinuum::emits[sector][ilayer][1136];
		}

		/*fprintf(FP2, "%le\t%le\t%le\t%le\t%le\n",path_prev,CRT::ray_intr[1136]*cn, CRT::ray_cont[1136]*cn, ems,CSolver::paths[ish].path);
		

		fclose(FP2);	*/

		if(!App::isStatMode)
		{
			double Z_H = 0.;
			double ems_O = 0.;

			if(sector > -1 && ilayer > -1)
			{
				Z_H = Abund::abundances[sector][ilayer][0];
				ems_O = CLine::emits[sector][ilayer][3] + 
				CLine::emits[sector][ilayer][4] +
				CLine::emits[sector][ilayer][5] +
				CLine::emits[sector][ilayer][10];
			}

			for(int i=0;i<Abund::nElements;i++)
			{
				double Z = 0.;

				if(sector > -1 && ilayer > -1)
				{
					Z = Abund::abundances[sector][ilayer][i];
				}
				Abund::abmass[i] += Z*Z_H;
				Abund::abemso[i] += Z*ems_O;
			}
			Abund::mass += Z_H;
			Abund::emso += ems_O;
		}

		path_prev = path_total;

	}
	
	//fclose(RT);
	/*for(int ic=0;ic<CContinuum::cellCount-1200; ic++)
	{
		//getting Emissivity

		printf("%le => %le\n", CContinuum::anu[ic], CRT::ray_cont[ic]);
	}*/

	
	//exit(1);

	return 1;
}