#include "basics.h"
#include "geom3d.h"
#include "solver.h"
#include "rt.h"
#include "geometry.h"
#include "continuum.h"
#include "lines.h"
#include "app.h"
#include "integration.h"

double *CRT::ray_cont;
double *CRT::ray_lines;
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
			double emct = 0.0;
			if(ilayer>=0)
			{
				emct = CContinuum::emits[sector][ilayer][ic];
				locEmiss = (CContinuum::emits[sector][ilayer][ic]/(4.0*M_PI))*CSolver::paths[ish].path;
			}
			
			CRT::ray_cont[ic]  *= dRFac;

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
			CRT::ray_cont[ic] *= opfac;
			CRT::ray_cont[ic] += locEmiss*dtau;
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

			CRT::ray_lines[il] += locEmiss;
			
			/*if(il==4)
			{
				printf("dR: %le\t; RL: %le; opfac: %le\n", dRFac, CRT::ray_lines[il],opfac);
			}*/
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