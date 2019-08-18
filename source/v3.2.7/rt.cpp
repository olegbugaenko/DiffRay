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
#include "grain_temp.h"
#include "physics.h"

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

	CSolver::getPoints();
	CRT::calc_cont(mode);
	
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

	double *face_opacities;
	double *total_opacities;

	face_opacities = new double[CLine::linesCount];
	total_opacities = new double[CLine::linesCount];
	
	int lowp = CSolver::npoints-1;

	if(mode == 0)
	{
		for(int ish=CSolver::npoints-1;ish>=0;--ish)
		{
			int ilayer = CSolver::paths[ish].layer;
			if(ilayer == -1)
			{
				lowp = ish;
				break;
			}
		}
	}

	for(int il=0;il<CLine::linesCount;il++)
	{
		CRT::ray_lines[il] = 0.0;
		CLine::cumulative[il] = 0.0;
		face_opacities[il] = 0.0;
		total_opacities[il] = CRT::TauAbsToPoint(il,lowp);
	}

	int nshells = CSolver::npoints-1;
	//printf("Reaady calc cont: %ld",CSolver::npoints);
	double path_total = 0.;
	double path_prev = 0.0;

	bool iCenterPassed = false;
	double dRFacTotl = 1.0;

	int nump = 0;

	for(int ish=nshells;ish>=0;--ish)
	{
		int ilayer = CSolver::paths[ish].layer;
		int sector = CSolver::paths[ish].sector;

		
		//printf("%d of %d - %d; %d; \n", ish, nshells, ilayer,sector);
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


			for(int ic = 0; ic<CContinuum::cellCount; ic++)
			{
				CRT::ray_intr[ic] = CContinuum::in_fluxes[ic]*App::cov_fac;

			}
		}

		if(mode==0 && !iCenterPassed)
		{
			continue;
		}

		double dVeff = 0.;

		if(ilayer > 0)
		{
			double Rout = CGeometry::outer_radius[sector][ilayer];
			double Rin = CGeometry::outer_radius[sector][ilayer-1];
			double vin = pow(Rin/CGeometry::inRadius[sector],2)*Rin;
			double vout = pow(Rout/CGeometry::inRadius[sector],2)*Rout;
			dVeff = (vout - vin)/3.*App::cov_fac;
		}

		//printf("%d dVeff = %le; Hop=%le; (lowp = %d/%d)\n", ilayer, dVeff, total_opacities[0], lowp, CSolver::npoints);

		double dRFac = 1.0;

		path_total += CSolver::paths[ish].path;

		if(path_prev>0.)
		{
			dRFac = pow(path_prev/path_total,2.0);
			
		}

		//printf("OK1\n");
		
		dRFacTotl *= dRFac;
		double opfac = 1.0;
			
		
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
			
		}

		//printf("CONT_LOC: %le; dRFac: %le (%le/%le)^2\n",CRT::ray_cont[0],dRFac,path_prev,path_total);
		//printf("OK2\n");
		
		//lines
		for(int il=0;il<CLine::linesCount;il++)
		{
			//!!!!!!!!!Debug if ish = 5!!!!!!!!!!!!!
			if(ish == 6){
				//printf("OK2.1: sec %d; lay %d\n", sector, ilayer);
			}
			double locEmiss = 0.0;
			if(ilayer>=0)
			{
				if(CLine::emits[sector][ilayer][il] > 1.18e-36)
					locEmiss = (CLine::emits[sector][ilayer][il]/(4*M_PI))*CSolver::paths[ish].path;
			}
			if(ish == 6){
				//printf("OK2.2\n");
			}
			CRT::ray_lines[il]  *= dRFac;
			
			int lId = CLine::lineIds[il];
			if(ish == 6){
				//printf("OK2.3\n");
			}
			double en = CLine::linesEn[lId];

			if(en<=0)
			{
				en = 1.0;
			}

			opfac = 1.0;
			double tau_this = 0.0;
			int icell = CRT::get_mesh_ind(en);
			if(ish == 6){
				//printf("OK2.4\n");
			}
			if(ilayer>=0)
			{
				tau_this = CSolver::paths[ish].path*CContinuum::opacs[sector][ilayer][icell];
				opfac = exp(-tau_this);
			}
			CRT::ray_lines[il]  *= opfac;
			if(ish == 6){
				//printf("OK2.5\n");
			}
			if(mode == 0)
			{
				locEmiss *= 4*M_PI*App::cov_fac;
			}

			CRT::ray_lines[il] += locEmiss;

			face_opacities[il] += tau_this;

			if(ish == 6){
				//printf("OK2.6\n");
			}
			double dTau = total_opacities[il] - face_opacities[il]+tau_this;
			if(dTau<0)
			{
			    dTau = 1.e-36;
			    //printf("dTau is negative for line %d\n",i);
			    //printf("Face: %le; Total: %le; tau_this: %le\n",total_opacities[i],face_opacities[i],tau_this);
			}
			double line_opfc = CLine::e2(dTau);
			double line_inopfc = CLine::e2(face_opacities[il])*CLine::e2(total_opacities[il]);
	
			if(ish == 6){
				//printf("OK2.7\n");
			}
			if(ilayer>=0)
			{
				CLine::cumulative[il] += dVeff*(0.5*line_opfc + 0.5*line_inopfc)*CLine::emits[sector][ilayer][il];
				//printf("E: %le V: %le Tot: %le\n",dVeff, CLine::emits[sector][ilayer][il],CLine::cumulative[il]);
			}

			
		}

		

		if(!App::isStatMode)
		{
			double Z_H = 0.;
			double ems_O = 0.;

			if(sector > -1 && ilayer > -1 && ilayer < Abund::nrows)
			{
				//printf("Layer: %d\n");
				Z_H = CSolver::paths[ish].path*Abund::abundances[sector][ilayer][0];
				ems_O = CSolver::paths[ish].path*(CLine::emits[sector][ilayer][3] + 
				CLine::emits[sector][ilayer][4] +
				CLine::emits[sector][ilayer][5] +
				CLine::emits[sector][ilayer][10]);
			}

			double Z = 0.;

			for(int i=0;i<Abund::nElements;i++)
			{
				Z = 0.;

				if(sector > -1 && ilayer > -1 && ilayer < Abund::nrows)
				{
					
					Z = Abund::abundances[sector][ilayer][i]/Abund::abundances[sector][ilayer][0];
					Abund::abmass[i] += Z*Z_H;
					Abund::abemso[i] += Z*ems_O;
				}
				
			}
			Abund::mass += Z_H;
			Abund::emso += ems_O;
			
		}
		
		if(App::isStatMode || App::punchStatistics)
		{
			Abund::statistics[CLine::iStat][0] = path_total;
			CLine::statistics[CLine::iStat][0] = path_total;
			GrainTemp::statistics[CLine::iStat][0] = path_total;		
			
			for(int i=0;i<Abund::nElements;i++)
			{
				double Z = 0.;

				if(sector > -1 && ilayer > -1)
				{

					
					Z = Abund::abundances[sector][ilayer][i];
					Abund::statistics[CLine::iStat][i+1] = Z;
				}

				if(i==0 && Z > pow(10.0, -2.0) && App::iApp == 5)
				{
					printf("App%d peak: %le at %d;%d path=%le\n",App::iApp,Z,sector,ilayer,CSolver::paths[ish].path);
					printf("iRay: %d; point: %d of %d\n", CIntegration::nRay, ish, CSolver::npoints);
					printf("Struct data: \n");
					printf("%d: Sector: %d; Layer: %d; P=%le\n", ish, CSolver::paths[ish].sector, CSolver::paths[ish].layer,CSolver::paths[ish].path);
				}
				
			}
			

			for(int i=0;i<CLine::linesCount;i++)
			{
				double Z = 0.;

				if(sector > -1 && ilayer > -1 && ilayer < CLine::nrows)
				{
					
					Z = CLine::emits[sector][ilayer][i];
					CLine::statistics[CLine::iStat][i+1] = Z;
				}
				
			}
			
			
			for(int i=0;i<GrainTemp::nBins;i++)
			{
				double Z = 0.;

				if(sector > -1 && ilayer > -1)
				{
					//printf("S: %d/19, L: %d, B: %d (%d)\n", sector, ilayer, i, GrainTemp::nBins);
					Z = GrainTemp::gridTemps[sector][ilayer][i];
					/*if(i == GrainTemp::nBins-1)
					{
						printf("%d (%d) GT: %le\n", CLine::iStat, i, Z);
					}*/
					
					GrainTemp::statistics[CLine::iStat][i+1] = Z;
				}
				
			}
/*
			if(ilayer > -1 && sector > -1){
				printf("====\nA[%le]=%d of %d;%d;:%le\n",Abund::statistics[CLine::iStat][0],ilayer,Abund::nrows,sector,Abund::statistics[CLine::iStat][1]);
			}
			else
			{
				printf("====\nA[%le]=%d;%d:0.0\n",Abund::statistics[CLine::iStat][0],ilayer,sector);
			}*/
			
			CLine::iStat++;

		}

		path_prev = path_total;
		nump++;

	}
	return 1;
}

double CRT::TauAbsTotal(int iline)
{
	double en = CLine::linesEn[iline];

	if(en<=0)
	{
		en = 1.0;
	}

	int icell = CRT::get_mesh_ind(en);
	
	double tau = 1.e-40;

	if(icell<0 || icell>CContinuum::cellCount-1)
	    return 0.;
	
		
	for(int i=0;i<CSolver::npoints-2;i++)
	{
		int nlw = CSolver::paths[i].layer;
		int nsc = CSolver::paths[i].sector;
		
		tau += CContinuum::opacs[nsc][nlw][icell]*CSolver::paths[i].path;
	}
	
	return tau;
}


double CRT::TauAbsTotalCont(int icont)
{
	double tau = 1.e-40;
	
		
	for(int i=0;i<CSolver::npoints-2;i++)
	{
		int nlw = CSolver::paths[i].layer;
		int nsc = CSolver::paths[i].sector;
		
		tau += CContinuum::opacs[nsc][nlw][icont]*CSolver::paths[i].path;
	}
	
	return tau;
}


double CRT::TauAbsToPoint(int iline, int lowp)
{
	double en = CLine::linesEn[iline];

	if(en<=0)
	{
		en = 1.0;
	}

	int icell = CRT::get_mesh_ind(en);
	
	if(icell<0 || icell>CContinuum::cellCount-1)
	    return 0.;
	
	double tau = 1.e-40;
	
		
	for(int i=0;i<lowp;i++)
	{
		int nlw = CSolver::paths[i].layer;
		int nsc = CSolver::paths[i].sector;
		//printf("point: %d; mlw=%d\n",i,nlw);
		if(nlw>=0)
		{
			double dtau = CContinuum::opacs[nsc][nlw][icell]*CSolver::paths[i].path;
			
			tau += dtau;
		}
	}
	
	return tau;
}


double CRT::TauAbsToPointCont(int icont, int lowp)
{
	
	double tau = 1.e-40;
	
		
	for(int i=0;i<lowp;i++)
	{
		int nlw = CSolver::paths[i].layer;
		int nsc = CSolver::paths[i].sector;
		if(nlw>=0)
		{
			double dtau = CContinuum::opacs[nsc][nlw][icont]*CSolver::paths[i].path;
			if(dtau>1.e-40)
			tau += dtau;
		}
	}
	
	return tau;
}