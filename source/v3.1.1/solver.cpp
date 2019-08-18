#include "basics.h"
#include "geom3d.h"
#include "solver.h"
#include "geometry.h"
#include "app.h"

int CSolver::npoints = 0;
intersection CSolver::points[3000];
TPath CSolver::paths[3000];

int normalizeSector(int phs)
{
	if(phs>=3*App::nSectors)
		phs = 4*App::nSectors-1-phs;
	else
	if(phs>=2*App::nSectors)
		phs = phs - 2*App::nSectors;
	else
	if(phs>=App::nSectors)
		phs = 2*App::nSectors-1-phs;
	
	if(phs<0)
		phs = 0;

	return phs;	
}

int CSolver::getPoints()
{
	//For the first time just use demo-ray
	TRay ray = App::rayIntegration;
	//ray.angle.y += 0.2;
	int npoints = 0;

	intersection curpoint;
	curpoint.point = App::rayIntegration.start;
	curpoint.sec_p = 0;
	curpoint.sec_n = 0;
	curpoint.lay_p = -2;
	curpoint.lay_n = -2;
	curpoint.is = 0;
	curpoint.phi = 0;
	curpoint.R = 0.0;

	CSolver::points[npoints] = curpoint;

	npoints++;
	//printvec(App::rayToObj.start);
	//printvec(App::rayToObj.angle);


	vector3 *iPoint;
	iPoint = new vector3[2];
	//get all points from spheres
	for(int is=0; is<4*App::nSectors; is++)
	{
		int phys_sec = normalizeSector(is);

		//printf("%d - > %d\n", is, phys_sec);
		
		int geom_sec = phys_sec;

		int nlayers = CGeometry::nRadiuses[geom_sec];

		//printf("R_out[%d] = %le(%d)\n", geom_sec, CGeometry::outer_radius[geom_sec][nlayers-1],nla7yers);

		double sect_left  = is*M_PI/(2*App::nSectors) - M_PI/(4*App::nSectors);
		double sect_right = is*M_PI/(2*App::nSectors) + M_PI/(4*App::nSectors);
		
		double left  = normalize_phi(sect_left);
		double right = normalize_phi(sect_right);
		
		
		for(int ilayer=0;ilayer<nlayers;ilayer++)
		{
			double Radius = CGeometry::outer_radius[geom_sec][ilayer];
			
			intersectSphere(ray,Radius,iPoint);

			
			for(int k=0;k<2;k++)
			{
				if(module(iPoint[k])<1.e+80)
				{

					double phi = getPhi(iPoint[k]);

					phi = normalize_phi(phi);
					
					if(checkAngleWithin(left, right, phi))
					{
						intersection curpoint;
						curpoint.point = iPoint[k];
						curpoint.sec_p = geom_sec;
						curpoint.sec_n = geom_sec;
						curpoint.lay_p = ilayer-1;
						curpoint.lay_n = ilayer;
						curpoint.is = is;
						curpoint.phi = getPhi(iPoint[k]);
						curpoint.R = calc_dist(iPoint[k],App::rayIntegration.start);
						CSolver::points[npoints] = curpoint;

						npoints++;
					}
					
				}
			}
		
			//delete[] iPoint;
		
		}
		
		vector3 iPointS = intersectRayWithPhi(ray,left);
		double is_phi = getPhi(iPointS);
		
		if(fabs(is_phi-left)<=1.e-12)
		{
			intersection curpoint;
			curpoint.point = iPointS;
			curpoint.sec_p = normalizeSector(is);
			curpoint.sec_n = normalizeSector(is-1);

			double pointR = module(iPointS);
			int clayer = -5;
			//printf("R: %le\n", pointR);
			for(int ilayer=0;ilayer<nlayers;ilayer++)
			{
				double Radius = CGeometry::outer_radius[geom_sec][ilayer];
				if(pointR<Radius)
				{
					clayer = ilayer-1;
					break;
				}
			}

			if(clayer>-5)
			{
				curpoint.lay_p = clayer;
				curpoint.lay_n = clayer;
				curpoint.R = calc_dist(iPointS, App::rayIntegration.start);
				CSolver::points[npoints] = curpoint;

				npoints++;
			}
		}
	}
	
	delete[] iPoint;

	/*if(fabs(App::rayToObj.angle.y + 0.6108652)<0.001)
	{
		printvec(ray.angle);
		printf("NP: %d\n", npoints);
		for(int ispt=0;ispt<npoints && ispt<10;ispt++)
		{
			printvec(CSolver::points[ispt].point);
			printf("Sec: %ld->%ld; Layer: %ld->%ld; %d\n", CSolver::points[ispt].sec_p,CSolver::points[ispt].sec_n,CSolver::points[ispt].lay_p,CSolver::points[ispt].lay_n,CSolver::points[ispt].is);
		}
	}*/
	//printf("NP: %d\n", npoints);
	CSolver::npoints = npoints;
	
	FILE *FP;

	//sorting
	int i,j;
	
	intersection temp;

	for(i=1;i<=npoints;i++)
	{
	 	for(j=0;j<npoints-i;j++)
	 	{
	 		if(CSolver::points[j].R>CSolver::points[j+1].R)
	 		{
	 			temp=CSolver::points[j];
	 			CSolver::points[j]=CSolver::points[j+1];
	 			CSolver::points[j+1]=temp;
	 		}
	 	}
	}

	/*if(fabs(App::rayIntegration.angle.y + 0.6108652)<0.001)
	{
		char fname[255];

		printvec(ray.start);
		printvec(ray.angle);

		sprintf(fname, "sorted_points_%lf.dat",App::rayIntegration.angle.z);
		FP = fopen(fname,"w+");
		for(i=0;i<npoints;i++)
		{
			fprintf(FP,"Sector: %d/%d (%d - %lf); Layer: %d/%d; Point: %le;%le;%le; R=%le\n", CSolver::points[i].sec_p, CSolver::points[i].sec_n,CSolver::points[i].is,CSolver::points[i].phi, CSolver::points[i].lay_p, CSolver::points[i].lay_n, CSolver::points[i].point.x, CSolver::points[i].point.y, CSolver::points[i].point.z,CSolver::points[i].R);
			//printvec(CSolver::points[i].point);
		}
		fclose(FP);
	}*/

	for(i=0; i<npoints-1;i++)
	{
		TPath cpath;
		if(i>=npoints-1 || (CSolver::points[i].lay_p<=CSolver::points[i+1].lay_p && CSolver::points[i].lay_p>-2))
		{
			cpath.layer = CSolver::points[i+1].lay_p;
			cpath.sector = CSolver::points[i+1].sec_p;
		}
		else
		{
			cpath.layer = CSolver::points[i].lay_p;
			cpath.sector = CSolver::points[i].sec_p;
		}
		
		
		cpath.path = calc_dist(CSolver::points[i].point, CSolver::points[i+1].point);
		CSolver::paths[i] = cpath;
	}

}