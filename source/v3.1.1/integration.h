#include "basics.h"
#include "geom3d.h"
#include "solver.h"
#include "rt.h"
#include "geometry.h"
#include "continuum.h"
#include "lines.h"

class CIntegration
{
public:
	static int mode;
	static double *flux_lines;
	static double *flux_continua;
	static int nPhi,nTheta;
	static double lumfactor;
	static int InitIntegration(int mode=0,int nPrecision=10);
	static int doCalc();
	static int finish();
	static long long memPerGeom;
	static long long memPerRT;
	
};