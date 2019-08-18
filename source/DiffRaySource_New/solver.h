#pragma once;
#include "basics.h"
#include "geom3d.h"

struct intersection
{
	vector3 point;
	int sec_p,sec_n;
	int lay_p,lay_n;
	double R;

	int is;
	double phi;
};

struct  TPath
{
	int layer;
	int sector;
	double path;
};

class CSolver
{
public:
	static intersection points[3000];
	static TPath paths[3000];
	static int npoints;

	static int getPoints();
};