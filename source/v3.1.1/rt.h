#pragma once
#include "geom3d.h"
#include "solver.h"

class CRT 
{
public:
	static double *ray_lines;
	static double *ray_cont;
	static double dS;

	static int calc_ray(int mode=0);
	static int calc_cont(int mode=0);
	static int get_mesh_ind(double energy);
};