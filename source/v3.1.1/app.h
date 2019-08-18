#pragma once
#include "basics.h"
#include "geom3d.h"

class App{
public:
	static char *model_name;
	static char *data_suffix;
	static int nSectors;
	static vector3 obj_rotation;
	static double distance;
	static double phi;
	static double theta;
	static double phi_width;
	static double theta_width;
	static double cov_fac;
	static bool CalcCont;
	static bool CalcLines;
	static bool CalcOpac;
	static bool CalcAbund;
	static double age;
	static char output_dir[255];
	static char input_dir[255];

	static TRay rayToObj; //ray to center
	static TRay rayIntegration;
	static bool init();
	static bool setRay(double phi, double theta);
	static bool readCommands();
};