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
	static bool CalcIRBands;
	static double age;
	static char output_dir[255];
	static char output_dir_in[255];
	static char input_dir[255];

	static TRay rayToObj; //ray to center
	static TRay rayIntegration;
	static int AppMode;
	static bool isStatMode;

	static Apperture appertures[10];
	static int nApps;
	static int iApp;
	static double AppDPhi;
	static double AppDTheta;

	static bool init();
	static bool setRay(double phi, double theta);
	static bool readCommands();
	static bool addApperture(double phi, double theta, double dphi, double dtheta);
	static bool initAperture(int i);
};