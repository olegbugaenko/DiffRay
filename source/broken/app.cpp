#include "basics.h"
#include "app.h"
#include "geom3d.h"
#include "math.h"
#include "reader.h"
#include <dirent.h>
#include <sys/stat.h>


char *App::model_name = "wew";
char *App::data_suffix = "Age60.00Myr.dat";
int App::nSectors = 20; //set to 20
vector3 App::obj_rotation = vector3(
		0, //that means nothing
		//0.225*M_PI,//0.376094213*M_PI,//M_PI/4, //Ange between sector 0 midline and direction to observer
		(12.5/40)*M_PI,
		0//0.360278922*M_PI  //Angle between plane, perpendicular to symmethry axis and direction to observer 
	);

double App::distance = 1.e+23;
double App::phi = 0;//0.376094213*M_PI;
double App::theta = 0;//0.360278922*M_PI;
double App::phi_width = 2*M_PI;
double App::theta_width = M_PI;
double App::cov_fac = 0.025;
double App::age = 60.0;
double App::AppDPhi = 0.0;
double App::AppDTheta = 0.0;

char App::output_dir[255];
char App::output_dir_in[255];
char App::input_dir[255];

bool App::CalcCont = true;
bool App::CalcOpac = false;
bool App::CalcLines = true;
bool App::CalcAbund = true; 
bool App::CalcIRBands = true;
bool App::isStatMode = true;

TRay App::rayToObj;
TRay App::rayIntegration;

int App::AppMode = 0;

int App::iApp = 0;

Apperture App::appertures[10];
int App::nApps = 0;

bool App::addApperture(double phi, double theta, double dphi, double dtheta){
	
	App::AppMode = 1;
	
	App::appertures[App::nApps].reset(phi,theta,dphi,dtheta);
	
	App::nApps++;
}

bool App::init()
{

	vector3 start_pos = spherical2desc(vector3(App::distance, M_PI + App::phi, -App::theta));

	vector3 angle = vector3(0, M_PI + App::phi, App::theta);

	TRay rayToObj = TRay(start_pos, vector3(0,App::phi,App::theta));

	rayToObj = rotateRay(rayToObj,multiply(1.0,App::obj_rotation));

	App::rayToObj = rayToObj;

}

bool App::initAperture(int i)
{
	App::iApp = i;
	Apperture a = App::appertures[i];
	App::AppDPhi = a.phi;
	App::AppDTheta = a.theta;
	App::phi_width = a.phi_width;
	App::theta_width = a.theta_width;

	//char path[500];

	sprintf(App::output_dir,"%s/app_%d",App::output_dir_in,i);

	printf("App to %s\n",App::output_dir);

	DIR *st;

	st = opendir(App::output_dir);
	if (st == NULL || !st)
    {
        /* Directory does not exist. EEXIST for race condition */
        mkdir(App::output_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

	App::init();
}

bool App::setRay(double phi, double theta)
{
	vector3 start_pos = spherical2desc(vector3(App::distance, M_PI + phi, -theta));

	vector3 angle = vector3(0, M_PI + phi, theta);

	TRay rayToObj = TRay(start_pos, vector3(0, phi, theta));

	rayToObj = rotateRay(rayToObj,multiply(1.0,App::obj_rotation));

	/*printf("RAY SET TO:\n");
	printf("--- %le; %le -----\n", phi, theta);
	printvec(rayToObj.start);
	printvec(rayToObj.angle);*/

	App::rayToObj = rayToObj;
}

bool App::readCommands()
{
	FILE *F;

	F = fopen("commands.ini","r");
	char **commands;
	
	commands = CReader::read_commands(F);
	while(commands != NULL && commands[0] != NULL)
	{
		printf("------------%s------------------\n",commands[0]);
		if(strcmp(commands[0],"distance") == 0)
		{
			double factor = 1.0;

			printf("distance\n");
			if(strcmp(commands[1],"meter") == 0)
			{
				factor = 100;
			}
			else
			if(strcmp(commands[1],"km") == 0)
			{
				factor = 1.e+5;
			}

			double distance;

			printf("%s\n", commands[2]);

			sscanf(commands[2],"%lf",&distance);

			printf("DIST: %le",distance);
			
			App::distance = pow(10.0,distance)*factor;
		}
		else
		if(strcmp(commands[0],"object") == 0)
		{
			printf("object\n");
			if(strcmp(commands[1],"phi") == 0)
			{
				sscanf(commands[2],"%lf",&App::phi);
			}
			else
			if(strcmp(commands[1],"theta") == 0)
			{
				sscanf(commands[2],"%lf",&App::theta);
			}
			else
			if(strcmp(commands[1],"age") == 0)
			{
				sscanf(commands[2],"%lf",&App::age);
			}		
		}
		else
		if(strcmp(commands[0],"apperture") == 0)
		{
			printf("aperture\n");
			if(strcmp(commands[1],"width") == 0)
			{
				if(strcmp(commands[2],"full") == 0)
				{
					App::phi_width = 2.0*M_PI;
					App::theta_width = M_PI;
				}
			}
			else if(strcmp(commands[1],"add") == 0)
			{
				printf("Apperture added\n");
				double phi, theta, phi_width, theta_width;
				sscanf(commands[2],"%lf",&phi);
				sscanf(commands[3],"%lf",&theta);
				sscanf(commands[4],"%lf",&phi_width);
				sscanf(commands[5],"%lf",&theta_width);
				App::addApperture(phi,theta,phi_width,theta_width);
			}
		}
		else
		if(strcmp(commands[0],"print") == 0)
		{
			printf("print\n");
			if(strcmp(commands[1],"directory") == 0)
			{
				sscanf(commands[2],"%s",&App::output_dir);
				sscanf(commands[2],"%s",&App::output_dir_in);
			}
		}
		else
		if(strcmp(commands[0],"input") == 0)
		{
			printf("input\n");
			if(strcmp(commands[1],"location") == 0 || strcmp(commands[1],"directory") == 0) //for compability with DiffRay v2.0
			{
				sscanf(commands[2],"%s",&App::input_dir);
			}
			printf("DIR: %s\n",App::input_dir);
		}
		
		commands = CReader::read_commands(F);
	};
	
	fclose(F);

	App::setRay(App::phi, App::theta);

	return true;
}
