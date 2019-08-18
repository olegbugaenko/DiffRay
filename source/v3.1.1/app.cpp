#include "basics.h"
#include "app.h"
#include "geom3d.h"
#include "math.h"
#include "reader.h"


char *App::model_name = "wew";
char *App::data_suffix = "Age60.00Myr.dat";
int App::nSectors = 20; //set to 20
vector3 App::obj_rotation = vector3(
		0, //that means nothing
		0,//M_PI/4, //Ange between sector 0 midline and direction to observer
		0  //Angle between plane, perpendicular to symmethry axis and direction to observer 
	);

double App::distance = 1.e+23;
double App::phi = 0.*M_PI/400;
double App::theta = 0.*M_PI/6;
double App::phi_width = 2*M_PI;
double App::theta_width = M_PI;
double App::cov_fac = 0.025;
double App::age = 60.0;
char App::output_dir[255];
char App::input_dir[255];

bool App::CalcCont = true;
bool App::CalcOpac = true;
bool App::CalcLines = false;
bool App::CalcAbund = false; 

TRay App::rayToObj;
TRay App::rayIntegration;

bool App::init()
{

	vector3 start_pos = spherical2desc(vector3(App::distance, M_PI + App::phi, -App::theta));

	vector3 angle = vector3(0, M_PI + App::phi, App::theta);

	TRay rayToObj = TRay(start_pos, vector3(0,App::phi,App::theta));

	rayToObj = rotateRay(rayToObj,multiply(1.0,App::obj_rotation));

	App::rayToObj = rayToObj;

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
		}
		else
		if(strcmp(commands[0],"print") == 0)
		{
			printf("print\n");
			if(strcmp(commands[1],"directory") == 0)
			{
				sscanf(commands[2],"%s",&App::output_dir);
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
