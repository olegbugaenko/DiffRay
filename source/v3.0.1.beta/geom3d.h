#pragma once;
#include "basics.h"

struct vector3
{
	double x;
	double y;
	double z;

	vector3(double x=0., double y=0., double z=0.) 
        : x(x), y(y), z(z)
    {
    }

	vector3& operator=(const vector3& a)
    {
        x=a.x;
        y=a.y;
        z=a.z;
        return *this;
    }

    // addop. doesn't modify object. therefore const.
    vector3 operator+(const vector3& a) const
    {
        return vector3(a.x+x, a.y+y, a.z+z);
    }

    vector3 operator-(const vector3& a) const
    {
        return vector3(x-a.x, y-a.y, z-a.z);
    }

    double operator*(const vector3 a)
    {
    	return a.x*x+a.y*y+a.z*z;
    }

    // equality comparison. doesn't modify object. therefore const.
    bool operator==(const vector3& a) const
    {
        return (x == a.x && y == a.y && z == a.z);
    }

    

};

int sign(double val)
{
	if(val>=0.0)
		return 1;
	else
		return -1;
}

void printvec(vector3 a)
{
	printf("[%le,%le,%le]\n", a.x, a.y, a.z);
}

double module(vector3 a)
{
	return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
}

double getPhi(vector3 point)
{
	double dxy = sqrt(point.x*point.x+point.y*point.y);
	
	if(dxy<1.e-40*point.x)
		return 0;

	if(point.y>0)
		return acos(point.x/dxy);
	else
		return 2.0*M_PI -  acos(point.x/dxy);
}

double getTheta(vector3 point)
{
	double dxy = module(point);
	if(dxy>1.e-40*fabs(point.z))
		return asin(point.z/dxy);
	else
		return M_PI*sign(point.z)/2.0;
}

//a = (dr, dphi, dtheta)
vector3 rotate(vector3 resv, vector3 a)
{

	double r = module(resv);
	double phia = getPhi(resv);
	double theta = getTheta(resv);

	printf("PHI: %lf\n", phia);
	printf("THETA: %lf\n", theta);
	printf("R: %lf\n", r);
	printvec(a);
	printvec(resv);

	vector3 rsv;
	rsv.x = (r+a.x)*cos(phia+a.y)*cos(theta+a.z);
	rsv.y = (r+a.x)*sin(phia+a.y)*cos(theta+a.z);
	rsv.z = (r+a.x)*sin(theta+a.z);
	
	return rsv;
}

struct TRay
{
	vector3 start;
	vector3 angle;
};

TRay createRay(vector3 start, vector3 angle)
{
	TRay result;
	result.start = start;
	result.angle = angle;
	return result;
}

vector3* intersectSphere(TRay ray, double radius)
{
	vector3 *res;
	res = new vector3[2];
	//square root solver

	double b = 2*(ray.start.x*cos(ray.angle.y)*cos(ray.angle.z)+ray.start.y*sin(ray.angle.y)*cos(ray.angle.z)+ray.start.z*sin(ray.angle.z));
	double c = pow(ray.start.x,2)+pow(ray.start.y,2)+pow(ray.start.z,2)-radius*radius;	

	double D = b*b-4*c;

	res[0] = vector3(1.e+80,1.e+80,1.e+80);
	res[1] = vector3(1.e+80,1.e+80,1.e+80);
	
	
	if(D>=0.0)
	{
		double r1 = (-b + sqrt(D))/(2.0);
		double r2 = (-b - sqrt(D))/(2.0);
		
	
		if(r1>0)
		{
			res[0] = vector3(ray.start.x+r1*cos(ray.angle.y)*cos(ray.angle.z),ray.start.y+r1*sin(ray.angle.y)*cos(ray.angle.z),ray.start.z+r1*sin(ray.angle.z));
		}

		if(r2>0)
		{
			res[1] = vector3(ray.start.x+r2*cos(ray.angle.y)*cos(ray.angle.z),ray.start.y+r2*sin(ray.angle.y)*cos(ray.angle.z),ray.start.z+r2*sin(ray.angle.z));
		}
	}
	return res;
}

struct TPlane
{
	double A;
	double B;
	double C;
	double D;
};

vector3 getPointOnLine(vector3 line, vector3 start, double par)
{
	vector3 res;
	res.x = start.x+line.x*par;
	res.y = start.y+line.y*par;
	res.z = start.z+line.z*par;
	return res;
}

vector3 normalize(vector3 p)
{
	double dis = module(p);
	return vector3(p.x/dis,p.y/dis,p.z/dis);
}

vector3 rotateOverAxis(vector3 p, vector3 a, double angle)
{
	//Translator matrix should apply to normalized only axis
	vector3 axis = normalize(a);

	vector3 result;

	double cosin = cos(angle);
	double m_cosin = 1 - cos(angle);
	double sinus = sin(angle);

	result.x = (cosin+m_cosin*axis.x*axis.x)*p.x        +(m_cosin*axis.x*axis.y-sinus*axis.z)*p.y         +(m_cosin*axis.x*axis.z+sinus*axis.y)*p.z;
	result.y = (m_cosin*axis.x*axis.y+sinus*axis.z)*p.x   +(cosin+m_cosin*axis.y*axis.y)*p.y                +(m_cosin*axis.y*axis.z-sinus*axis.x)*p.z;
	result.z = (m_cosin*axis.z*axis.x-sinus*axis.y)*p.x +(m_cosin*axis.z*axis.y+sinus*axis.x)*p.y         +(cosin+m_cosin*axis.z*axis.z)*p.z;

	return result;
}

TPlane getFromPoints(vector3 m1, vector3 m2, vector3 m3)
{
	TPlane result;

	result.A = ((m2.y-m1.y)*(m3.z-m1.z) - (m3.y-m1.y)*(m2.z-m1.z));
	result.B = ((m2.z-m1.z)*(m3.x-m1.x) - (m3.z-m1.z)*(m2.x-m1.x));
	result.C = ((m2.x-m1.x)*(m3.y-m1.y) - (m3.x-m1.x)*(m2.y-m1.y));
	result.D = -1.*(result.A*m1.x + result.B*m1.y+result.C*m1.z);

	return result;
}

void printplane(TPlane a)
{
	printf("[%le*x+%le*y+z*%le+%le]\n", a.A, a.B, a.C, a.D);
}

TPlane getFromPhi(vector3 axis, double phi)
{
	vector3 pointA,pointB,pointC;

	double phic = getPhi(axis);
	double theta = getTheta(axis);

	pointA = getPointOnLine(axis, vector3(0.,0.,0.),0.0);
	pointB = getPointOnLine(axis, vector3(0.,0.,0.),10.0);

	double theta_N = theta;

	if(theta>M_PI/4.0)
		theta_N = -M_PI/6.0;
	else
		theta_N =  M_PI/6.0;

	vector3 rot_point = vector3(0.0,0.0,theta_N);

	if(module(pointB)>0)
		pointC = rotate(pointB,rot_point);
	else
		pointC = rotate(pointA,rot_point);

	
	printf("AXIS:\n");
	printvec(axis);
	printvec(pointC);

	//TPlane test = getFromPoints(pointA,pointB,pointC);
	printf("TEST:\n");
	//printplane(test);

	//now point c should be rotated over axis pointA - pointB
	pointC = rotateOverAxis(pointC,axis,phi);
	printvec(pointC);

	return getFromPoints(pointA,pointB,pointC);
}

