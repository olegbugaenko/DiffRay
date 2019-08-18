#pragma once
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

    vector3 operator-() const
    {
        return vector3(-x, -y, -z);
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

vector3 multiply(double M, vector3 v);

int sign(double val);

void printvec(vector3 a);

double module(vector3 a);

double getPhi(vector3 point);

double getTheta(vector3 point);

//a = (dr, dphi, dtheta);

vector3 spherical2desc(vector3 sph);

vector3 rotate(vector3 resv, vector3 a);

vector3 normalize(vector3 p);

vector3 rotateOverAxis(vector3 p, vector3 axis, double angle);


struct TRay
{
	vector3 start;
	vector3 angle;

	TRay(vector3 start = vector3(0,0,0), vector3 angle = vector3())
		: start(start), angle(angle)
		{

		}
};

TRay createRay(vector3 start, vector3 angle);

TRay rotateRay(TRay initial, vector3 angle);

vector3* intersectSphere(TRay ray, double radius, vector3 res[2]);


struct TPlane
{
	double A,B,C,D;

	TPlane(double A=0., double B=0., double C=0., double D=0.) 
        : A(A), B(B), C(C), D(D)
    {
    }
};

TPlane createPlaneFromPoints(vector3 M1, vector3 M2, vector3 M3);

TPlane getPlaneByAxisAndPhi(vector3 axis, double phi);

void printplane(TPlane P);

struct THalfplane
{
	TPlane    plane;
	vector3   lpoint;
	int side_des;
};

THalfplane create(TPlane P, vector3 point1, int side_des = 0);


THalfplane createFromAngles(double axis_phi, double axis_theta, double rot_phi);

vector3 intersectHalfplane(TRay ray, THalfplane sector);

//Closer to real geometry of object:
//Center of object always at zero
//So in fact half-plane is only defined by phi

TPlane makeCenteredByPhi(double phi);

vector3 intersectRayWithPhi(TRay ray, double phi);

double normalize_phi(double phi);

double calc_dist(vector3 a, vector3 b);

bool checkAngleWithin(double left, double right, double angle);