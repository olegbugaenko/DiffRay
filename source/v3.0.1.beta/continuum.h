#pragma once
#include "basics.h"
#include "reader.h"

class CContinuum{
	public:
	static int cellCount;
	static int nrows;
	static double ***emits;
	static double ***opacs;
	static double ***trans;
	static double *anu;
	static double endSkip;
	static int skipInterval;
	static int nsectors;
	static int iSector;
	//Init sectorial data
	static int initSectorials(int nsectors);
	//Read contmesh
	static int readMesh(char fname[500]);
	//Assign Mesh Size
	static int setMeshSize(int cells);
	//Put Mesh cell enerhy
	static int putMeshEnergy(int row, int col, double val);

	//Read cont. emissivity
	static int readEmissivity(char fname[500],int iSector);

	//Put cont. emissivity
	static int putContEmissivity(int row, int col, double val);

	//Read cont. opacity
	static int readOpacity(char fname[500],int iSector);
	
	//Put cont. opacity
	static int putContOpacity(int row, int col, double val);

	//Read transitions
	static int readTransitions(char fname[500],int iSector);
	
	//Put transitions
	static int putTransitions(int row, int col, double val);

	static int getRadii(int nrows);
};