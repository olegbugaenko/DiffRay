#pragma once
#include "basics.h"
#include "reader.h"

class CLine{
	public:
	static int linesCount;
	static int nrows;
	static double ***emits;
	static int nsectors;
	static int iSector;
	static int lineIds[100];
	static int initSectorials(int nsectors);
	static int readLines(char fname[500],int iSector);
	static int putLine(int nlines);
	static int putLineVal(int raw, int col, double val);
	static int getRadii(int nrows);
	static int getLineList(int irow,char str[1000]);
	static int readDatabase();

	static char linesCapDB[100000][14];
	static double *linesEn;
	static int ndatabase;
	static int readen;
};