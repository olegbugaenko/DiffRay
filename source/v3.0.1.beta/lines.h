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
	static int initSectorials(int nsectors);
	static int readLines(char fname[500],int iSector);
	static int putLine(int nlines);
	static int putLineVal(int raw, int col, double val);
	static int getRadii(int nrows);
};