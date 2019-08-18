#include "basics.h"
int** CBasics::int2Arr(int **arr,int sizex, int sizey, int defaultval)
{
	arr = new int*[sizex];
	for(int x=0;x<sizex;x++)
	{
		arr[x] = new int[sizey];
		for(int y=0;y<sizey;y++)
		{
			arr[x][y] = defaultval;
		}
	}
	
	//for(int i=0;i<4;i++)
	//{
	//	printf("%d:%d -> %d\n",i,i,arr[i][i]);
	//}

	return arr; 
};

double** CBasics::float2Arr(double **arr,int sizex, int sizey, double defaultval)
{
	arr = new double*[sizex];
	for(int x=0;x<sizex;x++)
	{
		arr[x] = new double[sizey];
		for(int y=0;y<sizey;y++)
		{
			arr[x][y] = defaultval;
		}
	}
	
	//for(int i=0;i<4;i++)
	//{
	//	printf("%d:%d -> %d\n",i,i,arr[i][i]);
	//}

	return arr; 
};

int CBasics::throwError(const char *string)
{
	printf("%s\n", string);
	exit(1);
}