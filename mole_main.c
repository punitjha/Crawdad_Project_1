#include <iostream> //
#include <fstream> // for r/w to a file
#include <string>
#include <sstream>
#include <cstring>
#include <vector> //for using the vectors in the program
#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include "molecule.h"
using namespace std;

/* ***********************************************************************************************************/	
int main()
{
	molecule acet("acetal.txt", 0);
	acet.print_geom();
	R_bonds = new double* [natom];
	for(int i=0; i < natom; i++)
	R_bonds[i] = new double[natom];
	for (int i=0; i< natom; i++)
	for (int j=1; j< natom; j++)
	{
		R_bonds[i][j]=acet.bond(i,j);
		cout<<R_bonds[i][j]<<endl;
	}

}
