//Header
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <vector>
#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
/***********************************************************************************************************/	
using namespace std;

class molecule
	{
	public:
		int natom;
		int charge;
		int *zvals;
		double **geom;
		string point_group;


//the functions are defined here
		void print_geom();
		void rotate(double phi);
		void translate(double x, double y, double z);
		double bond(int atom1, int atom2);
		double angle(int i, int j, int k);
		double compo(int i, int j, int cartesian);
		double cross_x(int i, int j, int k, int l);
		double cross_y(int i, int j, int k, int l);
		double cross_z(int i, int j, int k, int l);
		double outplane( int i, int j, int k, int l);
		double torsion(int atom1, int atom2, int atom3, int atom4);
		void print(double **mat, int row, int col);	 
/***************************************************************************************************************/
//default constructor
		molecule();
		molecule (const char *filename, int q);


/************************************************************************************************************/	
 //destructor
		~molecule();

/***********************************************************************************************************/	


	private:
//Member variables
	

};
