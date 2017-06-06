#include "molecule.h"

/********************************************************/	
//default constructor
	molecule::molecule (const char *filename, int q)
	{
        // open filename
  		ifstream myreadfile1 (filename);
		charge = q;
		if (myreadfile1.is_open())
		{
			cout<<"The file is now open"<<endl;
			// read the number of atoms 
			myreadfile1 >> natom;
			zvals = new int[natom];
			geom = new double* [natom];
			for(int i=0; i < natom; i++)
				geom[i] = new double[2];
			for(unsigned int i=0; i < natom; ++i)
				myreadfile1 >> zvals[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];
		}
		else  
			cout << "Unable to open file"<<endl;
	}

/********************************************************/	

//setter functions

	void molecule::print_geom()
	{
	    for(int i=0; i < natom; i++)
	    printf("%d %8.11f %8.11f %8.11f\n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);
	}
	
/******************************************************************/
	void molecule::rotate(double phi)
	{
	}
/*****************************************************************/
	void molecule::print(double **mat, int row, int col)
	{
		for(int i=0; i<row; i++)
		{
			for(int j=0; j<col; j++)
			{
			printf("%8.9f "  "",mat[i][j]);
			}
		cout<<endl;
		}
	}
/*****************************************************************/

	/*****************************************************************/


	void molecule::translate(double x, double y, double z)
	{
		for(int i=0; i < natom; i++) 
		{
		     geom[i][0] += x;
		     geom[i][1] += y;
		     geom[i][2] += z;
		}
	}
/******************************************************************/
    
	double molecule::bond(int a, int b)
	{
	 return sqrt( (geom[a][0]-geom[b][0])*(geom[a][0]-geom[b][0])+(geom[a][1]-geom[b][1])*(geom[a][1]-geom[b][1])+(geom[a][2]-geom[b][2])*(geom[a][2]-geom[b][2]) );
	}
/******************************************************************/
	double molecule::compo( int cartesian,int i, int j)
	{
	return -(geom[i][cartesian]-geom[j][cartesian])/bond(i,j);
	}	
	
/******************************************************************/
	double molecule::angle(int i, int j, int k)
	{
		return acos(compo(0,j,i)*compo(0,j,k)+compo(1,j,i)*compo(1,j,k)+compo(2,j,i)*compo(2,j,k));
	}
/****************************************************************/
	double molecule::cross_x(int i, int j, int k, int l )
	{
		return ( (compo(1,k,j)*compo(2,k,l)-compo(2,k,j)*compo(1,k,l)));
	}
/****************************************************************/
	double molecule::cross_y(int i, int j, int k, int l )
	{
		return ( (compo(2,k,j)*compo(0,k,l)-compo(0,k,j)*compo(2,k,l)));
	}

/****************************************************************/
	double molecule::cross_z(int i, int j, int k, int l )
	{
		return ( (compo(0,k,j)*compo(1,k,l)-compo(1,k,j)*compo(0,k,l)));
	}



/****************************************************************/

	double molecule::outplane(int i, int j, int k, int l)
	{
		double theta= (((cross_x(i,j,k,l))*compo(0,k,i)+ cross_y(i,j,k,l)*compo(1,k,i)+cross_z(i,j,k,l)*compo(2,k,i))/sin(angle(j,k,l)));
		if (theta <-1.0)
			return asin(-1.0);
	 	else if (theta > 1.0)
			return asin(1.0);
		else return asin(theta);
	}



/****************************************************************/


	double molecule::torsion(int atom1, int atom2, int atom3, int atom4)
	{
		return acos(cross_x(l,j,i,k)*cross()) 
	}

/*******************************************************************/	
//destructor
	molecule::~molecule()
	{
		delete[] zvals;
		for(int i=0; i < natom; i++)
			delete[] geom[i];
		delete[] geom;
			
	}

      
/*******************************************************************************/	
//general functions/display
