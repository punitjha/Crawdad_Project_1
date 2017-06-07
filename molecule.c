#include "molecule.h"


/**************************************************************************************/	
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
				myreadfile1 >> zvals[i] >> geom[i][0] >> geom[i][1] >> 
					    geom[i][2];
		}
		else  
			cout << "Unable to open file"<<endl;
	}

/******************************************************************************************/	
//setter functions
/******************************************************************************************/	

	void molecule::print_geom()
	{
	    for(int i=0; i < natom; i++)
	    printf("%d %8.5f %8.5f %8.5f\n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);
	}
	
/******************************************************************/
	void molecule::print(double **mat, int row, int col)
	{
		for(int i=0; i<row; i++)
		{
			for(int j=0; j<col; j++)
			{
			printf("%8.5f "  "",mat[i][j]);
			}
		cout<<endl;
		}
	}
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
	 return sqrt( (geom[a][0]-geom[b][0])*(geom[a][0]-geom[b][0])+(geom[a][1]-geom[b][1])*
			 (geom[a][1]-geom[b][1])+(geom[a][2]-geom[b][2])*(geom[a][2]-geom[b][2]) );
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
		double theta= (((cross_x(i,j,k,l))*compo(0,k,i)+ cross_y(i,j,k,l)*compo(1,k,i)+
					cross_z(i,j,k,l)*compo(2,k,i))/sin(angle(j,k,l)));
		if (theta <-1.0)
			return asin(-1.0);
	 	else if (theta > 1.0)
			return asin(1.0);
		else return asin(theta);
	}
/****************************************************************/
	double molecule::torsion(int i, int j, int k, int l)
	{

//the cross product components

		double ijjk_x=(compo(1,i,j)*compo(2,j,k)-compo(2,i,j)*compo(1,j,k));
		double ijjk_y=(compo(2,i,j)*compo(0,j,k)-compo(0,i,j)*compo(2,j,k));
		double ijjk_z=(compo(0,i,j)*compo(1,j,k)-compo(1,i,j)*compo(0,j,k));
		double jkkl_x=(compo(1,j,k)*compo(2,k,l)-compo(2,j,k)*compo(1,k,l));
		double jkkl_y=(compo(2,j,k)*compo(0,k,l)-compo(0,j,k)*compo(2,k,l));
		double jkkl_z=(compo(0,j,k)*compo(1,k,l)-compo(1,j,k)*compo(0,k,l));
	
//the dot product components
		double dot_x=ijjk_x*jkkl_x;
		double dot_y=ijjk_y*jkkl_y;
		double dot_z=ijjk_z*jkkl_z;
		double theta=(dot_x+dot_y+dot_z)/(sin(angle(i,j,k))*sin(angle(j,k,l)));
	
//checking the numerical precision and the sign of torsion angle
		
		if (theta <-1.0)
			theta= acos(-1.0);
	 	else if (theta > 1.0)
			theta= acos(1.0);
		else theta= acos(theta);
		double cross_x=ijjk_y*jkkl_z-ijjk_z*jkkl_y;
		double cross_y=ijjk_z*jkkl_x-ijjk_x*jkkl_z;
		double cross_z=ijjk_x*jkkl_y-ijjk_y*jkkl_x;
		double norm=cross_x*cross_x+cross_y*cross_y+cross_z*cross_z;
		cross_x/=norm;
		cross_y/=norm;
		cross_z/=norm;
		double sign=1.0;
		double dot=cross_x*compo(0,j,k)+cross_y*compo(1,j,k)+cross_z*compo(2,j,k);
		if (dot < 0.0)
			sign=-1.0;
		return theta*sign;

	}
/****************************************************************/	
//destructor
/****************************************************************/	
	molecule::~molecule()
	{
		delete[] zvals;
		for(int i=0; i < natom; i++)
			delete[] geom[i];
		delete[] geom;
			
	}
/****************************************************************/	
