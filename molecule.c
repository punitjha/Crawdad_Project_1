#include "molecule.h"


/* ***********************************************************************************************************/	
//default constructor
	molecule::molecule (const char *filename, int q)
		{
        // open filename
  		ifstream myreadfile1 (filename);
		charge = q;
		if (myreadfile1.is_open())
		{
		cout<<"The file is now open"<<endl;
		// read the number of atoms from filename
        myreadfile1 >> natom;

		zvals = new int[natom];
  		geom = new double* [natom];
		for(int i=0; i < natom; i++)
    	geom[i] = new double[3];

  		for(unsigned int i=0; i < natom; i++)
    	myreadfile1 >> zvals[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];
		}
		else  
		cout << "Unable to open file"<<endl;
		
		}
		

		

/* ***********************************************************************************************************/	

//setter functions


	
	void molecule::print_geom()
	{
    for(int i=0; i < natom; i++)
    printf("%d %8.9f %8.9f %8.9f\n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);
	}
	
/* *****************************************************************/
    void molecule::rotate(double phi){
	
}

    void molecule::translate(double x, double y, double z)
	{
	for(int i=0; i < natom; i++) 
	 	{
     geom[i][0] += x;
     geom[i][1] += y;
     geom[i][2] += z;
     	}
	 }
/* *****************************************************************/
    
	double molecule::bond(int a, int b)
	{
	 return sqrt( (geom[a][0]*geom[b][0])*(geom[a][0]*geom[b][0])+(geom[a][1]*geom[b][1])*(geom[a][1]*geom[b][1])+(geom[a][2]*geom[b][2])*(geom[a][2]*geom[b][2]) );
	}



/* *****************************************************************/
    double molecule::angle(int atom1, int atom2, int atom3){
	
}

    double molecule::torsion(int atom1, int atom2, int atom3, int atom4){
	
}


/* ******************************************************************/	
//destructor
molecule::~molecule()
{
  delete[] zvals;
	for(int i=0; i < natom; i++)
    delete[] geom[i];
  	delete[] geom;
	for(int i=0; i < natom; i++) 
	delete[] R_bonds[i];
	delete[] R_bonds;
	
}

      
/* ***********************************************************************************************************/	
//general functions/display
