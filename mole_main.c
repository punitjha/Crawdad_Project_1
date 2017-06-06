#include "molecule.h"
#include "mass.h"
using namespace std;

/**************************************************************/	
int main()
{
	molecule acet("acetal.txt", 0);
	cout<<"We now print the geometry"<<endl;
	acet.print_geom();
	double **R_bonds=new double* [acet.natom];
	for(int i=0; i<acet.natom; i++)
	{
		R_bonds[i]=new double [acet.natom];
	}
	for(int i=0; i<acet.natom; i++)
	{
		for(int j=0; j<acet.natom; j++)
		{
			R_bonds[i][j]=acet.bond(i,j);
			//cout<<R_bonds[i][j]<<"   ";
		}
	}
	cout<<endl<<"Printing the bond distances"<<endl;
	acet.print(R_bonds, acet.natom,acet.natom);
	double **ex=new double* [acet.natom];
	double **ey=new double* [acet.natom];
	double **ez=new double* [acet.natom];
	for(int i=0; i< acet.natom;i++)
	{
		ex[i]=new double[acet.natom];
		ey[i]=new double[acet.natom];
		ez[i]=new double[acet.natom];
	}
	for( int i=0; i<acet.natom; i++)
	{
		for(int j=0; j<i ;j++)
		{
			ex[i][j]=ex[j][i]=-(acet.geom[i][0]-acet.geom[j][0])/R_bonds[i][j];
			ey[i][j]=ey[j][i]=-(acet.geom[i][1]-acet.geom[j][1])/R_bonds[i][j];
			ez[i][j]=ez[j][i]=-(acet.geom[i][2]-acet.geom[j][2])/R_bonds[i][j];
		}	
	}
	int arlen=acet.natom;
	//cout<<endl<< " We now print the x-components "<<endl;
	//acet.print(ex,arlen,arlen);

//three dimensional array to store the bond angles
	double ***phi = new double** [acet.natom];
	for(int i=0; i< acet.natom; i++)
	{
		phi[i]=new double* [acet.natom];
		for(int j=0; j<acet.natom; j++)
		{
			phi[i][j]=new double[acet.natom];
		}
	}
	cout<<endl<<"Printing the bond angles"<<endl;
	for(int i=0; i< acet.natom; i++)
	{
		for (int j=0; j<i; j++)
		{
			for(int k=0; k<j; k++)
			{
				if (R_bonds[i][j]<4.0 && R_bonds[j][k] < 4.0)
				{
				     phi[i][j][k]=acet.angle(i,j,k);
					printf("%2d %2d %2d %10.4f\n",i,j,k,phi[i][j][k]*(180.0/acos(-1.0)));
												}
											}
										}
									}
		
	
//now we calculate the out of plane angles
	cout<<endl<<"Out of plane angles"<<endl;
	for (int i=0; i<acet.natom; i++)
	{
		for (int k=0; k<acet.natom; k++)
		{
			for(int j=0; j<acet.natom; j++)
			{
				for(int l=0; l<j; l++)
				{
					if (R_bonds[i][k] < 4.0 && R_bonds[j][k]< 4.0 && R_bonds [k][l]< 4.0 && i!=j && i!=k && i!= l && j!=k && k!=l ) 
					{ 
						printf(" %2d %2d %2d %2d %10.5f\n", i,j,k,l,acet.outplane(i,j,k,l)*(180.0/acos(-1.0)));
					}
				  }
			 }
		}
	}
//now we calculate the torsion/dihedral angles
	cout<<endl<<"The torsion angles"<<endl;
	for (int i=0; i<acet.natom; i++)
	{
		for (int j=0; j<i; j++)
		{
			for(int k=0; k<j; k++)
			{
				for(int l=0; l<k; l++)
				{
					if (R_bonds[i][j] < 4.0 && R_bonds[j][k]< 4.0 && R_bonds [k][l]< 4.0 ) 
					{ 
						printf(" %2d %2d %2d %2d %10.5f\n", i,j,k,l,acet.torsion(i,j,k,l)*(180.0/acos(-1.0)));
					}
				  }
			 }
		}
	}

//center of mass translation
	double x_cm=0.0;
	double y_cm=0.0;
	double z_cm=0.0;
	double sum=0.0;
	for (int i=0; i<acet.natom; i++)
	{
		x_cm=x_cm+masses[acet.zvals[i]]*acet.geom[i][0];
		y_cm=y_cm+masses[acet.zvals[i]]*acet.geom[i][1];
		z_cm=z_cm+masses[acet.zvals[i]]*acet.geom[i][2];
		sum=sum+masses[acet.zvals[i]];
	}
	cout<<endl<<"X_COM"<<' '<<x_cm/sum<<endl;
	cout<<"Y_COM"<<' '<<y_cm/sum<<endl;
	cout<<"Z_COM"<<' '<<z_cm/sum<<endl;

//translating the molecule to the centre of mass
	acet.translate(x_cm/sum, y_cm/sum, z_cm/sum);
	cout<<endl<<"Translating it to the centre of mass"<<endl;
	acet.print_geom();
	



//Calculating the principle moments of inertia

	double



//deleting the memory allocated to the bond vector
	for(int i=0; i < acet.natom; i++) 
		delete[] R_bonds[i];
	delete[] R_bonds;
//deleting the memory allocated to the angle vector
	for(int i=0; i<acet.natom; i++)
	{
		delete[] ex[i];
		delete[] ey[i];
		delete[] ez[i];
	}
	delete[] ex;
	delete[] ey;
	delete[] ez;
//deleting the memory allocated to the phi vector
	for (int i=0; i<acet.natom; i++)
	{
	}
















}
