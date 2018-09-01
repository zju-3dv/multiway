


#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <vector>

using namespace std;// BOOOOOOOOOOOOOOOOOOOH!!!

int main()
{

//declaration of functions
	int BruteSearch(double* ,double* ,int,int,double* );
   void GenerateRandom(double *,int ,int );
   void BruteKSearch(double* ,double * ,int ,int ,int,double*,int*);
	void  BruteRSearch(double* ,double*  ,double ,int ,int,vector<int>* );
	void  BruteRSearchWithDistance(double* ,double*  ,double ,int ,int,vector<int>*,vector<double>* );
    vector<int>* BruteBoxSearch(double* ,double*,int,int);

	return 0;
}


 
int BruteSearch(double *p,double *qp ,int N,int dim,double*mindist)
{  //Ricerca bruta del nearest neighbour in 2D 
	//Calcola tutte le distanze e prende la minima
	// -future improvement> search with plane splitting filter
	

	 int idc,i,j;
	 
	 double dist;
    *mindist=HUGE_VAL;

     for (i=0;i<N;i++)
		{  
           dist=0;
		   for(j=0;j<dim && dist<*mindist;j++)
		   {dist=dist+(p[i+N*j]-qp[j])*(p[i+N*j]-qp[j]);}

		   if (dist<*mindist)
		   {*mindist=dist;
		   idc=i;}
	   }

	 return idc;
}	   



void BruteKSearch(double* p,double *qp ,int k,int N,int dim,double* distances,int * idc)
{  //Ricerca bruta del nearest neighbour in 2D 
	//Calcola tutte le distanze e prende la minima
	// -future improvement> search with plane splitting filter
	double dist;
		int i,j,count;
	  double mindist=HUGE_VAL;

	  for (i=0;i<k;i++)
	  { distances[i]=HUGE_VAL;}

     for (i=0;i<N;i++)
		{  
		   dist=0;
		   for(j=0;j<dim && dist<mindist;j++)
		   {dist=dist+(p[i+N*j]-qp[j])*(p[i+N*j]-qp[j]);}

		   if (dist<mindist)
		   {
			   count=0;
               for (j=1;j<k;j++)
			   {                  
				   if (dist<=distances[j])
						{count=count+1;}
                        else
						{break;}
			   }        
              for (j=0;j<count;j++)
			  {
                        idc[j]=idc[j+1];
                        distances[j]=distances[j+1];
			  } 
            idc[count]=i;
            distances[count]=dist;
            mindist=distances[0];
		   }
		   
	 }

}
	
			   




void GenerateRandom(double *p,int n,int dim) 
{//function to generate random points in 0-1 range
	//* pinter to vectorized array
	//n number of points
    // dimensions of points

// loop to generate random points not normalized
    int i;//counter;
	int tempmax=0;//temporary integer maximum random point;
	srand( (unsigned)time( NULL ) );
    
	

    for( i = 0;   i <n*dim;i++ )
	{		
		p[i]=rand();
		if (p[i]>tempmax)
		{tempmax=p[i];}
	}

	

    //loop to normalize
    for( i = 0;   i <n*dim;i++ )
	{
		p[i]=p[i]/tempmax;
	}
	
}





void BruteRSearch(double *p,double *qp ,double r,int N,int dim,vector<int>* idcv)
{  //Ricerca bruta del nearest neighbour in 2D 
	//Calcola tutte le distanze e prende la minima
	// -future improvement> search with plane splitting filter
	

	 int i,j;
	 double mindist=r*r;
	 double dist;


     for (i=0;i<N;i++)
		{  
           dist=0;
		   for(j=0;j<dim && dist<mindist;j++)
		   {
			   dist=dist+(p[i+N*j]-qp[j])*(p[i+N*j]-qp[j]);
		   }

		   if (dist<=mindist)
		   {
		   idcv->push_back(i);
		   }
	   }
     
}



vector<int>* BruteBoxSearch(double *p,double*box,int N,int dim)
{  //Ricerca bruta del nearest neighbour in 2D 
	//Calcola tutte le distanze e prende la minima
	// -future improvement> search with plane splitting filter
	

	 int i,j;
	vector<int> idbox;
	bool in;
    

     for (i=0;i<N;i++)
		{  
          in=true;
		   for(j=0;j<dim && in;j++)
		   {
			if (p[i+N*j]<box[1] && p[i+N*j]>box[2] && p[i+N*j]<box[3] && p[i+N*j]>box[4]);//checkfor outside points
		   {in=false;break;}
           }

		  if (in)  
		  {
			  idbox.push_back(i);
		  }
	   }

	 return &idbox;
}



void  BruteRSearchWithDistance(double *p,double *qp ,double r,int N,int dim,vector<int>* idcv,vector<double>* distvect)
{  //Ricerca bruta del nearest neighbour in 2D 
	//Calcola tutte le distanze e prende la minima
	// -future improvement> search with plane splitting filter
	

	 int i,j;
	 double mindist=r*r;
	 double dist;


     for (i=0;i<N;i++)
		{  
           dist=0;
		   for(j=0;j<dim && dist<mindist;j++)
		   {
			   dist=dist+(p[i+N*j]-qp[j])*(p[i+N*j]-qp[j]);
		   }

		   if (dist<=mindist)
		   {
		   idcv->push_back(i);
		   distvect->push_back(sqrt(dist));
		   }
	   }
     
}