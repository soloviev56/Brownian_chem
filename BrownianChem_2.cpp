/* Brownian dynamics simulation for the system of different types particles with chemical reactions  */
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cmath>
# include <fstream> 
#include <random>
# include <string>



using namespace std;

#define PI 3.141592653589793

void timestamp ( );
double cpu_time ( );


  
void funInterac0(double parameters[], double unidelRij[3], double sqRij, double Fij[], double &Eij );
void funInterac1(double parameters[], double unidelRij[3], double sqRij, double Fij[], double &Eij );

typedef void funInterac(double parameters[], double delRij[3], double sqRij, double Fij[], double &Eij  );
funInterac *funInteracArr[] = { funInterac0, funInterac1 };


int funReact0(double parameters[], double sqRij, double Eij);
int funReact1(double parameters[], double sqRij, double Eij);

typedef int funReact(double parameters[], double sqRij, double Eij);
funReact *funReactArr[] = { funReact0, funReact1 };



int main ( int argc, char *argv[] )
{
  double ctime, tmp1, tmp2, tmp3;
  char* modelFilename;
  unsigned int particklTipesNum;
  unsigned int reactionsNumber;
  unsigned int totalParticklNumber,existed;
  double delRij[3],delRij2,sqRij,xi,Eij,hFi[3], totEnergy, totSqrij;
  double *Fij;
  double rand1,rand2,kBT,kBTS,randhBr[3];
  int reactionResult=0;

  unsigned int id1,id2;
	const double kB=1.38064852E-09; // [fg A^2  /(ps^2 )/K]




  struct partickleTipe
  {
    char name[20];
    unsigned int partTypeID;
    unsigned int maxNumber;
    unsigned int initNumber;
    double diameter;
    double viscosityPar;
    double randomForceScale;

  };
  const int maxParametersNumber=10;
  double pars[maxParametersNumber];

  struct interactParams
  {
    unsigned int funcInteractIndex;
    unsigned int funcParametersNumber;
    double funcParameters[maxParametersNumber];

  };

  struct reactionsParams
  {
    unsigned int reactionTipe;
    unsigned int reagentsNumber;
    unsigned int productsNumber;
    unsigned int reagents[3];
    unsigned int products[3];
    unsigned int funcReactionIndex;
    unsigned int funcReactParametersNumber;
    double funcReactionParams[maxParametersNumber];
  
  };

  struct particleArr
  {
    unsigned int partTypeID;
    bool existence;
    double coord[3];
    double force[3];
    double delBrouncoord[3];
  };


  timestamp ( );

// Start timing  
  ctime = cpu_time ( );


//Read model parameters file

  if ( 1 == argc )
  {
    cout << "No model file name is given" <<endl;
    exit(1);
  }
  else
  {
    modelFilename = argv[1];
  }

  ifstream modelFile( modelFilename );
  if(!modelFile)
	{
		cout<<"Error in opening model file\n";
		exit(1);
	}
  cout << "Model file name is: " << modelFilename << "\n" <<"\n" << endl;
//Read partickle tipes
  modelFile>>particklTipesNum;
  partickleTipe particklTipes[particklTipesNum];
  unsigned int numSqij[particklTipesNum][particklTipesNum],nS;
  double avgSqij[particklTipesNum][particklTipesNum];
  double avgEij[particklTipesNum][particklTipesNum];

  for(int i=0; i < particklTipesNum; i++)
  {
    modelFile>>
      particklTipes[i].name>>
      particklTipes[i].partTypeID>>
      particklTipes[i].maxNumber>>
      particklTipes[i].initNumber>>
      particklTipes[i].diameter>>
      particklTipes[i].viscosityPar>>
      particklTipes[i].randomForceScale;
  }
// Read interactions parameters
  //interactParams particklParams[particklTipesNum+particklTipesNum*(particklTipesNum-1)/2];
  interactParams particklParams[particklTipesNum][particklTipesNum];
  int l=0;
  for(int i=0; i<particklTipesNum; i++)
    {
      for(int j=0; j<= i; j++)
      {
        modelFile>>id1>>id2>>
          particklParams[id1][id2].funcInteractIndex>>
          particklParams[id1][id2].funcParametersNumber;

        for (int k=0; k<particklParams[id1][id2].funcParametersNumber; k++)	
	      {
	        modelFile >>  particklParams[id1][id2].funcParameters[k];
	      }
        l++;
      }
    }

//Read reactions
  modelFile>>reactionsNumber;
  reactionsParams reactions[reactionsNumber];
  

  for(int i=0; i<reactionsNumber; i++)
  {
    modelFile>>
      reactions[i].reactionTipe>>
      reactions[i].reagentsNumber>>
      reactions[i].productsNumber;
      
	    for (int k=0; k<reactions[i].reagentsNumber; k++)	
	    {
	        modelFile >>  reactions[i].reagents[k];
	    }
	    for (int k=0; k<reactions[i].productsNumber; k++)	
	    {
	        modelFile >>  reactions[i].products[k];
	    }
      modelFile>>
        reactions[i].funcReactionIndex >>
        reactions[i].funcReactParametersNumber;
      for (int k=0; k<reactions[i].funcReactParametersNumber; k++)	
	      {
	        modelFile >>  reactions[i].funcReactionParams[k];
	      }

  }
//Read general parameters
  double T,eta;
  double delTime;
  int timeStepsNumber;
  int dataCollectNumber;
  int relaxStepsNumber;
  double boxSize[3];
  bool periodic;
  string xyzFileName,xyzFileCoords;
  string csvFileName;

//Time step, Number of steps, Data collection number of steps, equilibration
	modelFile>>T>>eta;
  modelFile>>delTime>>timeStepsNumber>>dataCollectNumber>>relaxStepsNumber; 
  modelFile>>periodic>>boxSize[0]>>boxSize[1]>>boxSize[2]; //Size of the box
  modelFile>>xyzFileName>>csvFileName;
  ofstream xyzFile( xyzFileName ); 
  ofstream csvFile( csvFileName );
  kBT=kB*T;    
/////////////////////////////////////////////////////////////////////////

// Display model parameters
  cout << "Partickle tipes:\n"; 
  for(int i=0; i < particklTipesNum; i++)
  {
    cout<<
      particklTipes[i].name<<" "<<
      particklTipes[i].partTypeID<<" "<<
      particklTipes[i].maxNumber<<" "<<
      particklTipes[i].initNumber<<" "<<
      particklTipes[i].diameter<<" "<<
      particklTipes[i].viscosityPar<<" "<<
      particklTipes[i].randomForceScale<<" "<<
      "\n"<<endl;
  }
  cout << "Partickle interactions:\n";
  l=0;
  for(int i=0; i<particklTipesNum; i++)
    {
      for(int j=0; j<= i; j++)
      {
        cout<<
          i<<" "<<
          j<<" "<<
          particklParams[i][j].funcInteractIndex<<" "<<
          particklParams[i][j].funcParametersNumber<<" "<<
          endl;

        for (int k=0; k<particklParams[i][j].funcParametersNumber; k++)	
	      {
	        cout <<  particklParams[i][j].funcParameters[k]<<" ";
	      }
        cout<<"\n" << endl;
        l++;
      }
    }
  cout<<"Reactions:\n";
  cout << reactionsNumber<<endl;
  for(int i=0; i<reactionsNumber; i++)
  {
    cout<<
      reactions[i].reactionTipe<<" "<<
      reactions[i].reagentsNumber<<" "<<
      reactions[i].productsNumber<<endl;
      
	    for (int k=0; k<reactions[i].reagentsNumber; k++)	
	    {
	       cout <<  reactions[i].reagents[k]<<" ";
	    }
      cout<< endl;
	    for (int k=0; k<reactions[i].productsNumber; k++)	
	    {
	        cout <<   reactions[i].products[k]<<" ";
	    }
      cout << "\n" <<reactions[i].funcReactionIndex<<" " <<reactions[i].funcReactParametersNumber << endl;
      for (int k=0; k<reactions[i].funcReactParametersNumber; k++)	
	      {
	        cout <<  reactions[i].funcReactionParams[k]<<" ";
	      }
      cout<<"\n" << endl;
  }
	cout << "T= "<<T<<" eta= "<< eta<<endl;
  cout<<"Time step, Number of steps, Data collection steps, equilibration steps:\n"<<
    delTime<<" "<<timeStepsNumber<<" "<<dataCollectNumber<<" "<<relaxStepsNumber<<"\n"<<endl;
  cout<<"Size of the box: \n"<<boxSize[0]<<" "<<boxSize[1]<<" "<<boxSize[2]<<"\n"<<endl; 

/////////////////////////////////////////////////////////////////////////////
  totalParticklNumber=0;
  for(int i=0; i < particklTipesNum; i++)
  {
    totalParticklNumber += particklTipes[i].maxNumber; 
  } 

  cout << "totalParticklNumber= "<< totalParticklNumber << endl;
  
  particleArr partickles[totalParticklNumber];
//Initiate partickles
  srand(unsigned(time(0)));
  l=0;
  existed=0;

  for(int i=0; i<particklTipesNum; i++)
  {
    for(int j=0; j<particklTipes[i].maxNumber; j++)
    {
      partickles[l].partTypeID=particklTipes[i].partTypeID;
      if(j<particklTipes[i].initNumber)
      {
        partickles[l].existence=1;
        for(int k=0; k<3; k++)
        {
          partickles[l].coord[k]=boxSize[k]*random()/((double)RAND_MAX);
          partickles[l].force[k]=0.0;
          partickles[l].delBrouncoord[k]=0.0;
        }
        existed++;
      }
      else 
      {
        partickles[l].existence=0;
      }
      l++;
    }
  }
  xyzFile << existed << "\n" <<endl;
  for(int i=0; i<totalParticklNumber; i++)
  {
    if(partickles[i].existence==1)
    {
      xyzFile << particklTipes[partickles[i].partTypeID].name<<" "
      << partickles[i].coord[0]<<" "<< partickles[i].coord[1]<<" "<< partickles[i].coord[2]<<endl;
    } 

  }

   
//Main loop
  Fij = new double[3];

  for(int i=0; i< particklTipesNum; i++)
  {
    for(int j=0; j< particklTipesNum; j++)
    {
      numSqij[i][j]=0;
      avgSqij[i][j]=0.0;
      avgEij[i][j]=0.0;
    }
  }  
  for(int step=1; step<=timeStepsNumber; step++)
  {
    reactionResult = 0;

    for(int i=0; i<totalParticklNumber; i++)
    {
      if(partickles[i].existence==1)
      {
        for(int k=0; k<3; k++)
        {
          partickles[i].force[k]=0.0;
          partickles[i].delBrouncoord[k]=0.0;
          hFi[k]=0.0;
        }
        for(int j=0; j<i; j++) //delRij
        {
          if(partickles[j].existence==1)
          {
            delRij2=0;
            for(int k=0; k<3; k++)
            {
              delRij[k]=partickles[i].coord[k]-partickles[j].coord[k];
              if(periodic)
              {
                delRij[k] -= rint(delRij[k]/boxSize[k])*boxSize[k];
              }

              delRij2 += delRij[k]*delRij[k];
            }
            sqRij=sqrt(delRij2);
            id1=particklParams[partickles[i].partTypeID][partickles[j].partTypeID].funcInteractIndex;
            id2=particklParams[partickles[i].partTypeID][partickles[j].partTypeID].funcParametersNumber;
            for(int k=0; k< id2; k++)
            {
              pars[k]=particklParams[partickles[i].partTypeID][partickles[j].partTypeID].funcParameters[k];
            }
              funInteracArr[id1](pars,delRij,sqRij,Fij,Eij);
//Reactions
              for(int m=0; m<reactionsNumber; m++)
              {
                if(reactions[m].reagentsNumber == 1)
                {
                  if(partickles[i].partTypeID == reactions[m].reagents[0] && partickles[i].existence==1)
                  { 
                    reactionResult=funReactArr[reactions[m].funcReactionIndex](reactions[m].funcReactionParams,sqRij,Eij);
                    if(reactionResult == 1)
                    {
                      partickles[i].existence=0;
                      for(int n=0;n<reactions[m].productsNumber;n++)
                      {                     
                        for(int p=0; p<totalParticklNumber; p++)                       
                        {
                          if(partickles[p].existence ==0 && partickles[p].partTypeID==reactions[m].products[n])
                          {
                            partickles[p].existence = 1;
                            for(int k=0; k<3; k++)
                            {
                              rand1=random()/((double)RAND_MAX);
                              rand2=random()/((double)RAND_MAX);
                              kBT=particklTipes[partickles[p].partTypeID].randomForceScale;
                              randhBr[k]=pow( -2.*(2.*sqRij)*log(rand1+0.000000001) , 0.5 ) * cos(2.*PI*rand2);
                              partickles[p].coord[k] = partickles[i].coord[k] + randhBr[k];
                            }
                            cout<<"step "<<step<<" reaction number "<< m
                                << " product ID "<<reactions[m].products[n] <<endl;
                             n++;
                          }                          
                        }                          
                      }
                    }  
                  }
                }
                if(reactions[m].reagentsNumber == 2)
                {
                  if(partickles[i].existence==1 && partickles[j].existence==1)
                  {
                    if((partickles[i].partTypeID == reactions[m].reagents[0] && partickles[j].partTypeID == reactions[m].reagents[1])
                      || (partickles[i].partTypeID == reactions[m].reagents[1] && partickles[j].partTypeID == reactions[m].reagents[0]))
                    {
                      reactionResult=funReactArr[reactions[m].funcReactionIndex](reactions[m].funcReactionParams,sqRij,Eij);
                      if(reactionResult == 1)
                        {
                          partickles[i].existence=0;
                          partickles[j].existence=0;
                          for(int n=0;n<reactions[m].productsNumber;n++)
                          {    
                            for(int p=0; p<totalParticklNumber; p++)                       
                            {
                              if(partickles[p].existence ==0 && partickles[p].partTypeID==reactions[m].products[n])
                              {
                                partickles[p].existence = 1;
                                for(int k=0; k<3; k++)
                                {
                                  rand1=random()/((double)RAND_MAX);
                                  rand2=random()/((double)RAND_MAX);
                                  randhBr[k]=pow( -2.*(2.*sqRij)*log(rand1+0.000000001) , 0.5 ) * cos(2.*PI*rand2);
                                  partickles[p].coord[k] = partickles[i].coord[k] + randhBr[k];
                                }
                                cout<<"step "<<step<<" reaction number "<< m
                                  << " product ID "<<reactions[m].products[n] <<endl;
                                n++;
                              }                            
                            }
                          }
                        }
                      }
                    }
                  }
                }                        
                if(reactionResult == 0)
                { 
                  xi=3*PI*eta*(particklTipes[partickles[i].partTypeID].diameter+particklTipes[partickles[j].partTypeID].diameter)/2.0;               
                  for(int k=0; k<3; k++)
                  {
                    hFi[k]+=delTime*Fij[k]/xi;
                  }          
                  for(int k=0; k<3; k++)
                  {    
                    rand1=random()/((double)RAND_MAX);
                    rand2=random()/((double)RAND_MAX);
                    kBTS=kBT*particklTipes[partickles[i].partTypeID].randomForceScale;
                    randhBr[k]=pow( -2.*(2.*delTime*kBTS/xi)*log(rand1+0.000000001) , 0.5 ) * cos(2.*PI*rand2);
                    partickles[i].coord[k] = partickles[i].coord[k]+ hFi[k] + randhBr[k];
  
                    rand1=random()/((double)RAND_MAX);
                    rand2=random()/((double)RAND_MAX);
                    kBTS=kBT*particklTipes[partickles[i].partTypeID].randomForceScale;
                    randhBr[k]=pow( -2.*(2.*delTime*kBTS/xi)*log(rand1+0.000000001) , 0.5 ) * cos(2.*PI*rand2);
                    partickles[j].coord[k] = partickles[j].coord[k]- hFi[k] + randhBr[k];
  
                    if(periodic)
                    {
                      partickles[i].coord[k] -= rint(partickles[i].coord[k]/boxSize[k]-0.5)*boxSize[k];
                      partickles[j].coord[k] -= rint(partickles[j].coord[k]/boxSize[k]-0.5)*boxSize[k];
                    }
                  }

                  numSqij[partickles[i].partTypeID][partickles[j].partTypeID]++;
                  nS=numSqij[partickles[i].partTypeID][partickles[j].partTypeID];
                  avgSqij[partickles[i].partTypeID][partickles[j].partTypeID]=
                    (nS-1)*avgSqij[partickles[i].partTypeID][partickles[j].partTypeID]/(double)nS+sqRij/(double)nS;
                  avgEij[partickles[i].partTypeID][partickles[j].partTypeID]=
                    (nS-1)*avgEij[partickles[i].partTypeID][partickles[j].partTypeID]/(double)nS+Eij/(double)nS;

                }           
          }                  
        }        
      }
    }
    if(step > relaxStepsNumber)
    {
      if(step % dataCollectNumber == 0)

      {
        existed=0;
        xyzFileCoords="";
        for(int i=0; i<totalParticklNumber; i++)
        {
          if(partickles[i].existence==1)
          {
            existed++;
            xyzFileCoords = xyzFileCoords + particklTipes[partickles[i].partTypeID].name +" "
            + to_string(partickles[i].coord[0])+" "+to_string(partickles[i].coord[1])+" "+to_string(partickles[i].coord[2])+"\n";
          } 

        }
        
        xyzFile << existed <<endl;
        xyzFile << " " << step <<endl;
        xyzFile << xyzFileCoords <<endl;

        csvFile<< step<<" "<< avgEij[0][0]<< " "<< avgSqij[0][0]<< " "<< avgSqij[1][0]<< " "<< avgSqij[2][0] <<endl;
      }
    }



  }
  delete [] Fij;




//
//  Report timing.
//
  ctime = cpu_time ( ) - ctime;
  cout << "  Elapsed cpu time " << ctime << " seconds.\n";
//

  timestamp ( );
  return 0;
}



double cpu_time ( )
//    Reports the elapsed CPU time.
//    Output, double CPU_TIME, the current total elapsed CPU time in second.
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
//****************************************************************************80

void timestamp ( )
//    Prints the current YMDHMS date as a time stamp.
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << "\n"<< time_buffer << "\n"<<endl;

  return;
# undef TIME_SIZE
}
//****************************************************************************80


void funInterac0(double parameters[], double delRij[3], double sqRij, double Fij[], double &Eij )
{

 
// Lennard-Jones potential energy  

  double eps=parameters[0];
  double sigma=parameters[1];
  double rcoff=parameters[2];

  if(sqRij>rcoff)
  {
    Eij=0.0;
//		cout << sqRij<<endl;	
    for(int k=0; k<3; k++) Fij[k]=0.0;
    return;   
  }


  double r=0.000001;

  double tmp=sigma/sqRij;
  double tmp2=tmp*tmp;
  double tmp6=tmp2*tmp2*tmp2;
  double tmp12=tmp6*tmp6;
  Eij=4*eps*(tmp12-tmp6);
//	cout << Eij<<endl;
  for(int i=0; i<3; i++)
  {
    Fij[i]=Eij*delRij[i]/sqRij;
  }
 
  return;
}

void funInterac1 (double parameters[], double delRij[3], double sqRij, double Fij[], double &Eij)
{ 
  double eps=parameters[0];
  double sigma=parameters[1];
  double rcoff=parameters[2];

  if(sqRij>rcoff)
  {
    Eij=0.0;
    for(int k=0; k<3; k++) Fij[k]=0.0;
    return;   
  }

  double r=0.000001;

  double tmp=sigma/sqRij;
  double tmp2=tmp*tmp;
  double tmp6=tmp2*tmp2*tmp2;
  double tmp12=tmp6*tmp6;
  Eij=4*eps*(tmp12-tmp6);

  for(int i=0; i<3; i++)
  {
    Fij[i]=Eij*delRij[i]/sqRij;
  }

  return;
}

int funReact0(double parameters[], double sqRij, double Eij)
{ 
  double Energy=parameters[0];
  double probability=parameters[1];
  double randprob;
  randprob = random()/((double)RAND_MAX);

  if(Eij>Energy && randprob > probability)
  {
    return(1);
  }
  else return(0);  
}

int funReact1(double parameters[], double sqRij, double Eij)
{ 
  double probability=parameters[0];
  double Energy=parameters[1];

  if(Eij>Energy && random()/((double)RAND_MAX)>probability)
  {
    return(1);
  }
  else return(0);  
}




