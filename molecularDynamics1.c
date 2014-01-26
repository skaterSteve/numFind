#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double f(double y)
{
  return -y - y*y*y;
}

int main(void)
{
  int N = 20;
  double v[N];
  int k, num;
  double sum;
  int taken[N];
  double E = (double)N/2.0;   //energy is constant

  srand(time(NULL));

  for( k = 0; k < N; k++ )
  {
    v[k] = -1.0;
    taken[k] = 0;
  }

  for( k = 0; k < N/2; k++)
  {
    while(1)
    {
	num = rand()%N;
	if (taken[num] == 0)
	{
	    taken[num] = 1;
	    break;
	}
    }
    v[num] = 1.0;
  }


  sum = 0;
  for(k = 0; k < N; k++)
    {
      printf("v[%2d] = %10.4f\n",k,v[k]);
      sum += v[k];
    }
  printf("sum = %10.4f\n", sum);

  printf("\n\n\n");
  /*------------------------------------------*/

  int time, T = 500000;
  double y[N];     //define equilibrium position
  double h = 0.01; //define step size

  
  for(k = 0; k < N; k++)
  {
     y[k] = 0;           //initially at equil = 0 point
  }
  

  for (time = 0; time < T; time++)    //position verlet
  {
    // printf("\ntime = %d:\n", time+1);
    //printf("k    v[k]\n");
      for (k = 0; k < N; k++)
      {       
	y[k] += 0.5 * h * v[k];
      }

     
      v[0] += h * ( f( y[0] - y[1] ) - f( y[N-1] - y[0] ) );
      for (k = 1; k < N-1; k++)
      {
	v[k] += h * ( f( y[k] - y[k+1] ) - f( y[k-1] - y[k] ) );
      }
      v[N-1] += h * ( f( y[N-1] - y[0] ) - f( y[N-2] - y[N-1] ) );


      for (k = 0; k < N; k++)
      {
	y[k] += 0.5 * h * v[k];

	//	printf("%d    %f\n", k, v[k]);   //print evry dt
      }	
 

  }  

  
  printf("After %d timesteps:\n", time);
  printf("k   v[k]\n-----------\n");
  for( k = 0; k < N; k++)
    {
      printf("%d    %f\n", k+1, v[k]);
    }
  
  double max = v[0];
  double min = v[0];


  for(k = 0; k < N; k++)
    {
      if(v[k] > max)
	{
	  max = v[k];
       
	}
    }

  for(k = 0; k < N; k++)
    {
      if (v[k] < min)
	{
	  min = v[k];
	}
    }
  printf("\n\n");
  printf("MAX vel is: %f\n", max);
  printf("MIN vel is: %f\n", min);

  /*    begin binning algorithm    */
  int m;
  int count[N/2+1];
  double bin[N/2+1];
 
  printf("bin#   binvalue   counts\n");
  double store[N];

  for( k = 0; k< N; k++)
    {
      store[k] = v[k];
      //  printf("%f\n", store[k]);

      for(m=0; m< N/2+1; m++)
	{
	   bin[m] = min + 2 * m * (max - min)/(double)N;
	   count[m] = 0;
	   if(bin[m]>=store[k] && bin[m+1] < store[k])
	     {
	       count[m] = count[m] +1;
	     }
	}
      if( bin[N/2] <= store[k] && bin[N/2-1] > store[k] )
	{
	  count[N/2]=count[N/2]+1;
	}
    }
  printf("bin#   counts\n");
  for( m = 0; m < N/2+1; m++)
    {
      printf("%d   %d\n", m, count[m]);
    }

    for(m=0; m< N/2+1; m++)
      {
	printf("%d   %f\n", m, bin[m]);   
      }


  /*
  for(k=0; k < N-1; k++)
    {
        for(i=0; i < N/2; i++)
	  {
	    count[i] = 0;
	    bin[i] = min + 2 * i * (max - min)/(double)N;
	    //	    printf("%d   %f\n", i, bin[i]);    
	  }
      
	if( v[k] < bin[i+1] && v[k] >= bin[i])
	  {
	    count[i]= count[i]+1;
	  }
	else;
    }
  
  for(i=0; i < N/2; i++)
    {
      printf("%d   %d\n", i,  count[i]);
    }
  
  */
  return 0;
}
