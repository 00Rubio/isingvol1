#include <stdio.h>
#include <math.h>
#include "gsl_rng.h"

//N la cambiaré para cada prueba
#define N 16

gsl_rng *tau;

int main ()
{
//Defino cada variable que utilice
int orden,i,j,k,l,n,m,pmc,t,r;
int s[N][N];
double E,p,ex,min,eps,T;
double sumas,Mn,Mnp,Et,Etp,E2,E2p,En,Enp,Cnp,f[N],frp,fpsum,fc[N];
FILE *f1, *f2, *f3, *f4;

//Inicializo algunas variables
Mn=0;
En=0;
t=0;
pmc=0;
sumas=0;

//gsl
extern gsl_rng *tau;
int semilla= 2346548967;
tau=gsl_rng_alloc(gsl_rng_taus);
gsl_rng_set(tau,semilla);

//Abrimos fichero de escritura
f1=fopen("Magnetizacion16.txt","a");
f2=fopen("Energia_media16.txt","a");
f3=fopen("Calor_específico16.txt","a");
f4=fopen("Funcion_correlacion16.txt","a");

//Pedimos si queremos un estado inicial ordenado o desordenado
orden=1;


//Escribo en cada fichero a qué N corresponde

    fprintf(f1,"N=%i\n",N);
    fprintf(f2,"N=%i\n",N);
    fprintf(f3,"N=%i\n",N);
    fprintf(f4,"N=%i\n",N);


//Comienzo el bucle variando la temperatura para generar todo de una vez
    for(T=1.5;T<=3.7;T+=0.2)
    {

//Inicializo algunas variables
Mn=0;
En=0;
t=0;
pmc=0;
sumas=0;
E2=0;
Et=0;
orden=1;

//GENERAMOS LAS CONDICIONES INICIALES
if(orden==1)
	{
		for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				s[i][j]=1;
			}
		}
	}

else{
	for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				k=gsl_rng_uniform_int(tau,2);
				if(k==0)
					{s[i][j]=-1;}
				else 
					{s[i][j]=1;}
			}
		}
	}

k=0;




    //COMENZAMOS EL BUCLE

    
        for(i=0;i<pow(N,2)*1000000;i++)
        {
	        //Elegimos un punto al azar
	        n=gsl_rng_uniform_int(tau,N);
	        m=gsl_rng_uniform_int(tau,N);

	        //Calculamos la variación de energía en ese punto con las condiciones de contorno
	        if((n==0) && (m==0))
	        {
		        E=2*s[n][m]*(s[n+1][m]+s[N-1][m]+s[n][m+1]+s[n][N-1]);
	        }
	        else if((n==0) && (m==N-1))
	        {
		        E=2*s[n][m]*(s[n+1][m]+s[N-1][m]+s[n][0]+s[n][m-1]);
	        }
	        else if((n==N-1) && (m==0))
	        {
		        E=2*s[n][m]*(s[0][m]+s[n-1][m]+s[n][m+1]+s[n][N-1]);
	        }
	        else if((n==N-1) && (m==N-1))
	        {
		        E=2*s[n][m]*(s[0][m]+s[n-1][m]+s[n][0]+s[n][m-1]);
	        }
	        else if((m==0) && (n!=N-1) && (n!=0))
	        {
		        E=2*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][m+1]+s[n][N-1]);
	        }
	        else if((m==N-1) && (n!=N-1) && (n!=0))
	        {
		        E=2*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][0]+s[n][m-1]);
	        }
	        else if((n==0) && (m!=N-1) && (m!=0))
	        {
		        E=2*s[n][m]*(s[n+1][m]+s[N-1][m]+s[n][m+1]+s[n][m-1]);
	        }
	        else if((n==N-1) && (m!=N-1) && (m!=0))
	        {
		        E=2*s[n][m]*(s[0][m]+s[n-1][m]+s[n][m+1]+s[n][m-1]);
	        }
	        else
	        {
		        E=2*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][m+1]+s[n][m-1]);
	        }



	        //Evaluamos
	        ex=exp(-E/T);
	        if(ex<1){min=ex;}
	        else{min=1;}

	        //Generamos epsilon entre 0 y 1
	        eps=gsl_rng_uniform(tau);

	        //Vemos si cambiamos o no
	        if(eps<min){s[n][m]=-s[n][m];}



        //CALCULO LAS MAGNITUDES CADA 100 PMC
	        l++;
	        if(l==pow(N,2)){pmc++; l=0;}
	        if(pmc%100==0)
	        {

		        //Contabilizo cada 100pmc
		        t++;

		        //Cálculos para la magnetización promedio
		        for(j=0;j<N;j++)
			        {
				        for(k=0;k<N;k++)
				        {
					        sumas=sumas+s[j][k];
				        }
			        }					
			        Mn=Mn+((1/pow(N,2))*fabs(sumas));
                    sumas=0;

		        //Cálculos para la Energía media
		        Et=Et+E;
		        E2=E2+pow(E,2);
			        
		        //Cálculos para la función correlación

                frp=0;
                fpsum=0;

		        for(r=0;r<N/2;r++)
		        {
		            for(j=0;j<N/2;j++)
		            {
			            for(k=0;k<N;k++)
			            {

					            f[r]=s[j][k]*s[j+r][k];
					            frp=frp+(f[r]/t);
				            
			            }
                    fpsum=fpsum+frp;
		            }

                }
 
	        }


        }

	    //CALCULO LOS PROMEDIOS DE LAS MAGNITUDES	
	    Mnp=Mn/t;

	    Enp=Et/(t*2*N);

	    E2p=E2/t;
	    Cnp=(1/(N*N*T))*(E2p-pow(Etp,2));


	    for(r=0;r<N/2;r++)
	    {
		    fc[r]=(1/pow(N,2))*fpsum;	
	    }	
	    

	    //Saco a un fichero los resultados obtenidos
	    fprintf(f1,"%lf\t%lf\n",T,Mnp);
	    fprintf(f2,"%lf\t%lf\n",T,Enp);
	    fprintf(f3,"%lf\t%lf\n",T,Cnp);
		fprintf(f4,"%lf\t%lf\n",T,fc[1]);

    	
    }



//Cerramos fichero de escritura
fclose(f1);
fclose(f2);
fclose(f3);
fclose(f4);


return 0;
}
