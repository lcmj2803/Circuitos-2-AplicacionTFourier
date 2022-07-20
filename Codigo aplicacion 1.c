#include <stdio.h>
#include <math.h>
#define TWOPI 6.28318530717959
int main (void)
{
double fs = 10000; /* fs que mencionamos  [Hz] */
double f1 = 1234; /* la señal numero 1 [Hz] */
double amp1 = 2.82842712474619; /* 2 Vrms */
double f2 = 2500.2157; /* la señal 2 [Hz] */
double amp2 = 1; /* 0.707 Vrms */
double ulsb = 1e-3; /* Value of 1 LSB in Volt */
int i;
double t, u, ur;
for (i = 0; i < 1000000; i++)
{
t = (double) i / fs;
u = amp1 * sin (TWOPI * f1 * t) + amp2 * sin (TWOPI * f2 * t);
ur = floor (u / ulsb + 0.5) * ulsb; 
printf ("%10.6f %8.5f\n", t, ur); 
fwrite (&ur, sizeof (double), 1, stdout); 
}
return 0;
}

/* un par de procedimientos no estan explicados porque no tenian que ver con la transformada 
en el propio documento se menciona cuando se utiliza */