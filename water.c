#include "stdheader.h"
#include "water.h"

Water createwater(double x, double y, double z)
{
int Nx=1,Ny=1,Nz=1;
int xcm=0;ycm=0;zcm=-0.0616571936;
Water water;
water.xo=0;
water.yo=0;
water.zo=0;
water.xh1=  0.952;
water.yh1=0;
water.zh1=-0.551;
water.xh2=-0.952;
water.yh2=0;
water.zh2=-0.551;
water.eps[0] = 0.2104;
water.eps[1]=0;
water.eps[2]=0;
water.sigma[1] = 3.06647338784;
water.sigma[2]=0,
water.sigma[3]=0;
water.zo=water.zo-zcm;
water.zh1=water.zh1+zcm;
water.zh2=water.zh2+zcm;
printf("%lf %lf %lf",water.zo, 
}


