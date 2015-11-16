#include <math.h>
#include <stdlib.h>
#include <vector>

struct Bisection {

static void get_eig2(int row_num,std::vector<RealD> &ALPHA,std::vector<RealD> &BETA, std::vector<RealD> & eig)
{
  int i,j;
  std::vector<RealD> evec1(row_num+3);
  std::vector<RealD> evec2(row_num+3);
  RealD eps2;
  ALPHA[1]=0.;
  BETHA[1]=0.;
  for(i=0;i<row_num-1;i++) {
    ALPHA[i+1] = A[i*(row_num+1)].real();
    BETHA[i+2] = A[i*(row_num+1)+1].real();
  }
  ALPHA[row_num] = A[(row_num-1)*(row_num+1)].real();
  bisec(ALPHA,BETHA,row_num,1,row_num,1e-10,1e-10,evec1,eps2);
  bisec(ALPHA,BETHA,row_num,1,row_num,1e-16,1e-16,evec2,eps2);

  // Do we really need to sort here?
  int begin=1;
  int end = row_num;
  int swapped=1;
  while(swapped) {
    swapped=0;
    for(i=begin;i<end;i++){
      if(mag(evec2[i])>mag(evec2[i+1]))	{
	swap(evec2+i,evec2+i+1);
	swapped=1;
      }
    }
    end--;
    for(i=end-1;i>=begin;i--){
      if(mag(evec2[i])>mag(evec2[i+1]))	{
	swap(evec2+i,evec2+i+1);
	swapped=1;
      }
    }
    begin++;
  }

  for(i=0;i<row_num;i++){
    for(j=0;j<row_num;j++) {
      if(i==j) H[i*row_num+j]=evec2[i+1];
      else H[i*row_num+j]=0.;
    }
  }
}

static void bisec(std::vector<RealD> &c,   
		  std::vector<RealD> &b,
		  int n,
		  int m1,
		  int m2,
		  RealD eps1,
		  RealD relfeh,
		  std::vector<RealD> &x,
		  RealD &eps2)
{
  std::vector<RealD> wu(n+2);

  RealD h,q,x1,xu,x0,xmin,xmax; 
  int i,a,k;

  b[1]=0.0;
  xmin=c[n]-fabs(b[n]);
  xmax=c[n]+fabs(b[n]);
  for(i=1;i<n;i++){
    h=fabs(b[i])+fabs(b[i+1]);
    if(c[i]+h>xmax) xmax= c[i]+h;
    if(c[i]-h<xmin) xmin= c[i]-h;
  }
  xmax *=2.;

  eps2=relfeh*((xmin+xmax)>0.0 ? xmax : -xmin);
  if(eps1<=0.0) eps1=eps2;
  eps2=0.5*eps1+7.0*(eps2);
  x0=xmax;
  for(i=m1;i<=m2;i++){
    x[i]=xmax;
    wu[i]=xmin;
  }

  for(k=m2;k>=m1;k--){
    xu=xmin;
    i=k;
    do{
      if(xu<wu[i]){
	xu=wu[i];
	i=m1-1;
      }
      i--;
    }while(i>=m1);
    if(x0>x[k]) x0=x[k];
    while((x0-xu)>2*relfeh*(fabs(xu)+fabs(x0))+eps1){
      x1=(xu+x0)/2;

      a=0;
      q=1.0;
      for(i=1;i<=n;i++){
	q=c[i]-x1-((q!=0.0)? b[i]*b[i]/q:fabs(b[i])/relfeh);
	if(q<0) a++;
      }
      //			printf("x1=%e a=%d\n",x1,a);
      if(a<k){
	if(a<m1){
	  xu=x1;
	  wu[m1]=x1;
	}else {
	  xu=x1;
	  wu[a+1]=x1;
	  if(x[a]>x1) x[a]=x1;
	}
      }else x0=x1;
    }
    x[k]=(x0+xu)/2;
  }
}
}
