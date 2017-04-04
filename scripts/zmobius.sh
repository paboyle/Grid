#!/bin/bash
fn=$1

grep "double zmobius_" $fn |
awk 'BEGIN{ m["zmobius_b_coeff"]=0; m["zmobius_c_coeff"]=1; }{ val[m[substr($2,0,15)]][substr($2,17)+0]=$4; }END{

    ls=length(val[0])/2;

    print "ls = " ls

    bmc=-111;

    for (s=0;s<ls;s++) {
      br[s] = val[0][2*s + 0];
      bi[s] = val[0][2*s + 1];
      cr[s] = val[1][2*s + 0];
      ci[s] = val[1][2*s + 1];

      t=br[s] - cr[s];
      if (bmc == -111)
        bmc=t;
      else if (bmc != t)
        print "Warning: b-c is not constant!";

      omegar[s] = (-1.0 + 2.0* br[s])/(4.0*bi[s]**2.0 + (1.0 - 2.0* br[s])**2);
      omegai[s] = - 2.0* bi[s]/(4.0*bi[s]**2.0 + (1.0 - 2.0* br[s])**2);
    }

    print "b-c = " bmc

    for (s=0;s<ls;s++) {
      printf( "omega.push_back( std::complex<double>(%.15g,%.15g) );\n",omegar[s],omegai[s]);
    }

}'
