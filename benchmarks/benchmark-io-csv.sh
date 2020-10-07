#!/usr/bin/env bash

awkscript='
BEGIN{
  i = 0;
  print "local L,std read (MB/s),std write (MB/s),Grid Lime read (MB/s),Grid Lime write (MB/s)"
}

/Benchmark std write/{
  i    = 0; 
  mode = "stdWrite";
} 

/Benchmark std read/{
  i    = 0; 
  mode = "stdRead"
} 

/Benchmark Grid C-Lime write/{
  i    = 0; 
  mode = "gridWrite";
} 

/Benchmark Grid C-Lime read/{
  i    = 0; 
  mode = "gridRead";
}

/Local volume/{
  match($0, "[0-9]+\\^4");
  l[i] = substr($0, RSTART, RLENGTH-2);
}

/MB\/s/{
  match($0, "[0-9.eE]+ MB/s");
  p = substr($0, RSTART, RLENGTH-5);
  if (mode == "stdWrite")
  {
    sw[i] = p;
  }
  else if (mode == "stdRead")
  {
    sr[i] = p;
  }
  else if (mode == "gridWrite")
  {
    gw[i] = p;
  }
  else if (mode == "gridRead")
  {
    gr[i] = p;
  }
  i++;
}

END{
  s = 0
  for (a in l)
  {
    s++;
  }
  for (j = 0; j < s; j++)
  {
    printf("%s,%s,%s,%s,%s\n", l[j], sr[j], sw[j], gr[j], gw[j]);
  }
  printf("\n");
}
'

if (( $# != 1 )); then
    echo "usage: `basename $0` <log file>" 1>&2
    exit 1
fi
LOG=$1

awk "${awkscript}" ${LOG} 
