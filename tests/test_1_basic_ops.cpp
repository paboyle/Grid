#include "catch.hpp"

#include "GridSIMD.hpp"


TEST_CASE ("Reading input with Python")
{

  cout << "Hello!" << endl;
  cout << "n SIMD lanes is " << n_simd_lanes << endl;


  const long N = 1000;

  vReal1 vec1;
  vReal1 vec2;
  vReal1 vec3;
  vReal1 vec4;

  vec1.resize (N);
  vec2.resize (N);
  vec3.resize (N);
  vec4.resize (N);



  long index = 0;


  for (long n = 0; n < N; n++)
  {
    index++;

    for (int lane = 0; lane < n_simd_lanes; lane++)
    {
      vec1[n].putlane(index, lane);
      vec2[n].putlane(index, lane);
    }
  }


  for (long n = 0; n < N; n++)
  {
    vec3[n] = vec1[n] + vec2[n];
    vec4[n] = vec3[n] / vec1[n];
  }


  for (long n = 0; n < N; n++)
  {
    for (int lane = 0; lane < n_simd_lanes; lane++)
    {
      double v3 = vec3[n].getlane(lane);
      double v4 = vec4[n].getlane(lane);

    CHECK (v4 == 2.0);
    }
  }

}
