---
title : "API Documentation"
author_profile: false
excerpt: "Coordinates"
header:
  overlay_color: "#5DADE2"
permalink: /docs/API/coordinates.html
sidebar:
  nav : docs
---

The Grid is define on a N-dimensional set of integer coordinates. 

The maximum dimension is eight, and indexes in this space make use of the `Coordinate` class.
The coordinate class shares a similar interface to `std::vector<int>`, but contains all data within the
object, and has a fixed maximum length (template parameter).

**Example**:

```c++
const int Nd=4;
Coordinate point(Nd);

for(int i=0;i<Nd;i++) 
   point[i] = 1;

std::cout<< point <<std::endl;
point.resize(3);
std::cout<< point <<std::endl;
```

This enables the coordinates to be manipulated without heap allocation or thread contention,
and avoids introducing STL functions into GPU code, but does so at the expense of introducing
a maximum dimensionality. This limit is easy to change (`lib/util/Coordinate.h`).
