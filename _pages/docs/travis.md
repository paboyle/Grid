---
layout: single
title : "Documentation"
author_profile: false
excerpt: "Travis Continous Integration framework"
header:
  overlay_color: "#5DADE2"
permalink: /docs/travis/
sidebar:
  nav : docs
---
{% include base_path %}
<!-- {% include toc icon="gears" title="Helpers" %} -->

For continous testing of every commit we use the [Travis CI framework](https://travis-ci.org/). 

The current status is 

| Branch       |  Status |
|----------    |  ------ |
| [Master](https://travis-ci.org/paboyle/Grid)       |  [![Build Status](https://travis-ci.org/paboyle/Grid.svg?branch=master)](https://travis-ci.org/paboyle/Grid) |
| [Develop](https://travis-ci.org/paboyle/Grid)       | [![Build Status](https://travis-ci.org/paboyle/Grid.svg?branch=develop)](https://travis-ci.org/paboyle/Grid) |
| [Release 0.7.0](https://github.com/paboyle/Grid/tree/release/v0.7.0)       | [![Build Status](https://travis-ci.org/paboyle/Grid.svg?branch=release/v0.7.0)](https://github.com/paboyle/Grid/tree/release/v0.7.0) |


### Automated tests

Travis will test the compilation workflow for single and double precision version on the following compilers for each commit:

 - clang 3.7.0 on Ubuntu 14.04
 - clang 3.8.0 on Ubuntu 14.04
 - gcc 5.4.1 on Ubuntu 14.04
 - gcc 4.9.4 on Ubuntu 14.04.1
 - Apple LLVM version 8.1.0 (clang-802.0.42) on OSX (x86_64-apple-darwin16.5.0)

Due to the limitations of the Travis virtual machines, the archictecture is limited to SSE4 and few tests. 

May 2017: a new server using [TeamCity](https://www.jetbrains.com/teamcity/specials/teamcity/teamcity.html?gclid=CjwKEAjwutXIBRDV7-SDvdiNsUoSJACIlTqlygt_V8-PqWvjV23oAj8wf2suNmct9-sFfplBFYctzBoCnTvw_wcB&gclsrc=aw.ds.ds&dclid=COOh9rPt6dMCFYOmUQodkpwLfQ) is being setup for extensive testing of every commit. 


{% include paginator.html %}