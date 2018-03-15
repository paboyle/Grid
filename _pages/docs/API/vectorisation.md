---
title : "API Documentation"
author_profile: false
excerpt: "Vectorisation"
header:
  overlay_color: "#5DADE2"
permalink: /docs/API/vectorisation.html
sidebar:
  nav : docs
---


Internally, Grid defines a portable abstraction SIMD vectorisation, via the following types:

* `vRealF`

* `vRealD`

* `vComplexF`

* `vComplexD`

These have the usual range of arithmetic operators and functions acting upon them. They do not form
part of the API, but are mentioned to (partially) explain the need for controlling the
layout transformation in lattice objects. 