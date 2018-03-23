---
title : "Documentation"
author_profile: false
excerpt: "Running on a Intel Knights Landing"
header:
  overlay_color: "#5DADE2"
permalink: /docs/running_knl/
sidebar:
  nav : docs
---
{% include base_path %}
{% include toc icon="gears" title="Contents" %}

These are few suggestions in order to get the best performance on the Intel Knights Landing (KNL). 



### Bind the memory allocation to the MCDRAM NUMA node

The KNL  has two memory systems, the DDR4 (~90 GFlops/s) and the High Bandwidth Memory (MCDRAM, ~400 Gflops/s).
Each of the two memory system is attached to a different [NUMA context](https://software.intel.com/en-us/articles/optimizing-applications-for-numa).

On a KNL node the command `numactl --hardware` will report which NUMA context is connected to the faster MCDRAM. 
A typical report looks like this 

``` text
node 0 size: 98178 MB
node 0 free: 92899 MB
node 1 cpus:
node 1 size: 16384 MB
node 1 free: 15926 MB
```

In this case the node 1 is related to the 16GB MCDRAM (this is the typical situation on KNLs)

To bind the memory allocation to NUMA node 1 use

``` bash
numactl --membind 1 ./your-executable
```

### Controlling threading

The number of threads can be set in GRID at runtime by the flag

``` text
--threads <#threads>
```

A finer control can be achieved using the environment variable `KMP_HW_SUBSETS` (or the deprecated `KMP_PLACE_THREADS`).

From the [Intel developer guide](https://software.intel.com/en-us/node/684313):

>The KMP_HW_SUBSETS variable controls the hardware resource that will be used by the program. This variable specifies the number of sockets to use, how many cores to use per socket and how many threads to assign per core. For example, on Intel® Xeon Phi™ coprocessors, while each coprocessor can take up to four threads, specifying fewer than four threads per core may result in a better performance. While specifying two threads per core often yields better performance than one thread per core, specifying three or four threads per core may or may not improve the performance. This variable enables you to conveniently measure the performance of up to four threads per core.

A typical setting for the best performance on a single node is to use **62 cores with 1 threads per code**, on the bash shell this is set by

``` bash
export KMP_HW_SUBSETS=62c,1t
```

### Using the optimised Wilson Dslash kernels

Beside the generic implementation using stencils, GRID has optimised version of the Dslash kernels (for Wilson and DWF fermions). 

Flags at runtime can be used for the optimised paths

| Flag                  | Description                                                           |
| --------------------- | --------------------------------------------------------------------- |
| `--dslash-generic`  	| This is the default option and used the implementation with stencils 	|
| `--dslash-unroll`   	| This explicitly unroll the colour loops. It is tied to `Nc=3`        	|
| `--dslash-asm`      	| This is specific for AVX512-F architectures and `Nc=3`               	|



The information included in this page has been updated on *November 2016* and it is valid for the release version 0.6.0.
{: .notice}


{% include paginator.html %}