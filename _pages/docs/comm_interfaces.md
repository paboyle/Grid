---
title : "Documentation"
author_profile: false
excerpt: "Supported communication interfaces"
header:
  overlay_color: "#5DADE2"
permalink: /docs/comm_interfaces/
sidebar:
  nav : docs
---
{% include base_path %}

The following options can be used with the `--enable-comms=` option to target different communication interfaces:

| `<comm>`      | Description                                  |
| ------------- | -------------------------------------------- |
| `none`        | no communications                            |
| `mpi`         | MPI communications using [MPI-3 shared memory](https://software.intel.com/sites/default/files/managed/eb/54/An_Introduction_to_MPI-3.pdf) |
| `mpi-auto`    | MPI communications with compiler CXX but clone flags from MPICXX |

For the MPI interfaces the optional `-auto` suffix instructs the `configure` scripts to determine all the necessary compilation and linking flags. This is done by extracting the informations from the MPI wrapper specified in the environment variable `MPICXX` (if not specified `configure` will scan though a list of default names). The `-auto` suffix is not supported by the Cray environment wrapper scripts. Use the standard wrappers ( `CXX=CC` ) set up by Cray `PrgEnv` modules instead. 

### Shared memory communications support

The following options can be use with the `--enable-shm=` option to target different shared memory behaviours (default `shmopen`):

| `<shm>`       | Description                                                       |
| ------------- | ----------------------------------------------------------------- |
| `shmnone`     | uses anonymous spaces, use only for 1 MPI rank per node           |
| `shmopen`     | uses `shm_open` to allocate a shared memory space for inter socket communications. It used a unique file name in `/dev/shm` associated to the user. |
| `hugetlbfs`   | optional [libhugetlbfs](https://github.com/libhugetlbfs/libhugetlbfs) support to map the shared memory allocation into huge 2M pages                |


### Other flags

`--enable-shmpath=<path>`  to select the sharem memory map base path for [libhugetlbfs](https://github.com/libhugetlbfs/libhugetlbfs).

{% include paginator.html %}