SUBDIRS = . 

include Make.inc

# Precision variants of the two multigrid drivers.  The coarse-sector and
# dense-inversion precisions are compile-time instantiations, so each is its
# own binary built from the same source; the fine-level precision is a
# run-time parameter and needs no variant.  Per-target flags give each
# variant its own object file.
bin_PROGRAMS += Example_pvdagm_multigrid_fp32coarse \
                Example_pvdagm_multigrid_fp32coarse_fp32dense \
                Example_hdcg_multigrid_fp32coarse

Example_pvdagm_multigrid_fp32coarse_SOURCES  = Example_pvdagm_multigrid.cc
Example_pvdagm_multigrid_fp32coarse_CPPFLAGS = -DCOARSE_SINGLE
Example_pvdagm_multigrid_fp32coarse_LDADD    = $(top_builddir)/Grid/libGrid.a

Example_pvdagm_multigrid_fp32coarse_fp32dense_SOURCES  = Example_pvdagm_multigrid.cc
Example_pvdagm_multigrid_fp32coarse_fp32dense_CPPFLAGS = -DCOARSE_SINGLE -DGRID_DENSE_INVERSE_SINGLE
Example_pvdagm_multigrid_fp32coarse_fp32dense_LDADD    = $(top_builddir)/Grid/libGrid.a

Example_hdcg_multigrid_fp32coarse_SOURCES  = Example_hdcg_multigrid.cc
Example_hdcg_multigrid_fp32coarse_CPPFLAGS = -DCOARSE_SINGLE
Example_hdcg_multigrid_fp32coarse_LDADD    = $(top_builddir)/Grid/libGrid.a
