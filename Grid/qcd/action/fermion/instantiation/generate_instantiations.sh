#!/bin/sh

STAG_IMPL_LIST=" \
	   StaggeredImplF \
	   StaggeredImplD  "

STAG5_IMPL_LIST=""

WILSON_IMPL_LIST=" \
	   WilsonImplF \
	   WilsonImplD \
	   WilsonAdjImplF \
	   WilsonAdjImplD \
	   WilsonTwoIndexSymmetricImplF \
	   WilsonTwoIndexSymmetricImplD \
	   WilsonTwoIndexAntiSymmetricImplF \
	   WilsonTwoIndexAntiSymmetricImplD \
	   GparityWilsonImplF \
	   GparityWilsonImplD "

COMPACT_WILSON_IMPL_LIST=" \
	   WilsonImplF \
	   WilsonImplD "

DWF_IMPL_LIST=" \
	   WilsonImplF \
	   WilsonImplD \
	   ZWilsonImplF \
	   ZWilsonImplD "

GDWF_IMPL_LIST=" \
	   GparityWilsonImplF \
	   GparityWilsonImplD "

IMPL_LIST="$STAG_IMPL_LIST  $WILSON_IMPL_LIST $DWF_IMPL_LIST $GDWF_IMPL_LIST"

for impl in $IMPL_LIST
do
  echo $impl
  mkdir -p $impl
cat > $impl/impl.h <<EOF
#define IMPLEMENTATION $impl
EOF

done

CC_LIST="WilsonCloverFermionInstantiation WilsonFermionInstantiation WilsonKernelsInstantiation WilsonTMFermionInstantiation"

for impl in $WILSON_IMPL_LIST
do
for f in $CC_LIST
do
  ln -f -s ../$f.cc.master $impl/$f$impl.cc
done
done

CC_LIST="CompactWilsonCloverFermionInstantiation"

for impl in $COMPACT_WILSON_IMPL_LIST
do
for f in $CC_LIST
do
  ln -f -s ../$f.cc.master $impl/$f$impl.cc
done
done

CC_LIST=" \
  CayleyFermion5DInstantiation \
  ContinuedFractionFermion5DInstantiation \
  DomainWallEOFAFermionInstantiation  \
  MobiusEOFAFermionInstantiation \
  PartialFractionFermion5DInstantiation \
  WilsonFermion5DInstantiation \
  WilsonKernelsInstantiation "

for impl in $DWF_IMPL_LIST $GDWF_IMPL_LIST
do
for f in $CC_LIST
do
  ln -f -s ../$f.cc.master $impl/$f$impl.cc
done
done

# overwrite the .cc file in Gparity directories
for impl in $GDWF_IMPL_LIST
do
  ln -f -s ../WilsonKernelsInstantiationGparity.cc.master $impl/WilsonKernelsInstantiation$impl.cc
done


CC_LIST=" \
  ImprovedStaggeredFermion5DInstantiation \
  ImprovedStaggeredFermionInstantiation \
  NaiveStaggeredFermionInstantiation \
  StaggeredKernelsInstantiation "

for impl in $STAG_IMPL_LIST
do
for f in $CC_LIST
do
  ln -f -s ../$f.cc.master $impl/$f$impl.cc
done
done

CC_LIST=" \
  ImprovedStaggeredFermion5DInstantiation \
  StaggeredKernelsInstantiation "

