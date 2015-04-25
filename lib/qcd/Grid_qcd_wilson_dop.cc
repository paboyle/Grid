#ifnfdef GRID_QCD_WILSON_DOP_H
#define  GRID_QCD_WILSON_DOP_H

#include <Grid.h>

namespace Grid {
namespace QCD {


const std::vector<int> WilsonMatrix::directions   ({0,1,2,3, 0, 1, 2, 3,0});
const std::vector<int> WilsonMatrix::displacements({1,1,1,1,-1,-1,-1,-1,0});

  // Should be in header?
static const int WilsonMatrix::Xp = 0;
static const int WilsonMatrix::Yp = 1;
static const int WilsonMatrix::Zp = 2;
static const int WilsonMatrix::Tp = 3;
static const int WilsonMatrix::Xm = 4;
static const int WilsonMatrix::Ym = 5;
static const int WilsonMatrix::Zm = 6;
static const int WilsonMatrix::Tm = 7;
static const int WilsonMatrix::X0 = 8;
static const int WilsonMatrix::npoint=9;


WilsonMatrix::WilsonMatrix(LatticeGaugeField &_Umu,int _mass)
  : Stencil((&Umu._grid,npoint,0,directions,displacements),
    mass(_mass),
    Umu(_Umu)
{
  // Allocate the required comms buffer
  grid = _Umu._grid;
  comm_buf.resize(Stencil._unified_buffer_size);
}
void WilsonMatrix::multiply(const LatticeFermion &in, LatticeFermion &out)
{
  
}
void WilsonMatrix::Dhop(const LatticeFermion &in, LatticeFermion &out)
{
  Stencil.HaloExchange(in,comm_buf);

  for(int ss=0;ss<_grid->oSites();ss++){

    int offset,local;

    vSpinColourVector result;
    vHalfSpinColourVector UChi;

    // Xp
    offset = Stencil._offsets [Xp][ss];
    local  = Stencil._is_local[Xp][ss];
    if ( local ) {
      Uchi = U[]*spProjXp(in._odata[offset]);
    } else {
      Uchi = U[]*comm_buf._odata[offset]
    }
    result = ReconXp(Uchi);

    // Yp
    offset = Stencil._offsets [Yp][ss];
    local  = Stencil._is_local[Yp][ss];
    if ( local ) {
      Uchi = U[]*spProjYp(in._odata[offset]);
    } else {
      Uchi = U[]*comm_buf._odata[offset]
    }
    result+= ReconYp(Uchi);

    // Zp
    offset = Stencil._offsets [Zp][ss];
    local  = Stencil._is_local[Zp][ss];
    if ( local ) {
      Uchi = U[]*spProjZp(in._odata[offset]);
    } else {
      Uchi = U[]*comm_buf._odata[offset]
    }
    result+= ReconZp(Uchi);

    // Tp
    offset = Stencil._offsets [Tp][ss];
    local  = Stencil._is_local[Tp][ss];
    if ( local ) {
      Uchi = U[]*spProjTp(in._odata[offset]);
    } else {
      Uchi = U[]*comm_buf._odata[offset]
    }
    result+= ReconTp(Uchi);

    // Xm
    offset = Stencil._offsets [Xm][ss];
    local  = Stencil._is_local[Xm][ss];
    if ( local ) {
      Uchi = U[]*spProjXm(in._odata[offset]);
    } else {
      Uchi = U[]*comm_buf._odata[offset]
    }
    result+= ReconXm(Uchi);

    // Ym
    offset = Stencil._offsets [Ym][ss];
    local  = Stencil._is_local[Ym][ss];
    if ( local ) {
      Uchi = U[]*spProjYm(in._odata[offset]);
    } else {
      Uchi = U[]*comm_buf._odata[offset]
    }
    result+= ReconYm(Uchi);

    // Zm
    offset = Stencil._offsets [Zm][ss];
    local  = Stencil._is_local[Zm][ss];
    if ( local ) {
      Uchi = U[]*spProjZm(in._odata[offset]);
    } else {
      Uchi = U[]*comm_buf._odata[offset]
    }
    result+= ReconZm(Uchi);

    // Tm
    offset = Stencil._offsets [Tm][ss];
    local  = Stencil._is_local[Tm][ss];
    if ( local ) {
      Uchi = U[]*spProjTm(in._odata[offset]);
    } else {
      Uchi = U[]*comm_buf._odata[offset]
    }
    result+= ReconTm(Uchi);

    out._odata[ss] = result;
  }

}
void WilsonMatrix::Dw(const LatticeFermion &in, LatticeFermion &out)
{

}
void WilsonMatrix::MpcDag   (const LatticeFermion &in, LatticeFermion &out)
{

}
void WilsonMatrix::Mpc      (const LatticeFermion &in, LatticeFermion &out)
{

}
void WilsonMatrix::MpcDagMpc(const LatticeFermion &in, LatticeFermion &out)
{

}
void WilsonMatrix::MDagM    (const LatticeFermion &in, LatticeFermion &out)
{

}


}}
#endif
