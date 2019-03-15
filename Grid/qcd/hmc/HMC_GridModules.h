/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/HMC_GridModules.h

Copyright (C) 2015
Copyright (C) 2016

Author: Guido Cossu <guido.cossu@ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef HMC_GRID_MODULES
#define HMC_GRID_MODULES

namespace Grid {

// Resources
// Modules for grids

// Introduce another namespace HMCModules?

class GridModuleParameters: Serializable{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(GridModuleParameters,
  std::string, lattice,
  std::string, mpi);

  std::vector<int> getLattice() const {return strToVec<int>(lattice);}
  std::vector<int> getMpi()     const {return strToVec<int>(mpi);}


  void check() const {
    if (getLattice().size() != getMpi().size() ) {
      std::cout << GridLogError
                << "Error in GridModuleParameters: lattice and mpi dimensions "
                   "do not match"
                << std::endl;
      exit(1);
    }
  }

  template <class ReaderClass>
  GridModuleParameters(Reader<ReaderClass>& Reader, std::string n = "LatticeGrid"):name(n) {
    read(Reader, name, *this);
    check();
  }

  // Save on file
  template< class WriterClass>
  void save(Writer<WriterClass>& Writer){
    check();
    write(Writer, name, *this);
  }
private:
    std::string name;
};

// Lower level class
class GridModule {
 public:
  GridCartesian* get_full() {
    std::cout << GridLogDebug << "Getting cartesian in module"<< std::endl;
    return grid_.get(); }
  GridRedBlackCartesian* get_rb() {
    std::cout << GridLogDebug << "Getting rb-cartesian in module"<< std::endl;
    return rbgrid_.get(); }

  void set_full(GridCartesian* grid) { grid_.reset(grid); }
  void set_rb(GridRedBlackCartesian* rbgrid) { rbgrid_.reset(rbgrid); }
  void show_full_decomposition(){ grid_->show_decomposition(); }
  void show_rb_decomposition(){ rbgrid_->show_decomposition(); }

 protected:
  std::unique_ptr<GridCartesian> grid_;
  std::unique_ptr<GridRedBlackCartesian> rbgrid_;

};

////////////////////////////////////
// Classes for the user
////////////////////////////////////
// Note: the space time grid should be out of the QCD namespace
template <class vector_type>
class GridFourDimModule : public GridModule
{
public:
  GridFourDimModule()
  {
    using namespace QCD;
    set_full(SpaceTimeGrid::makeFourDimGrid(
        GridDefaultLatt(), 
        GridDefaultSimd(4, vector_type::Nsimd()),
        GridDefaultMpi()));
    set_rb(SpaceTimeGrid::makeFourDimRedBlackGrid(grid_.get()));
  }

  GridFourDimModule(const std::vector<int> tweak_simd)
  {
    using namespace QCD;
    if (tweak_simd.size() != 4)
    {
      std::cout << GridLogError
                << "Error in GridFourDimModule: SIMD size different from 4" 
                << std::endl;
      exit(1);
    }

    // Checks that the product agrees with the expectation
    int simd_sum = 1;
    for (auto &n : tweak_simd)
      simd_sum *= n;
    std::cout << GridLogDebug << "TweakSIMD: " << tweak_simd << "  Sum: " << simd_sum << std::endl;

    if (simd_sum == vector_type::Nsimd())
    {
      set_full(SpaceTimeGrid::makeFourDimGrid(
          GridDefaultLatt(), 
          tweak_simd, 
          GridDefaultMpi()));
      set_rb(SpaceTimeGrid::makeFourDimRedBlackGrid(grid_.get()));
    }
    else
    {
      std::cout << GridLogError 
                << "Error in GridFourDimModule: SIMD lanes must sum to " 
                << vector_type::Nsimd() 
                << std::endl;
    }
  }

  GridFourDimModule(const GridModuleParameters Params)
  {
    using namespace QCD;
    std::vector<int> lattice_v = Params.getLattice();
    std::vector<int> mpi_v = Params.getMpi();
    if (lattice_v.size() == 4)
    {
      set_full(SpaceTimeGrid::makeFourDimGrid(
          lattice_v, 
          GridDefaultSimd(4, vector_type::Nsimd()),
          mpi_v));
      set_rb(SpaceTimeGrid::makeFourDimRedBlackGrid(grid_.get()));
    }
    else
    {
      std::cout << GridLogError
                << "Error in GridFourDimModule: lattice dimension different from 4"
                << std::endl;
      exit(1);
    }
  }
};

typedef GridFourDimModule<vComplex> GridDefaultFourDimModule;


}  // namespace Grid

#endif  // HMC_GRID_MODULES
