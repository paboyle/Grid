/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/algorithms//CoarseningLookupTable.h

    Copyright (C) 2015 - 2020

Author: Daniel Richtmann <daniel.richtmann@ur.de>

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

    See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#pragma once

NAMESPACE_BEGIN(Grid);

class CoarseningLookupTable {
public:

  /////////////////////////////////////////////
  // Type Definitions
  /////////////////////////////////////////////

  typedef uint64_t index_type;
  typedef uint64_t size_type;

  /////////////////////////////////////////////
  // Member Data
  /////////////////////////////////////////////

private:
  GridBase*                       coarse_;
  GridBase*                       fine_;
  bool                            isPopulated_;
  std::vector<Vector<index_type>> lutVec_;
  Vector<index_type*>             lutPtr_;
  Vector<size_type>               sizes_;
  Vector<index_type>              reverseLutVec_;

  /////////////////////////////////////////////
  // Member Functions
  /////////////////////////////////////////////

public:
  CoarseningLookupTable(GridBase* coarse, GridBase* fine)
    : coarse_(coarse)
    , fine_(fine)
    , isPopulated_(false)
    , lutVec_(coarse_->oSites())
    , lutPtr_(coarse_->oSites())
    , sizes_(coarse_->oSites())
    , reverseLutVec_(fine_->oSites()) {
    populate(coarse_, fine_);
  }

  template<class ScalarField,
           typename std::enable_if<is_lattice<ScalarField>::value, ScalarField>::type * = nullptr>
  CoarseningLookupTable(GridBase* coarse, ScalarField const& mask)
    : coarse_(coarse)
    , fine_(mask.Grid())
    , isPopulated_(false)
    , lutVec_(coarse_->oSites())
    , lutPtr_(coarse_->oSites())
    , sizes_(coarse_->oSites())
    , reverseLutVec_(fine_->oSites()){
    populate(coarse_, mask);
  }

  CoarseningLookupTable()
    : coarse_(nullptr)
    , fine_(nullptr)
    , isPopulated_(false)
    , lutVec_()
    , lutPtr_()
    , sizes_()
    , reverseLutVec_()
  {}

  // clang-format off
  accelerator_inline std::vector<Vector<index_type>> const& operator()()  const { return lutVec_; }     // CPU access (TODO: remove?)
  accelerator_inline index_type const* const*               View()        const { return &lutPtr_[0]; } // GPU access
  accelerator_inline size_type  const*                      Sizes()       const { return &sizes_[0]; }  // also needed for GPU access
  accelerator_inline index_type const*                      ReverseView() const { return &reverseLutVec_[0]; }
  // clang-format on

  bool isPopulated() const { return isPopulated_; }

  bool gridPointersMatch(GridBase* coarse, GridBase* fine) const {
    // NOTE: This is the same check that "conformable" does
    return (coarse == coarse_) && (fine == fine_);
  }

  void setGridPointers(GridBase* coarse, GridBase* fine) {
    coarse_ = coarse;
    fine_   = fine;
  }

  void populate(GridBase* coarse, GridBase* fine) {
    Lattice<iScalar<vComplex>> fullmask(fine);
    fullmask = 1.;
    populate(coarse, fullmask);
  }

  template<class ScalarField>
  void populate(GridBase* coarse, ScalarField const& mask) {
    if(!gridPointersMatch(coarse, mask.Grid())) {
      setGridPointers(coarse, mask.Grid());
    }

    int        _ndimension = coarse_->_ndimension;
    Coordinate block_r(_ndimension);

    size_type block_v = 1;
    for(int d = 0; d < _ndimension; ++d) {
      block_r[d] = fine_->_rdimensions[d] / coarse_->_rdimensions[d];
      assert(block_r[d] * coarse_->_rdimensions[d] == fine_->_rdimensions[d]);
      block_v *= block_r[d];
    }
    assert(block_v == fine_->oSites()/coarse_->oSites());

    lutVec_.resize(coarse_->oSites());
    lutPtr_.resize(coarse_->oSites());
    sizes_.resize(coarse_->oSites());
    reverseLutVec_.resize(fine_->oSites());
    for(index_type sc = 0; sc < coarse_->oSites(); ++sc) {
      lutVec_[sc].resize(0);
      lutVec_[sc].reserve(block_v);
      lutPtr_[sc] = &lutVec_[sc][0];
      sizes_[sc]  = 0;
    }

    typename ScalarField::scalar_type zz = {0., 0.,};

    auto mask_v = mask.View();
    thread_for(sc, coarse_->oSites(), {
      Coordinate coor_c(_ndimension);
      Lexicographic::CoorFromIndex(coor_c, sc, coarse_->_rdimensions);

      int sf_tmp, count = 0;
      for(int sb = 0; sb < block_v; ++sb) {
        Coordinate coor_b(_ndimension);
        Coordinate coor_f(_ndimension);

        Lexicographic::CoorFromIndex(coor_b, sb, block_r);
        for(int d = 0; d < _ndimension; ++d) coor_f[d] = coor_c[d] * block_r[d] + coor_b[d];
        Lexicographic::IndexFromCoor(coor_f, sf_tmp, fine_->_rdimensions);

        index_type sf = (index_type)sf_tmp;

        if(Reduce(TensorRemove(coalescedRead(mask_v[sf]))) != zz) {
          lutPtr_[sc][count] = sf;
          sizes_[sc]++;
          count++;
        }
        reverseLutVec_[sf] = sc; // reverse table will never have holes
      }
      lutVec_[sc].resize(sizes_[sc]);
    });

    isPopulated_ = true;

    std::cout << GridLogMessage << "Recalculation of coarsening lookup table finished" << std::endl;
  }
};

NAMESPACE_END(Grid);
