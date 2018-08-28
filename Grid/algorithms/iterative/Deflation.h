    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/ImplicitlyRestartedLanczos.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#ifndef GRID_DEFLATION_H
#define GRID_DEFLATION_H

namespace Grid { 

template<class Field>
class Guesser {
public:
  Guesser(void) = default;
  virtual ~Guesser(void) = default;
  virtual void operator()(const Field &src, Field &guess) = 0;
};

template<class Field>
class ZeroGuesser: public Guesser<Field> {
public:
  virtual void operator()(const Field &src, Field &guess) { guess = zero; };
};

template<class Field>
class SourceGuesser: public Guesser<Field> {
public:
  virtual void operator()(const Field &src, Field &guess) { guess = src; };
};

////////////////////////////////
// Fine grid deflation
////////////////////////////////
template<class Field>
class DeflatedGuesser: public Guesser<Field> {
private:
  const std::vector<Field> &evec;
  const std::vector<RealD> &eval;

public:

  DeflatedGuesser(const std::vector<Field> & _evec,const std::vector<RealD> & _eval) : evec(_evec), eval(_eval) {};

  virtual void operator()(const Field &src,Field &guess) {
    guess = zero;
    assert(evec.size()==eval.size());
    auto N = evec.size();
    for (int i=0;i<N;i++) {
      const Field& tmp = evec[i];
      axpy(guess,TensorRemove(innerProduct(tmp,src)) / eval[i],tmp,guess);
    }
    guess.checkerboard = src.checkerboard;
  }
};

template<class FineField, class CoarseField>
class LocalCoherenceDeflatedGuesser: public Guesser<FineField> {
private:
  const std::vector<FineField>   &subspace;
  const std::vector<CoarseField> &evec_coarse;
  const std::vector<RealD>       &eval_coarse;
public:
  
  LocalCoherenceDeflatedGuesser(const std::vector<FineField>   &_subspace,
				const std::vector<CoarseField> &_evec_coarse,
				const std::vector<RealD>       &_eval_coarse)
    : subspace(_subspace), 
      evec_coarse(_evec_coarse), 
      eval_coarse(_eval_coarse)  
  {
  }
  
  void operator()(const FineField &src,FineField &guess) { 
    int N = (int)evec_coarse.size();
    CoarseField src_coarse(evec_coarse[0]._grid);
    CoarseField guess_coarse(evec_coarse[0]._grid);    guess_coarse = zero;
    blockProject(src_coarse,src,subspace);    
    for (int i=0;i<N;i++) {
      const CoarseField & tmp = evec_coarse[i];
      axpy(guess_coarse,TensorRemove(innerProduct(tmp,src_coarse)) / eval_coarse[i],tmp,guess_coarse);
    }
    blockPromote(guess_coarse,guess,subspace);
    guess.checkerboard = src.checkerboard;
  };
};



}
#endif
