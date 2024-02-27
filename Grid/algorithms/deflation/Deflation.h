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
class ZeroGuesser: public LinearFunction<Field> {
public:
  using LinearFunction<Field>::operator();
    virtual void operator()(const Field &src, Field &guess) { guess = Zero(); };
};
template<class Field>
class DoNothingGuesser: public LinearFunction<Field> {
public:
  using LinearFunction<Field>::operator();
  virtual void operator()(const Field &src, Field &guess) {  };
};
template<class Field>
class SourceGuesser: public LinearFunction<Field> {
public:
  using LinearFunction<Field>::operator();
  virtual void operator()(const Field &src, Field &guess) { guess = src; };
};

////////////////////////////////
// Fine grid deflation
////////////////////////////////
template<class Field>
class DeflatedGuesser: public LinearFunction<Field> {
private:
  const std::vector<Field> &evec;
  const std::vector<RealD> &eval;
  const unsigned int       N;

public:
  using LinearFunction<Field>::operator();

  DeflatedGuesser(const std::vector<Field> & _evec,const std::vector<RealD> & _eval)
  : DeflatedGuesser(_evec, _eval, _evec.size())
  {}

  DeflatedGuesser(const std::vector<Field> & _evec, const std::vector<RealD> & _eval, const unsigned int _N)
  : evec(_evec), eval(_eval), N(_N)
  {
    assert(evec.size()==eval.size());
    assert(N <= evec.size());
  } 

  virtual void operator()(const Field &src,Field &guess) {
    guess = Zero();
    for (int i=0;i<N;i++) {
      const Field& tmp = evec[i];
      axpy(guess,TensorRemove(innerProduct(tmp,src)) / eval[i],tmp,guess);
    }
    guess.Checkerboard() = src.Checkerboard();
  }
};

template<class FineField, class CoarseField>
class LocalCoherenceDeflatedGuesser: public LinearFunction<FineField> {
private:
  const std::vector<FineField>   &subspace;
  const std::vector<CoarseField> &evec_coarse;
  const std::vector<RealD>       &eval_coarse;
public:
  
  using LinearFunction<FineField>::operator();
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
    CoarseField src_coarse(evec_coarse[0].Grid());
    CoarseField guess_coarse(evec_coarse[0].Grid());    guess_coarse = Zero();
    blockProject(src_coarse,src,subspace);    
    for (int i=0;i<N;i++) {
      const CoarseField & tmp = evec_coarse[i];
      axpy(guess_coarse,TensorRemove(innerProduct(tmp,src_coarse)) / eval_coarse[i],tmp,guess_coarse);
    }
    blockPromote(guess_coarse,guess,subspace);
    guess.Checkerboard() = src.Checkerboard();
  };

  void operator()(const std::vector<FineField> &src,std::vector<FineField> &guess) {
    int Nevec = (int)evec_coarse.size();
    int Nsrc = (int)src.size();
    // make temp variables
    std::vector<CoarseField> src_coarse(Nsrc,evec_coarse[0].Grid());
    std::vector<CoarseField> guess_coarse(Nsrc,evec_coarse[0].Grid());    
    //Preporcessing
    std::cout << GridLogMessage << "Start BlockProject for loop" << std::endl;
    for (int j=0;j<Nsrc;j++)
    {
    guess_coarse[j] = Zero();
    std::cout << GridLogMessage << "BlockProject iter: " << j << std::endl;
    blockProject(src_coarse[j],src[j],subspace);
    }
    //deflation set up for eigen vector batchsize 1 and source batch size equal number of sources
    std::cout << GridLogMessage << "Start ProjectAccum for loop" << std::endl;
    for (int i=0;i<Nevec;i++)
    {
      std::cout << GridLogMessage << "ProjectAccum Nvec: " << i << std::endl;
      const CoarseField & tmp = evec_coarse[i];
      for (int j=0;j<Nsrc;j++)
      {
        axpy(guess_coarse[j],TensorRemove(innerProduct(tmp,src_coarse[j])) / eval_coarse[i],tmp,guess_coarse[j]);
      }
    }
    //postprocessing
    std::cout << GridLogMessage << "Start BlockPromote for loop" << std::endl;
    for (int j=0;j<Nsrc;j++)
    {
    std::cout << GridLogMessage << "BlockProject iter: " << j << std::endl;
    blockPromote(guess_coarse[j],guess[j],subspace);
    guess[j].Checkerboard() = src[j].Checkerboard();
    }
  };

  };



}
#endif
