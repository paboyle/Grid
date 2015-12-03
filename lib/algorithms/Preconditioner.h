#ifndef GRID_PRECONDITIONER_H
#define GRID_PRECONDITIONER_H

namespace Grid {

  template<class Field> class Preconditioner :  public LinearFunction<Field> { 
    virtual void operator()(const Field &src, Field & psi)=0;
  };

  template<class Field> class TrivialPrecon :  public Preconditioner<Field> { 
  public:
    void operator()(const Field &src, Field & psi){
      psi = src;
    }
    TrivialPrecon(void){};
  };

}
#endif
