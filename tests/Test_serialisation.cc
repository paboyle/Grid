#include <Grid.h>

using namespace Grid;

class myclass {
public:

  GRID_DECL_CLASS_MEMBERS(myclass,
			  int, x,
			  double, y,
			  bool , b,
			  std::string, name,
			  std::vector<double>, array,
			  std::vector<std::vector<double> >, twodimarray,
			  );

  myclass(){}
  myclass(int i) : array(4,5.1), twodimarray(3,std::vector<double>(2,1.23456)) {
    x=i;
    y=2*i;
    b=true;
    name="bother said pooh";
  }

};


uint16_t i16 = 1;
uint16_t u16 = 2;
uint32_t i32 = 3;
uint32_t u32 = 4;
uint64_t i64 = 5;
uint64_t u64 = 6;
float    f = M_PI;
double   d = 2*M_PI;
bool     b = false;

int main(int argc,char **argv)
{
  {
    Writer WR("bother.xml");

    push(WR,"BasicTypes");
    write(WR,"i16",i16);
    write(WR,"u16",u16);
    write(WR,"i32",i32);
    write(WR,"i32",u32);
    write(WR,"i64",i64);
    write(WR,"i64",u64);
    write(WR,"f",f);
    write(WR,"d",d);
    write(WR,"b",b);
    pop(WR);

    myclass obj(1234); // non-trivial constructor
    write(WR,"obj",obj);

  };

  Reader RD("bother.xml");

  myclass copy1;
  myclass copy2;
  myclass copy3;

  read(RD,"obj",copy1);
  std::cout << "Loaded "  << copy1<<std::endl;

  {
    BinaryWriter BWR("bother.bin");
    write(BWR,"discard",copy1 );
  }
  { 
    BinaryReader BRD("bother.bin");
    read (BRD,"discard",copy2 );
    std::cout<<copy2<<std::endl;
  }


  {
    TextWriter TWR("bother.txt");
    write(TWR,"discard",copy1 );
  }
  { 
    TextReader TRD("bother.txt");
    read (TRD,"discard",copy3 );
    std::cout<<copy3<<std::endl;
  }
  
}
