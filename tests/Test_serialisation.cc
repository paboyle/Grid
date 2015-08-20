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

  myclass() : array(4,5.0), twodimarray(3,std::vector<double>(2,2.0)) {
    x=1;
    y=2;
    b=false;
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

    myclass obj;
    write(WR,"obj",obj);

  };

  Reader RD("bother2.xml");

  myclass copy;

  read(RD,"obj",copy);
  std::cout << "Loaded "  << copy<<std::endl;
}
