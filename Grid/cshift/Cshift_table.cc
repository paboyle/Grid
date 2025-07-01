#include <Grid/GridCore.h>       
NAMESPACE_BEGIN(Grid);
std::vector<std::pair<int,int> > Cshift_table; 
deviceVector<std::pair<int,int> > Cshift_table_device; 
std::vector<int> Cshift_vector;
deviceVector<int> Cshift_vector_device; 
NAMESPACE_END(Grid);
