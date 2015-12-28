#include <Grid.h>

namespace Grid {

GridStopWatch Logger::StopWatch;
std::ostream  Logger::devnull(0);

GridLogger GridLogError      (1,"Error");
GridLogger GridLogWarning    (1,"Warning");
GridLogger GridLogMessage    (1,"Message");
GridLogger GridLogDebug      (1,"Debug");
GridLogger GridLogPerformance(1,"Performance");
GridLogger GridLogIterative  (1,"Iterative");
GridLogger GridLogIntegrator (1,"Integrator");

void GridLogConfigure(std::vector<std::string> &logstreams)
{
  GridLogError.Active(0);
  GridLogWarning.Active(0);
  GridLogMessage.Active(0);
  GridLogIterative.Active(0);
  GridLogDebug.Active(0);
  GridLogPerformance.Active(0);
  GridLogIntegrator.Active(0);

  for(int i=0;i<logstreams.size();i++){
    if ( logstreams[i]== std::string("Error")       ) GridLogError.Active(1);
    if ( logstreams[i]== std::string("Warning")     ) GridLogWarning.Active(1);
    if ( logstreams[i]== std::string("Message")     ) GridLogMessage.Active(1);
    if ( logstreams[i]== std::string("Iterative")   ) GridLogIterative.Active(1);
    if ( logstreams[i]== std::string("Debug")       ) GridLogDebug.Active(1);
    if ( logstreams[i]== std::string("Performance") ) GridLogPerformance.Active(1);
    if ( logstreams[i]== std::string("Integrator" ) ) GridLogIntegrator.Active(1);
  }
}

////////////////////////////////////////////////////////////
// Verbose limiter on MPI tasks
////////////////////////////////////////////////////////////
void Grid_quiesce_nodes(void)
{
#ifdef GRID_COMMS_MPI
  int me;
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  if ( me ) { 
    std::cout.setstate(std::ios::badbit);
  }
#endif
}

void Grid_unquiesce_nodes(void)
{
#ifdef GRID_COMMS_MPI
    std::cout.clear();
#endif
}


}

