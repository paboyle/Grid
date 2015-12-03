#ifndef GRID_INIT_H
#define GRID_INIT_H

namespace Grid {

  void Grid_init(int *argc,char ***argv);
  void Grid_finalize(void);
  // internal, controled with --handle
  void Grid_sa_signal_handler(int sig,siginfo_t *si,void * ptr);
  void Grid_debug_handler_init(void);
  void Grid_quiesce_nodes(void);
  void Grid_unquiesce_nodes(void);

  const std::vector<int> GridDefaultSimd(int dims,int nsimd);
  const std::vector<int> &GridDefaultLatt(void);
  const std::vector<int> &GridDefaultMpi(void);
  const int              &GridThreads(void)  ;
  void                    GridSetThreads(int t) ;

  // Common parsing chores
  std::string GridCmdOptionPayload(char ** begin, char ** end, const std::string & option);
  bool        GridCmdOptionExists(char** begin, char** end, const std::string& option);
  std::string GridCmdVectorIntToString(const std::vector<int> & vec);

  void GridParseLayout(char **argv,int argc,
		       std::vector<int> &latt,
		       std::vector<int> &simd,
		       std::vector<int> &mpi);


};
#endif
