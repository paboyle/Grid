#include <Grid/GridCore.h>

int                    Grid::BinaryIO::latticeWriteMaxRetry = -1;
Grid::BinaryIO::IoPerf Grid::BinaryIO::lastPerf;

// Target size of a single contiguous file extent under BINARYIO_AGGREGATE.
// 4MB is around the knee for Lustre; exposed so it can be swept at runtime.
uint64_t               Grid::BinaryIO::aggregateTargetBytes = 4*1024*1024;
