# CMakeLists.txt for test
#########################

# Include directory containing source
include_directories (../include)

# Create main lib (as advised in Catch docs, to speed up compilation)
add_library (TestMain
  testMain.cpp)

# Create executable to test io
add_executable (test_1_basic_ops.exe
  test_1_basic_ops.cpp)
target_link_libraries (test_1_basic_ops.exe
  GridSIMD
  TestMain)
