project(ME2CPP LANGUAGES CXX)
cmake_minimum_required(VERSION 3.12)
find_package(SHERPA-MC REQUIRED)
add_executable(test test.cpp)
target_link_libraries(test PRIVATE 
                           SHERPA::MEProcess 
                           SHERPA::ToolsOrg 
                           SHERPA::ToolsMath 
                           SHERPA::ToolsPhys 
                           SHERPA::SherpaMain 
                           SHERPA::PhasicProcess)
