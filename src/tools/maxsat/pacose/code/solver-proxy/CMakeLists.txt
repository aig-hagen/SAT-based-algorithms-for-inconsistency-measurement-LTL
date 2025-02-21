#cmake_minimum_required(VERSION 3.10)                 # CMake version check
project("SATSolverProxy")                                      # Create project "simple_example"

set(CMAKE_CXX_STANDARD 11)                           # Enable c++11 standard
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -g") #-g produce debug information, g3 additional debug information, -Wall activates Warnings, -Wextra activates additional warnings!
#set(CMAKE_EXE_LINKER_FLAGS "-static -Wall -Weffc++ -Wextra -Wswitch -Wsign-compare ${CMAKE_EXE_LINKER_FLAGS}" )
#set(CMAKE_C_FLAGS "-O3 -Wall -Weffc++ -Wextra -Wswitch -Wsign-compare ${CMAKE_C_FLAGS}")
#SET(Boost_USE_STATIC_LIBS ON)
#set (EXTRA_OPTIONS "-O3  -Wno-deprecated -fPIC -DLINUX")
#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-no-pie ") #-no-pie only for ANTOM needed
set (EXTRA_OPTIONS "-O3")

#set(CMAKE_BUILD_TYPE "Debug")                          #Possible values are empty, Debug, Release, RelWithDebInfo and MinSizeRel.
set(CMAKE_BUILD_TYPE Release)

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include_directories("$(PROJECT_SOURCE_DIR)")

link_directories(${PROJECT_SOURCE_DIR}/cryptoMiniSat565)
link_directories(${PROJECT_SOURCE_DIR}/glucose421)
link_directories(${PROJECT_SOURCE_DIR}/glucose3)
#link_directories(${PROJECT_SOURCE_DIR}/mapleLCMDistChronoBT)
#link_directories(${PROJECT_SOURCE_DIR}/antom)
link_directories(${PROJECT_SOURCE_DIR}/cadical)


set(SOURCE_FILES
    Main.cpp
    SATSolverProxy.cpp
#    MapleGlucoseSolverProxy.cpp
    Glucose421SolverProxy.cpp
    Glucose3SolverProxy.cpp
    CryptoMiniSat565SolverProxy.cpp
#    AntomSolverProxy.cpp
    CadicalSolverProxy.cpp
    ) # Add main.cpp file of project root directory as source file

add_library(SSolverProxy STATIC ${SOURCE_FILES})
add_executable(${PROJECT_NAME} ${SOURCE_FILES})               # Add executable target with source files listed in SOURCE_FILES variable

set(LIBS
    Glucose421
    Glucose3
    cryptominisat5
#    MapleGlucose
#    antomNDEBUG
    cadical
    pthread      #for ANTOM and CMS needed
#    rt          #only for ANTOM needed
    m4ri        #only for CMS static version needed
    )

target_link_libraries(${PROJECT_NAME} ${LIBS})
target_link_libraries(SSolverProxy ${LIBS})
