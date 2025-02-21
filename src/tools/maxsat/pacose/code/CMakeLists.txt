cmake_minimum_required(VERSION 3.10)                 # CMake version check
project("Pacose")

set(CMAKE_CXX_STANDARD 11)                           # Enable c++11 standard
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -g") #-g produce debug information, g3 additional debug information, -Wall activates Warnings, -Wextra activates additional warnings!
#set(CMAKE_EXE_LINKER_FLAGS "-static -Wall -Weffc++ -Wextra -Wswitch -Wsign-compare ${CMAKE_EXE_LINKER_FLAGS}" )
#set(CMAKE_C_FLAGS "-O3 -Wall -Weffc++ -Wextra -Wswitch -Wsign-compare ${CMAKE_C_FLAGS}")
#set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wswitch -Wsign-compare ${CMAKE_CXX_FLAGS}")
#set(CMAKE_CXX_FLAGS "-pg")

#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-no-pie ") #-no-pie only for ANTOM needed
set (EXTRA_OPTIONS "-O3")
#set (EXTRA_OPTIONS "-fPIC -no-pie")

#set(CMAKE_BUILD_TYPE "Debug")                          #Possible values are empty, Debug, Release, RelWithDebInfo and MinSizeRel.
set(CMAKE_BUILD_TYPE "Release")

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

add_subdirectory(solver-proxy)

include_directories(
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/maxSAT
    ${PROJECT_SOURCE_DIR}/maxSAT/DGPW
    ${PROJECT_SOURCE_DIR}/solver-proxy
    )

link_directories(${CMAKE_SOURCE_DIR}/solver-proxy/bin)
link_directories(${CMAKE_SOURCE_DIR}/solver-proxy/cryptoMiniSat565)
link_directories(${CMAKE_SOURCE_DIR}/solver-proxy/glucose421)
link_directories(${CMAKE_SOURCE_DIR}/solver-proxy/glucose3)
#link_directories(${CMAKE_SOURCE_DIR}/solver-proxy/mapleLCMDistChronoBT)
#link_directories(${CMAKE_SOURCE_DIR}/solver-proxy/antom)
link_directories(${CMAKE_SOURCE_DIR}/solver-proxy/cadical)
link_directories(${CMAKE_SOURCE_DIR}/bin)

set(SOURCE_FILES
    Main.cpp
    maxSAT/Pacose.cpp
    maxSAT/Encodings.cpp
    maxSAT/Greedyprepro.cpp
    maxSAT/DGPW/bucket.cpp
    maxSAT/DGPW/cascade.cpp
    maxSAT/DGPW/dgpw.cpp
    maxSAT/DGPW/multiplecascade.cpp
    maxSAT/DGPW/sorter.cpp
    ### NADJA cpp file
    ) 

#add_library(Pac STATIC ${SOURCE_FILES})

set(LIBS
    SSolverProxy
    z
    ### NADJA LIBRARY
#    pthread     #for ANTOM and CMS needed
#    rt          #only for ANTOM needed
#    m4ri        #only for CMS static version needed
    )

add_executable(${PROJECT_NAME} ${SOURCE_FILES})               # Add executable target with source files listed in SOURCE_FILES variable
target_link_libraries(${PROJECT_NAME} ${LIBS}  "-static")
#target_link_libraries(Pac ${LIBS})
