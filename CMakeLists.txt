cmake_minimum_required(VERSION 3.20)
project(
  edgeDeletionsHarmonicResistance
  VERSION 0.1.0
  LANGUAGES CXX)

# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_EXTENSIONS OFF)

  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()

include(FetchContent)

# networkit fork / local modifications in extern/networkit (remote is bernlu/networkit)
add_subdirectory(extern/networkit)

find_package(OpenMP REQUIRED)

# The library code is here
add_subdirectory(src)

# The executable code is here
add_subdirectory(apps)
