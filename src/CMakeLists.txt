cmake_minimum_required (VERSION 2.8)

project(liblaplus C)
set(CMAKE_C_FLAGS "-std=c99 -Wall -pedantic")
set(CMAKE_C_FLAGS_DEBUG "-g")

# Find GSL
include(FindGSL)
find_package(GSL REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})

set(LAPLUS_LIB_C_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/laplus.h)
set(LAPLUS_LIB_C_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/laplus.c)
add_library(laplus-static STATIC ${LAPLUS_LIB_C_FILES})
add_library(laplus SHARED ${LAPLUS_LIB_C_FILES})
target_link_libraries(laplus-static
  ${GSL_LIBRARIES})
target_link_libraries(laplus
  ${GSL_LIBRARIES})

install(TARGETS laplus DESTINATION lib)
install(FILES ${LAPLUS_LIB_C_HEADERS} DESTINATION include)
