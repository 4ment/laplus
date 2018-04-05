cmake_minimum_required (VERSION 2.8)

project(liblesplace C)
set(CMAKE_C_FLAGS "-std=c99 -Wall -pedantic")
set(CMAKE_C_FLAGS_DEBUG "-g")

# Find GSL
include(FindGsl)
find_package(Gsl REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})

set(LESPLACE_LIB_C_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/lesplace.h)
set(LESPLACE_LIB_C_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/lesplace.c)
add_library(lesplace-static STATIC ${LESPLACE_LIB_C_FILES})
add_library(lesplace SHARED ${LESPLACE_LIB_C_FILES})
target_link_libraries(lesplace-static
  ${GSL_LIBRARIES})
target_link_libraries(lesplace
  ${GSL_LIBRARIES})

install(TARGETS lesplace DESTINATION lib)
install(FILES ${LESPLACE_LIB_C_HEADERS} DESTINATION include)
