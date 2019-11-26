cmake_minimum_required(VERSION 3.1)

project(pitchshifterpp)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake_modules)

find_package(FFTW REQUIRED)
find_package(LibSndFile REQUIRED)

include_directories(
    include
    include/pitchshifterpp
    )

# for qt creator to be able to find the files for autcomplete and so on
file(GLOB_RECURSE SOURCES *.cpp *.cc *.cxx)
file(GLOB_RECURSE HEADERS *.h *.hpp *.hh)

add_executable(${PROJECT_NAME} ${HEADERS} src/pitchshifterpp.cpp)
target_link_libraries(${PROJECT_NAME} fftw3)
target_link_libraries(${PROJECT_NAME} sndfile)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${PROJECT_SOURCE_DIR}/sounds/sine1k.wav"
    $<TARGET_FILE_DIR:${PROJECT_NAME}>)

