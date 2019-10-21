cmake_minimum_required(VERSION 3.1)

project(pitchshifterpp)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake_modules)

if(UNIX)
    find_package(FFTW REQUIRED)
    find_package(LibSndFile REQUIRED)
elseif(WIN32)
endif()

include_directories(
    include
    include/pitchshifterpp
    )

if(WIN32)
    # even though this is a 3rd party lib, this include folder is required for sndfile.hh
    # to find sndfile.h with angled include for gcc (and mingw )
    include_directories(include/libsndfile)
endif()

# for qt creator to be able to find the files for autcomplete and so on
file(GLOB_RECURSE SOURCES *.cpp *.cc *.cxx)
file(GLOB_RECURSE HEADERS *.h *.hpp *.hh)

if(WIN32)
    link_directories(${CMAKE_SOURCE_DIR}/lib/fftw)
    link_directories(${CMAKE_SOURCE_DIR}/lib/libsndfile)
endif()

add_executable(${PROJECT_NAME} ${HEADERS} src/pitchshifterpp.cpp)


if(WIN32)
    target_link_libraries (${PROJECT_NAME}  libfftw3-3.dll)
    target_link_libraries (${PROJECT_NAME}  libsndfile-1.dll)
else()
    target_link_libraries(${PROJECT_NAME} fftw3)
    target_link_libraries(${PROJECT_NAME} sndfile)
endif()

if(WIN32)
    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${PROJECT_SOURCE_DIR}/lib/fftw/libfftw3-3.dll"
        $<TARGET_FILE_DIR:${PROJECT_NAME}>)

    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${PROJECT_SOURCE_DIR}/lib/libsndfile/libsndfile-1.dll"
        $<TARGET_FILE_DIR:${PROJECT_NAME}>)
endif()


add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${PROJECT_SOURCE_DIR}/sounds/sine1k.wav"
    $<TARGET_FILE_DIR:${PROJECT_NAME}>)

