# HDF5 libraries
find_package(HDF5)

if(!${HDF5_FOUND})
  message(FATAL_ERROR "${Red}The HDF5 libraries were not found. Please make sure you have HDF5 installed/loaded, and that the library and include directories can be found in, eg, CPLUS_INCLUDE_PATH and LD_LIBRARY_PATH.${ColorReset}")
else()
  message(STATUS " HDF5_C_LIBRARIES: ${HDF5_C_LIBRARIES}")
  set(HDF5_LINK_LIBRARY "${HDF5_LIBRARIES}")
  if(NOT DEFINED HDF5_C_LIBRARIES)
    message(STATUS " No HDF5_C_LIBRARIES. Trying HDF5_LIBRARIES: ${HDF5_LIBRARIES}")
    set(HDF5_LINK_LIBRARY "${HDF5_LIBRARIES}")
    if(NOT DEFINED HDF5_LIBRARIES)
      message(WARNING "${Yellow}The HDF5 libraries were not found. CMake will continue, but linking may not succeed. Please make sure you have HDF5 installed/loaded, and that the library and include directories can be found in, eg, CPLUS_INCLUDE_PATH and LD_LIBRARY_PATH.${ColorReset}")
    endif()
  endif()
  message(STATUS " HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")
  include_directories("${HDF5_INCLUDE_DIRS}")
endif()
