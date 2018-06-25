include_directories(.)
set(LIBS Su2Theta ${CMAKE_THREAD_LIBS_INIT} )
LIST(APPEND LIBS ${MPI_CXX_LIBRARIES})
LIST(APPEND LIBS ${HDF5_LIBRARY})

if(USE_GPU)
  LIST(APPEND LIBS ${CUDA_cufft_LIBRARY} ${CUDA_cudart_LIBRARY})
endif()

add_executable(Su2Test Su2Test.cpp)

target_link_libraries(Su2Test ${LIBS} )
