set(HEADERLIST

)

set(INCLUDEDIRS
  "${PROJECT_SOURCE_DIR}/headers"
)

set(SOURCELIST

    "${PROJECT_SOURCE_DIR}/src/main.cpp"

)

add_executable(TDSE main.cpp ${SOURCELIST})
target_include_directories(TDSE PUBLIC ${INCLUDEDIRS})
set_target_properties(TDSE PROPERTIES OUTPUT_NAME "simulation.exe")
target_link_libraries(TDSE PUBLIC PkgConfig::PETSC PkgConfig::SLEPC MPI::MPI_CXX)

