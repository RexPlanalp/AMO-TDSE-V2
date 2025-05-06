set(HEADERLIST

  "${PROJECT_SOURCE_DIR}/include/Angular.h"
  "${PROJECT_SOURCE_DIR}/include/Box.h"
  "${PROJECT_SOURCE_DIR}/include/misc.h"
  "${PROJECT_SOURCE_DIR}/include/Laser.h"
)

set(INCLUDEDIRS
  "${PROJECT_SOURCE_DIR}/include"
)

set(SOURCELIST

    "${PROJECT_SOURCE_DIR}/src/main.cpp"
    "${PROJECT_SOURCE_DIR}/src/Box.cpp"
    "${PROJECT_SOURCE_DIR}/src/Angular.cpp"
    "${PROJECT_SOURCE_DIR}/src/Laser.cpp"

)

add_executable(TDSE main.cpp ${SOURCELIST})
target_include_directories(TDSE PUBLIC ${INCLUDEDIRS})
set_target_properties(TDSE PROPERTIES OUTPUT_NAME "simulation.exe")
target_link_libraries(TDSE PUBLIC PkgConfig::PETSC PkgConfig::SLEPC MPI::MPI_CXX GSL::gsl GSL::gslcblas nlohmann_json::nlohmann_json)

