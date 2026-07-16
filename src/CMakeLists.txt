# src/CMakeLists.txt

# ── Generation of various libraries ──

# Basic utility library
add_library(nagamu_util SHARED
  ns_seed.cpp
  ns_myapp.cpp
  ns_angle_util.cpp
  ns_tuple_int.cpp
  ns_io_binary.cpp
  ns_mymacro.cpp
  ns_mylogger.cpp
  ns_param_constants.cpp
  cls_Angle.cpp
  cls_DepthResolutionSweep.cpp
  cls_PrefixSum2D.cpp
  ns_eigen_blas.cpp
  ns_detector_indexing.cpp
  ns_mix_error.cpp
  ns_iodir.cpp
  ns_stats_util.cpp
  ns_geom_util.cpp
  st_SignalNoiseStat.cpp
)

# Grid library
add_library(nagamu_grid SHARED
  cls_Grid1d.cpp
  cls_Grid2d.cpp
  cls_Grid3d.cpp
  cls_Grid1dXZ.cpp
)

# Terrain library
# NOTE: cls_Cuboid / cls_Grid2dCuboid* sources were renamed to cls_Pillar /
#       cls_Grid2dPillar* (the type models a vertical pillar, not a generic box).
add_library(nagamu_geom SHARED
  cls_Pillar.cpp
  cls_Voxel.cpp
  cls_Grid2dPillar.cpp
  cls_Grid2dPillarParameters.cpp
  cls_Grid2dVoxel.cpp
  cls_Grid3dVoxel.cpp
  cls_Grid3dVoxelParameters.cpp
  cls_AABB.cpp
  ns_dem_operator.cpp
)

# Detector library
add_library(nagamu_detector SHARED
  cls_Ray.cpp
  cls_DetectorElement.cpp
  cls_Grid2dBinGroup.cpp
  cls_Grid2dBinGroupParameters.cpp
  cls_DetectorPanel.cpp
  cls_DetectorPanelParameters.cpp
  cls_DetectorPanelArray.cpp
  cls_EfficiencyModel.cpp
)

# Path length calculation library
add_library(nagamu_path SHARED
  ns_pathcalc.cpp
  ns_calc_dNdD.cpp
  cls_PathLengthParameters.cpp
)

# Data table library
add_library(nagamu_data SHARED
  cls_Grid2dXYZ.cpp
  cls_FluxTable.cpp
  cls_FluxTableBuilder.cpp
  cls_NoiseParameters.cpp
)

# Map storage library
add_library(nagamu_map SHARED
  cls_VoxelUniqueIndexMapContainer.cpp
  cls_UqidManager.cpp
  cls_GroupManager.cpp
  cls_UqigAvailIndexer.cpp
  cls_DetectorIndexContainer.cpp
)

# Reconstruction library
add_library(nagamu_recon SHARED
  cls_NagaInv.cpp
  st_NagaInvLooper.cpp
  cls_NagaInvParameters.cpp
  cls_NagaInvManager.cpp
)

# cpp for modularizing the main execution function
add_library(nagamu_exec_module SHARED
  exemdl_init_json_logger_seed.cpp
  exemdl_load_parameters.cpp
  exemdl_build_geometry.cpp
  exemdl_build_detector.cpp
  exemdl_build_prior.cpp
  exemdl_calc_matrix.cpp
  exemdl_run_inversion.cpp
  exemdl_calc_rec_error.cpp
  exemdl_pipeline.cpp
  exemdl_sweep.cpp
)

# Other libraries
add_library(nagamu_etc SHARED
  cls_XYArray.cpp
  cls_RectAngularBinGroup.cpp
)

target_compile_definitions(nagamu_util PUBLIC SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE)
target_compile_definitions(nagamu_grid PUBLIC SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE)
target_compile_definitions(nagamu_geom PUBLIC SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE)
target_compile_definitions(nagamu_detector PUBLIC SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE)
target_compile_definitions(nagamu_path PUBLIC SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE)
target_compile_definitions(nagamu_data PUBLIC SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE)
target_compile_definitions(nagamu_map PUBLIC SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE)
target_compile_definitions(nagamu_recon PUBLIC SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE)
target_compile_definitions(nagamu_exec_module PUBLIC SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE)
target_compile_definitions(nagamu_etc PUBLIC SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE)

# ── Include directory settings for each library ──
foreach(lib nagamu_util nagamu_grid nagamu_geom nagamu_detector
          nagamu_path nagamu_data nagamu_map nagamu_recon nagamu_etc
          nagamu_exec_module )
  target_include_directories(${lib} PUBLIC
    ${PROJECT_SOURCE_DIR}/include
  )
  target_link_libraries(${lib} PUBLIC
    Eigen3::Eigen
    nlohmann_json::nlohmann_json
    fmt::fmt
    spdlog::spdlog
  )
endforeach()

# ── OpenMP link settings ──
foreach(lib nagamu_util nagamu_grid nagamu_geom nagamu_detector
              nagamu_path nagamu_data nagamu_map nagamu_recon nagamu_etc)
  target_link_libraries(${lib} PUBLIC OpenMP::OpenMP_CXX)
endforeach()

# ── BLAS/LAPACK link settings ──
target_link_libraries(nagamu_util PUBLIC ${BLAS_LINK_LIBRARIES})

# ── Dependencies between shared libraries ──
# nagamu_exec_module uses symbols from nagamu_recon (NagaInv, NagaInvManager, etc.)
target_link_libraries(nagamu_exec_module PUBLIC nagamu_recon)

# ── Generation of libmuog (changed to a STATIC library) ──
# * libmuog.cpp may be empty (dummy source)
add_library(libmuog STATIC libmuog.cpp)

target_link_libraries(libmuog
  PUBLIC
    nagamu_util
    nagamu_grid
    nagamu_geom
    nagamu_detector
    nagamu_path
    nagamu_data
    nagamu_map
    nagamu_recon
    nagamu_etc
    nagamu_exec_module
    Eigen3::Eigen
    nlohmann_json::nlohmann_json
    fmt::fmt
    spdlog::spdlog
    OpenMP::OpenMP_CXX
    ${BLAS_LINK_LIBRARIES}
)

target_include_directories(libmuog PUBLIC
  ${PROJECT_SOURCE_DIR}/include
)

# ── libmuog compile options/definitions ──
# IMPORTANT: EIGEN_DONT_ALIGN_STATICALLY must be consistent across ALL
# libraries (shared + static) to avoid ODR violation causing struct layout
# mismatch between .so and executable.  Previously only libmuog had this
# definition, so shared libraries used different Eigen alignment — resulting
# in corrupted AppParameters members at runtime (heap corruption / SEGV).
if(CMAKE_BUILD_TYPE MATCHES Debug)
  foreach(lib nagamu_util nagamu_grid nagamu_geom nagamu_detector
            nagamu_path nagamu_data nagamu_map nagamu_recon nagamu_etc
            nagamu_exec_module libmuog)
    target_compile_definitions(${lib} PUBLIC DEBUG)
    target_compile_definitions(${lib} PUBLIC EIGEN_DONT_ALIGN_STATICALLY)
  endforeach()
endif()

if(CMAKE_BUILD_TYPE MATCHES Release)
  target_compile_options(libmuog PUBLIC -O3)
  target_compile_definitions(libmuog PUBLIC SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  target_compile_options(libmuog PUBLIC -Wall -Wextra -Wconversion -Wfloat-conversion -Wsign-conversion)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  target_compile_options(libmuog PUBLIC -Wall -Wextra -Wconversion -fp-model strict)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  target_compile_options(libmuog PUBLIC /W4)
endif()

# Set spdlog_pch.hpp as the precompiled header
if(USE_PCH)
  foreach(target nagamu_util nagamu_grid nagamu_geom nagamu_detector nagamu_path nagamu_data nagamu_map nagamu_recon nagamu_etc libmuog)
    target_precompile_headers(${target} PRIVATE ${PROJECT_SOURCE_DIR}/include/spdlog_pch.hpp)
    target_link_libraries(${target} PUBLIC spdlog::spdlog fmt::fmt)    
  endforeach()
endif()
