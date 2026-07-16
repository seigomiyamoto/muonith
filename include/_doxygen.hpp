/// @file _doxygen.hpp
/// @brief Doxygen documentation structure and group definitions
/// @details Defines the main page and all documentation groups for organizing the codebase documentation.
/*! \mainpage
 * \brief MUONITH Project Documentation
 *
 * \section intro Introduction
 * This project provides a comprehensive ray tracing framework for muon radiography
 * and related computational physics simulations. The codebase is organized into
 * several major functional groups described below.
 *
 * \section groups Top-level Groups
 * - \ref RayTracingCore - Grid classes and spatial data structures
 * - \ref ray_tracing_functions - Ray traversal and intersection algorithms
 * - \ref Utils - General-purpose utilities and helper functions
 * - \ref Domain - Physical domain models (detectors, terrain)
 * - \ref Algorithm - Computational algorithms (matrices, flux calculations)
 * - \ref Parameters - Runtime configuration and parameter handling
 * - \ref Execution - Workflow orchestration and execution modules
 * - \ref Development - Debug macros and development tools
 *
 * \section subgroups Specialized Subgroups
 * - \ref basicTools - Basic utility functions
 * - \ref basicGridClasses - 1D/2D/3D grid classes
 * - \ref derivedGridClasses - Grid class extensions
 * - \ref geometryClasses - Geometric primitives (AABB, Ray)
 * - \ref detectorClasses - Detector geometry models
 * - \ref terrainClasses - Terrain representation
 * - \ref pathCalculation - Path length calculation for ray-terrain intersection
 * - \ref matrixClasses - Matrix operations and linear algebra
 * - \ref parameterClasses - RunCard parameter classes
 * - \ref idContainerClasses - Unique ID management
 * - \ref ExecModule - Simulation execution modules
 * - \ref macro_debug - Debugging and diagnostics macros
 * - \ref flux_tools - Muon flux tables and range calculations
 */

//==============================================================================
// Top-level group definitions
//==============================================================================

/*!
 * \defgroup RayTracingCore Ray Tracing Core Components
 * Grid classes, geometric primitives (AABB, Ray), and spatial data structures for ray tracing
 */

/*!
 * \defgroup ray_tracing_functions Ray Tracing Functions
 * Algorithms and utility routines for ray traversal and intersection tests
 */

/*!
 * \defgroup Utils Utilities
 * Miscellaneous utility namespaces and classes for various purposes
 */

/*!
 * \defgroup Domain Domain Model
 * Physical domain representations including detector geometry and terrain models
 */

/*!
 * \defgroup Algorithm Algorithm Implementation
 * Computational algorithms including matrix operations and muon flux calculations
 */

/*!
 * \defgroup Parameters Parameters and Configuration
 * Runtime configuration parameters loaded from RunCard and parameter handling utilities
 */

/*!
 * \defgroup Execution Execution Control
 * Execution modules and workflow orchestration for simulation tasks
 */

/*!
 * \defgroup Development Development Tools
 * Debug macros and development utilities for diagnostics and testing
 */

//==============================================================================
// Subgroup definitions (existing groups with hierarchy)
//==============================================================================

/*!
 * \defgroup basicTools Basic utility functions used throughout the codebase
 * \ingroup Utils
 */

/*!
 * \defgroup basicGridClasses Basic grid classes (Grid1d, Grid2d, Grid3d)
 * \ingroup RayTracingCore
 * Foundation is 1D Grid1d. \n
 * Grid2d has private members x_axis and y_axis of Grid1d. \n
 * Grid3d has private members x_axis, y_axis, and z_axis of Grid1d. \n
 */

/*!
 * \defgroup derivedGridClasses Classes derived from Grid1d, Grid2d, Grid3d
 * \ingroup basicGridClasses
 * \note This group has been merged into basicGridClasses. Maintained for backward compatibility.
 * Classes inheriting from Grid1d, Grid2d, Grid3d.
 */

/*!
 * \defgroup detectorClasses Detector-related classes
 * \ingroup Domain
 */

/*!
 * \defgroup terrainClasses Terrain-related classes
 * \ingroup Domain
 */

/*!
 * \defgroup matrixClasses Matrix-related classes and namespaces
 * \ingroup Algorithm
 */

/*!
 * \defgroup parameterClasses Parameter-related classes
 * \ingroup Parameters
 * Classes that receive and reflect information from RunCard
 */

/*!
 * \defgroup geometryClasses Geometry and collision detection classes (AABB, Ray)
 * \ingroup RayTracingCore
 */

/*!
 * \defgroup idContainerClasses Classes for assigning and storing unique IDs
 * \ingroup basicTools
 * \note This group has been merged into basicTools. Maintained for backward compatibility.
 */

/*!
 * \defgroup pathCalculation Path length calculation functions for terrain. Some support per-voxel operations.
 * \ingroup Domain
 */

/*!
 * \defgroup ExecModule Execution modules
 * \ingroup Execution
 */

/*!
 * \defgroup macro_debug Debug macros
 * \ingroup Development
 */

/*!
 * \defgroup flux_tools Flux tables for muon range calculations
 * \ingroup Algorithm
 */
