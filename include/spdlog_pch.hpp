/// @file spdlog_pch.hpp
/// @brief Precompiled header for spdlog logging library
///
/// @details
/// This header provides centralized inclusion of spdlog library headers for
/// precompilation, reducing compilation time across the project.
///
/// Key features:
/// - Sets SPDLOG_ACTIVE_LEVEL to TRACE for fine-grained logging control
/// - Includes spdlog core headers (common.h, spdlog.h)
/// - Intended for PCH (Precompiled Header) usage via build system
///
/// @note Thread-safety: spdlog is thread-safe by default. See spdlog documentation
/// for details on logger registration and concurrent usage.
///
/// @note This file should be included via the build system's PCH mechanism,
/// not directly in source files, to maximize compilation benefits.
#pragma once

#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#include <spdlog/common.h> // level_enum, SPDLOG_FUNCTION
#include <spdlog/spdlog.h>
