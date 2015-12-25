# mr-hrtf. A simple 3d sound head related transfer function (HRTF) filter
# Copyright (C) 2015 Richard Maxwell <jodi.the.tigger@gmail.com>
# This file is part of mr-hrtf
# mr-hrtf is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

# ----------------------------------------------------------------------------
# Project info
# ----------------------------------------------------------------------------

set(PROJECT_NAME "mr-hrtf")
project(${PROJECT_NAME} CXX C)

cmake_minimum_required(VERSION 3.3.0 FATAL_ERROR)

# ----------------------------------------------------------------------------
# Prerequisites
# ----------------------------------------------------------------------------

include_directories(${CMAKE_SOURCE_DIR})

# C++11/14 features I use.
set(NEEDED_FEATURES
    cxx_auto_type
    cxx_trailing_return_types
    cxx_right_angle_brackets
    cxx_nullptr)

# ----------------------------------------------------------------------------
# Source
# ----------------------------------------------------------------------------
set(SOURCE
    mr-hrtf/mr-hrtf.cpp
    mr-hrtf/mr-hrtf.h
    mr-hrtf/_kiss_fft_guts.h
    mr-hrtf/kiss_fft.c
    mr-hrtf/kiss_fft.h)

# ----------------------------------------------------------------------------
# Link
# ----------------------------------------------------------------------------

add_library(${PROJECT_NAME} ${SOURCE})
target_compile_features(${PROJECT_NAME} PRIVATE ${NEEDED_FEATURES})
target_link_libraries(${PROJECT_NAME} kiss_fft)
