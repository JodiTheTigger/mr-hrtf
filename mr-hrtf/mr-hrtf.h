// mr-hrtf. A simple 3d sound head related transfer function (HRTF) filter
// Copyright (C) 2015 Richard Maxwell <jodi.the.tigger@gmail.com>
// This file is part of mr-hrtf
// mr-hrtf is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Affero General Public License for more details.
// You should have received a copy of the GNU Affero General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>

#ifndef MR_HRTF_H
#define MR_HRTF_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C"
{
#endif

// -----------------------------------------------------------------------------
// Internal housekeeping
// -----------------------------------------------------------------------------
#ifndef MR_HRTF_REQUIRE
#include <assert.h>
#define MR_HRTF_REQUIRE(test) assert(test)
#endif

// -----------------------------------------------------------------------------
// Subject and Ear Descriptions (Anthropometry)
// -----------------------------------------------------------------------------

#define CIPIC_COUNT_X     17
#define CIPIC_COUNT_D     8
#define CIPIC_COUNT_THETA 2
#define CIPIC_COUNT_ALL   CIPIC_COUNT_X \
                          + (2 * (CIPIC_COUNT_D + CIPIC_COUNT_THETA))

typedef struct Mr_hrtf_anthropometry_raw
{
    float data[CIPIC_COUNT_ALL];
}
Mr_hrtf_anthropometry_raw;

typedef struct Mr_hrtf_anthropometry_pinnae_short
{
    float d[CIPIC_COUNT_D];
    float theta[CIPIC_COUNT_THETA];
}
Mr_hrtf_anthropometry_pinnae_short;

typedef struct Mr_hrtf_anthropometry_short
{
    float x[CIPIC_COUNT_X];

    Mr_hrtf_anthropometry_pinnae_short pinnae_left;
    Mr_hrtf_anthropometry_pinnae_short pinnae_right;
}
Mr_hrtf_anthropometry_short;


typedef struct Mr_hrtf_anthropometry_pinnae_long
{
    float cavum_concha_height;
    float cymba_concha_height;
    float cavum_concha_width;
    float fossa_height;
    float pinna_height;
    float pinna_width;
    float intertragal_incisure_width;
    float cavum_concha_depth;

    float pinna_rotation_angle;
    float pinna_flare_angle;
}
Mr_hrtf_anthropometry_pinnae_long;

typedef struct Mr_hrtf_anthropometry_long
{
    float head_width;
    float head_height;
    float head_depth;
    float pinna_offset_down;
    float pinna_offset_back;
    float neck_width;
    float neck_height;
    float neck_depth;
    float torso_top_width;
    float torso_top_height;
    float torso_top_depth;
    float shoulder_width;
    float head_offset_forward;
    float height;
    float seated_height;
    float head_circumference;
    float shoulder_circumference;

    Mr_hrtf_anthropometry_pinnae_long pinnae_left;
    Mr_hrtf_anthropometry_pinnae_long pinnae_right;
}
Mr_hrtf_anthropometry_long;

typedef union Mr_hrtf_anthropometry
{
    Mr_hrtf_anthropometry_raw   r;
    Mr_hrtf_anthropometry_short s;
    Mr_hrtf_anthropometry_long  l;
}
Mr_hrtf_anthropometry;

// -----------------------------------------------------------------------------
// HRTF Structures
// -----------------------------------------------------------------------------

typedef struct Mr_hrtf_location
{
    float azimuth;
    float elevation;
}
Mr_hrtf_location;

typedef struct Mr_hrtf_source
{
    Mr_hrtf_location location;
    short* source;
}
Mr_hrtf_source;

// -----------------------------------------------------------------------------
// Setup structures
// -----------------------------------------------------------------------------

enum Mr_hrtf_result
{
      MR_HRTF_RESULT_ALL_OK

    , MR_HRTF_RESULT_SETTINGS_NULL

    , MR_HRTF_RESULT_HRIM_CANT_OPEN
    , MR_HRTF_RESULT_HRIM_INVALID_SIZE
    , MR_HRTF_RESULT_HRIM_INVALID_VERSION

    , MR_HRTF_RESULT_HRIR_CANT_OPEN
    , MR_HRTF_RESULT_HRIR_INVALID_SIZE
    , MR_HRTF_RESULT_HRIR_INVALID_VERSION

    , MR_HRTF_RESULT_COUNT
};

// opaque handle
struct Mr_hrtf;
typedef struct Mr_hrtf Mr_hrtf;

struct Mr_hrtf_model_result
{
    Mr_hrtf_result  result;
    Mr_hrtf*        model;
};

struct Mr_hrtf_owned
{
    void* p;
    unsigned length_in_bytes;
};
//<! hrtf will free the array using free() when finished with it.

typedef void*         (*Mr_hrtf_malloc)(size_t);
typedef void          (*Mr_hrtf_free)(void*);
typedef Mr_hrtf_owned (*Mr_hrtf_ro_file_data)(Mr_hrtf_malloc, const char*);

typedef struct Mr_hrtf_settings
{
    unsigned               channel_count;
    // channel_count is the maximum amount of indivual sound streams to mix
    // into left and right outputs. It is assumed that when passing in the
    // sound streams that they stay on the same channel.

    unsigned               size_input;
    // how many input samples to process in one call to
    // mr_hrtf_mix_to_interleaved_lr or mr_hrtf_mix_to_seperate_lr.
    // If the amount of input samples are not enough toll fill in the input
    // buffer then it's assumed the buffer is padded to keep it the full
    // length.

    Mr_hrtf_anthropometry* anthropemetry;
    Mr_hrtf_anthropometry* anthropemetry_weights;
    // Set to nullptr if you just want to use the defaults.
    // Otherwise anthropemetry contains the data that you want the best
    // fit for (set unused values to 0.0f) and anthropemetry_weights if set
    // are how important each value is in choosing the result. 1.0f is normal,
    // 0.0f is ignored.

    Mr_hrtf_malloc       malloc;
    Mr_hrtf_free         free;
    Mr_hrtf_ro_file_data file_read;
    // Memory management and file io. If all are set to NULL, then standard
    // malloc, free, and an inbuilt file_read are used.
}
Mr_hrtf_settings;

// -----------------------------------------------------------------------------
// API (setup, cleanup)
// -----------------------------------------------------------------------------

Mr_hrtf_settings mr_hrtf_get_default_settings
(
      unsigned channel_count
    , unsigned size_input = 1024
);

Mr_hrtf_model_result mr_hrtf_get_model(const Mr_hrtf_settings* settings);
//<! Gets an opaque pointer to the hrtf engine setup for the passed in
//<! anthropemetry.

void mr_hrtf_free_model(Mr_hrtf* model);

// -----------------------------------------------------------------------------
// API (mixing)
// -----------------------------------------------------------------------------

void mr_hrtf_mix_to_interleaved_lr
(
      Mr_hrtf*               mix
    , const Mr_hrtf_source*  channels
    , unsigned               channel_count
    , short*                 destination
);

void mr_hrtf_mix_to_seperate_lr
(
      Mr_hrtf*              mix
    , const Mr_hrtf_source* channels
    , unsigned              channel_count
    , short*                left
    , short*                right
);

#ifdef __cplusplus
}
#endif

#endif // #define MR_HRTF_H
