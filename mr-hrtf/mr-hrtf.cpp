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

#include "mr-hrtf.h"

#include "kiss_fft.h"

#include <array>
#include <cstdint>
#include <cmath>

#include <stdio.h>

#define CIPIC_COUNT_HRIR_AZIMUTH    25
#define CIPIC_COUNT_HRIR_ELEVATION  50
#define CIPIC_COUNT_HRIR_SAMPLES    200
#define CIPIC_COUNT_SUBJECTS        45
#define CIPIC_COUNT_HRIR            \
            (CIPIC_COUNT_HRIR_AZIMUTH * CIPIC_COUNT_HRIR_ELEVATION)

static constexpr unsigned HRTF_LENGTH  = CIPIC_COUNT_HRIR_SAMPLES;
static constexpr auto     HRTF_OVERLAP = HRTF_LENGTH - 1;

static constexpr auto     VERSION_HRIM = 0u;
static constexpr auto     VERSION_HRIR = 0u;

static constexpr auto     INT16_FLOAT_MAX = ((1 << 15) - 1.0f);
static constexpr auto     INT16_FLOAT_MIN = (-1.0f * (1 << 15));

using Overlap_array = std::array<std::int16_t, HRTF_LENGTH>;
// Use length and waste 4 bytes so it's 16 byte aligned.

namespace {
// -----------------------------------------------------------------------------
// http://the-witness.net/news/2012/11/scopeexit-in-c11/

template <typename F>
struct Scoped_exit
{
    Scoped_exit(F f) : f(f) {}
    ~Scoped_exit() { f(); }
    F f;
};

template <typename F>
Scoped_exit<F> make_scoped_exit(F f)
{
    return Scoped_exit<F>(f);
};

#define MR_HRTF_STRING_JOIN2(arg1, arg2) MR_HRTF_DO_STRING_JOIN2(arg1, arg2)
#define MR_HRTF_DO_STRING_JOIN2(arg1, arg2) arg1 ## arg2
#define MR_HRTF_SCOPED_EXIT(code) \
    auto MR_HRTF_STRING_JOIN2(scope_exit_, __LINE__) = \
        make_scoped_exit([=](){code;})

// -----------------------------------------------------------------------------

template<typename T>
struct Span
{
    T*          p;
    std::size_t count;

    auto operator[](int i)       ->       T& { return p[i]; }
    auto operator[](int i) const -> const T& { return p[i]; }

    auto           data()            -> T* { return p; }
    auto constexpr data()      const -> T* { return p; }

    auto constexpr size()      const -> std::size_t { return count; }

    auto           begin()        -> T* { return p; }
    auto constexpr begin()  const -> T* { return p; }
    auto           end()          -> T* { return p + count; }
    auto constexpr end()    const -> T* { return p + count; }
};

template <typename T>
auto constexpr span(T* array, std::size_t array_size) -> Span<T>
{
    return
    {
          array
        , array_size
    };
}

template <typename T>
auto constexpr span(T& array) -> Span<typename T::value_type>
{
    return
    {
          array.data()
        , array.size()
    };
}

} // namespace

// -----------------------------------------------------------------------------

union Tag
{
    std::uint32_t dword;
    char          c[4];
};

struct Hrir_file
{
    using Hrir_pair =
        std::array<std::array<float, CIPIC_COUNT_HRIR_SAMPLES>, 2>;

    using Hrir_elevation = std::array<Hrir_pair, CIPIC_COUNT_HRIR_ELEVATION>;
    using Hrirs          = std::array<Hrir_elevation, CIPIC_COUNT_HRIR_AZIMUTH>;

    Tag                     tag;
    std::uint32_t           length;
    std::uint32_t           version;
    Mr_hrtf_anthropometry   anthropometry;
    Hrirs                   hrirs;
};

struct Hrim_map
{
    std::uint32_t           id;
    Mr_hrtf_anthropometry   anthropometry;
};

struct Hrim_file
{
    Tag                                        tag;
    std::uint32_t                              length;
    std::uint32_t                              version;
    std::array<Hrim_map, CIPIC_COUNT_SUBJECTS> anthropometry;
};

struct Mr_hrtf
{
    Mr_hrtf_settings settings;

    std::size_t  size_fft;

    kiss_fft_cfg fft;
    kiss_fft_cfg ifft;

    std::array<Mr_hrtf_location, CIPIC_COUNT_HRIR> locations;

    Span<kiss_fft_cpx> hrtfs;
    Span<kiss_fft_cpx> temps;
    Span<float>        mix;
    Span<int>          last_indicies;

    std::array<Span<Overlap_array>, 2> channel_overlaps;
};

enum Mr_hrtf_temp
{
      MR_HRTF_TEMP_INPUT_PRE_FFT
    , MR_HRTF_TEMP_INPUT_FFT
    , MR_HRTF_TEMP_OUTPUT_FFT
    , MR_HRTF_TEMP_OUTPUT_FFT_INVERSE
    , MR_HRTF_TEMP_OUTPUT_FFT_LAST
    , MR_HRTF_TEMP_OUTPUT_FFT_INVERSE_LAST

    , MR_HRTF_TEMP_COUNT
};

// -----------------------------------------------------------------------------
// Malloc and array math helpers
// -----------------------------------------------------------------------------

template<typename T>
auto malloc_span(unsigned count, Mr_hrtf_malloc mr_malloc) -> Span<T>
{
    auto bytes = sizeof(T) * count;
    return Span<T>{ static_cast<T*>(mr_malloc(bytes)), count};
}

template<typename T>
void free_span(Span<T>& span, Mr_hrtf_free mr_free)
{
    mr_free(span.data());
    span = {nullptr, 0};
}

enum Mr_hrtf_side
{
      MR_HRTF_SIDE_LEFT  = 0
    , MR_HRTF_SIDE_RIGHT = 1

    , MR_HRTF_SIDE_COUNT
};

auto get_hrtf
(
      Mr_hrtf& mix
    , Mr_hrtf_side side
    , unsigned index
)
-> kiss_fft_cpx*
{
    const auto size_fft = mix.size_fft;
          auto hrtfs    = mix.hrtfs;

    const auto hrtf_set_size = size_fft * CIPIC_COUNT_HRIR;

    index *= size_fft;
    index += hrtf_set_size * side;

    return &hrtfs[index];
}

auto get_temp(Mr_hrtf& mix, Mr_hrtf_temp index) -> kiss_fft_cpx*
{
    return &mix.temps[index * mix.size_fft];
}

auto get_mix(Mr_hrtf& mix, Mr_hrtf_side side) -> Span<float>
{
    return span
    (
          &mix.mix[side * mix.settings.size_input]
        , mix.settings.size_input
    );
}

auto default_ro_file_data
(
    Mr_hrtf_malloc mr_hrtf_malloc,
    const char*    filename
)
-> Mr_hrtf_owned
{
    Mr_hrtf_owned result =
    {
          0
        , 0
    };

    if ((!filename) || (!mr_hrtf_malloc))
    {
        return result;
    }

    auto* f = fopen(filename, "r");

    if (!f)
    {
        return result;
    }

    MR_HRTF_SCOPED_EXIT(fclose(f));

    if (!fseek(f, 0, SEEK_END))
    {
        auto bytes = ftell(f);
        rewind(f);

        auto data = mr_hrtf_malloc(bytes);
        fread(data, 1, bytes, f);

        result.p               = data;
        result.length_in_bytes = bytes;
    }

    return result;
}

inline auto operator-
(
    Mr_hrtf_anthropometry lhs
    , const Mr_hrtf_anthropometry& rhs
)
-> Mr_hrtf_anthropometry
{
    for (unsigned i = 0; i < CIPIC_COUNT_ALL; ++i)
    {
        lhs.r.data[i] -= rhs.r.data[i];
    }

    return lhs;
}

inline auto operator*
(
    Mr_hrtf_anthropometry lhs
    , const Mr_hrtf_anthropometry& rhs
)
-> Mr_hrtf_anthropometry
{
    for (unsigned i = 0; i < CIPIC_COUNT_ALL; ++i)
    {
        lhs.r.data[i] *= rhs.r.data[i];
    }

    return lhs;
}

inline auto is_zero(const Mr_hrtf_anthropometry& lhs) -> bool
{
    for (unsigned i = 0; i < CIPIC_COUNT_ALL; ++i)
    {
        if (lhs.r.data[i] != 0.0f)
        {
            return false;
        }
    }

    return true;
}

inline auto sum(const Mr_hrtf_anthropometry& lhs) -> float
{
    float result = 0.0f;

    for (unsigned i = 0; i < CIPIC_COUNT_ALL; ++i)
    {
        result += lhs.r.data[i];
    }

    return result;
}

auto closest_anthro
(
      const Hrim_file& anthro
    , const Mr_hrtf_anthropometry& wanted
    , const Mr_hrtf_anthropometry& weights

)
-> std::uint32_t
{
    std::uint32_t result = 0xFFFFFFFF;
    float base_error = 1000000000.0f;

    // values of zero are treated as invalid / missing
    for (unsigned i = 0; i < CIPIC_COUNT_ALL; ++i)
    {
        const auto& current = anthro.anthropometry[i].anthropometry;
        if (is_zero(current))
        {
            continue;
        }

        const auto error          = wanted - current;
        const auto error_weighted = error * weights;
        const auto error_squared  = error_weighted * error_weighted;
        const auto error_value    = sum(error_squared);

        if (error_value < base_error)
        {
            base_error = error_value;
            result = anthro.anthropometry[i].id;
        }
    }

    return result;
}

static const Mr_hrtf_anthropometry g_simple_mask =
{
    // x1, x3, x6, x12, d1, d3, d5, d6
    {
          1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f
        , 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f

        , 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f
        , 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f
    }
};

int get_closest_index
(
      Mr_hrtf_location location
    , Span<Mr_hrtf_location> locations
)
{
    auto error = 10000000.0f;
    int result = -1;

    const auto size = locations.size();
    for (unsigned i = 0; i < size; ++i)
    {
        const auto& test = locations[i];

        auto diff =
            (
                  (location.azimuth - test.azimuth)
                * (location.azimuth - test.azimuth)
            )
            +
            (
                  (location.elevation - test.elevation)
                * (location.elevation - test.elevation)
            );

        if (diff < error)
        {
            result = i;
            error  = diff;
        }
    }

    return result;
}

// -----------------------------------------------------------------------------
// API (mixing)
// -----------------------------------------------------------------------------

void mr_hrtf_mix_to_seperate_lr_work
(
      Mr_hrtf*              mix
    , const Mr_hrtf_source* sources
    , unsigned              channels
)
{
    MR_HRTF_REQUIRE(mix     != nullptr);
    MR_HRTF_REQUIRE(sources != nullptr);
    MR_HRTF_REQUIRE(channels <= mix->settings.channel_count);
    MR_HRTF_REQUIRE(mix->settings.size_input >= HRTF_LENGTH);

    MR_HRTF_REQUIRE
    (
        mix->size_fft >= (mix->settings.size_input + HRTF_LENGTH - 1)
    );

    const auto size_input = mix->settings.size_input;
    const auto size_fft   = mix->size_fft;

    auto input_pre_fft      = get_temp(*mix, MR_HRTF_TEMP_INPUT_PRE_FFT);
    auto input_fft          = get_temp(*mix, MR_HRTF_TEMP_INPUT_FFT);
    auto output_fft         = get_temp(*mix, MR_HRTF_TEMP_OUTPUT_FFT);
    auto output_fft_last    = get_temp(*mix, MR_HRTF_TEMP_OUTPUT_FFT_LAST);
    auto output_fft_inverse = get_temp(*mix, MR_HRTF_TEMP_OUTPUT_FFT_INVERSE);

    auto output_fft_inverse_last =
        get_temp(*mix, MR_HRTF_TEMP_OUTPUT_FFT_INVERSE_LAST);

    memset(mix->mix.data(), 0, 4 * mix->mix.size());
    auto out = std::array<Span<float>, 2>
    {
          get_mix(*mix, MR_HRTF_SIDE_LEFT)
        , get_mix(*mix, MR_HRTF_SIDE_RIGHT)
    };

    const auto fft_mult = []
    (
          const kiss_fft_cpx* lhs
        , const kiss_fft_cpx* rhs
        , kiss_fft_cpx* out
        , unsigned size_fft
    )
    {
        for (unsigned i = 0; i < size_fft; ++i)
        {
            out[i].r = (lhs[i].r * rhs[i].r) - (lhs[i].i * rhs[i].i);
            out[i].i = (lhs[i].r * rhs[i].i) + (lhs[i].i * rhs[i].r);
        }
    };

    for (unsigned channel_index = 0; channel_index < channels; ++channel_index)
    {
        const auto& source = sources[channel_index];

        const auto hrtf_index =
            get_closest_index(source.location, span(mix->locations));

        const auto last_index = mix->last_indicies[channel_index];

        mix->last_indicies[channel_index] = hrtf_index;

        const auto* in = source.source;

        if (!in)
        {
            continue;
        }

        for (int side = 0; side < 2; ++side)
        {
            const auto hrtf =
                get_hrtf(*mix, static_cast<Mr_hrtf_side>(side), hrtf_index);

            const auto hrtf_last = [&]() -> kiss_fft_cpx*
            {
                if ((hrtf_index != last_index) && (last_index >=0))
                {
                    return get_hrtf
                    (
                          *mix
                        , static_cast<Mr_hrtf_side>(side)
                        , last_index
                    );
                }

                return nullptr;
            }();

            auto& overlap = mix->channel_overlaps[side][channel_index];

            for (unsigned i = 0; i < HRTF_OVERLAP; ++i)
            {
                input_pre_fft[i].r = overlap[i];
                input_pre_fft[i].i = 0;
            }
            for (unsigned i = 0 ; i < size_input; ++i)
            {
                input_pre_fft[i + HRTF_OVERLAP].r = in[i];
                input_pre_fft[i + HRTF_OVERLAP].i = 0;
            }
            for (unsigned i = size_input + HRTF_OVERLAP; i < size_fft; ++i)
            {
                input_pre_fft[i].r = 0;
                input_pre_fft[i].i = 0;
            }

            const auto offset_input_overlap = size_input - HRTF_OVERLAP;
            for (unsigned i = 0; i < HRTF_OVERLAP; ++i)
            {
                overlap[i] = in[i + offset_input_overlap];
            }

            kiss_fft(mix->fft, input_pre_fft, input_fft);
            fft_mult(input_fft, hrtf, output_fft, size_fft);
            kiss_fft(mix->ifft, output_fft,  output_fft_inverse);

            // Mix with the hrtf of the previous location if needed
            if (!hrtf_last)
            {
                auto& out_side = out[side];
                for (unsigned i = 0; i < size_input; ++i)
                {
                    out_side[i] +=
                        output_fft_inverse[i + HRTF_OVERLAP].r / size_fft;
                }
            }
            else
            {
                fft_mult(input_fft, hrtf_last, output_fft_last, size_fft);
                kiss_fft(mix->ifft, output_fft_last, output_fft_inverse_last);

                auto& out_side = out[side];
                const auto step = 1.0f / size_input;
                for (unsigned i = 0; i < size_input; ++i)
                {
                    const auto lerp   = (i * step);
                    const auto lerp_i = 1.0f - lerp;
                    const auto mixed  =
                    (
                          (lerp   * output_fft_inverse[i + HRTF_OVERLAP].r)
                        + (lerp_i * output_fft_inverse_last[i + HRTF_OVERLAP].r)
                    )
                    / size_fft;

                    out_side[i] += mixed;
                }
            }
        }
    }
}

auto inline constexpr clip_to_short(float to_clip) -> std::int16_t
{
    return std::round
    (
        (to_clip < INT16_FLOAT_MAX)
            ? (to_clip > INT16_FLOAT_MIN)
                ? to_clip
                : INT16_FLOAT_MIN
            : INT16_FLOAT_MAX
    );
}

void mr_hrtf_mix_to_interleaved_lr
(
      Mr_hrtf*               mix
    , const Mr_hrtf_source*  channels
    , unsigned               channels_count
    , short*                 destination
)
{
    mr_hrtf_mix_to_seperate_lr_work(mix, channels, channels_count);

    auto out_left  = get_mix(*mix, MR_HRTF_SIDE_LEFT);
    auto out_right = get_mix(*mix, MR_HRTF_SIDE_RIGHT);

    const auto size_input = mix->settings.size_input;
    for (unsigned i = 0; i < size_input; ++i)
    {
        destination[i * 2 + 0] = clip_to_short(out_left[i]);
        destination[i * 2 + 1] = clip_to_short(out_right[i]);
    }
}

void mr_hrtf_mix_to_seperate_lr
(
      Mr_hrtf*              mix
    , const Mr_hrtf_source* channels
    , unsigned              channels_count
    , short*                left
    , short*                right
)
{
    MR_HRTF_REQUIRE(left     != nullptr);
    MR_HRTF_REQUIRE(right    != nullptr);

    mr_hrtf_mix_to_seperate_lr_work(mix, channels, channels_count);

    auto out_left  = get_mix(*mix, MR_HRTF_SIDE_LEFT);
    auto out_right = get_mix(*mix, MR_HRTF_SIDE_RIGHT);

    const auto size_input = mix->settings.size_input;
    for (unsigned i = 0; i < size_input; ++i)
    {
        left[i]  = clip_to_short(out_left[i]);
        right[i] = clip_to_short(out_right[i]);
    }
}

auto mr_hrtf_get_default_settings
(
      unsigned channel_count
    , unsigned size_input
)
-> Mr_hrtf_settings
{
    return
    {
          channel_count
        , size_input
        , nullptr
        , nullptr

        , nullptr
        , nullptr
        , nullptr
    };
}

auto mr_hrtf_get_model(const Mr_hrtf_settings* settings) -> Mr_hrtf_model_result
{
    if (!settings)
    {
        return {MR_HRTF_RESULT_SETTINGS_NULL, nullptr};
    }

    Mr_hrtf_settings config = *settings;

    if
    (
           (!config.malloc)
        && (!config.free)
        && (!config.file_read)
    )
    {
        config.malloc    = malloc;
        config.free      = free;
        config.file_read = default_ro_file_data;
    }

    // -------------------------------------------------------------------------

    auto raw_hrim = config.file_read(config.malloc, "cipic.hrim");
    if (!raw_hrim.p)
    {
        return {MR_HRTF_RESULT_HRIM_CANT_OPEN, nullptr};
    }
    MR_HRTF_SCOPED_EXIT(config.free(raw_hrim.p));

    if (raw_hrim.length_in_bytes != sizeof(Hrim_file))
    {
        return {MR_HRTF_RESULT_HRIM_INVALID_SIZE, nullptr};
    }

    const auto* hrim_file = static_cast<Hrim_file*>(raw_hrim.p);
    if (hrim_file->version != VERSION_HRIM)
    {
        return {MR_HRTF_RESULT_HRIM_INVALID_VERSION, nullptr};
    }

    // -------------------------------------------------------------------------

    std::uint32_t hrir_id = 3;

    if (settings->anthropemetry)
    {
        const auto weights =
            settings->anthropemetry_weights
            ? settings->anthropemetry_weights
            : &g_simple_mask;

        hrir_id = closest_anthro
        (
              *hrim_file
            , *settings->anthropemetry
            , *weights
        );
    }

    char file_name[32];

    sprintf(file_name, "subject_%03d.hrir", hrir_id);

    const auto raw_hrir = config.file_read(config.malloc, file_name);

    if (!raw_hrir.p)
    {
        return {MR_HRTF_RESULT_HRIR_CANT_OPEN, nullptr};
    }
    MR_HRTF_SCOPED_EXIT(config.free(raw_hrir.p));

    if (raw_hrir.length_in_bytes != sizeof(Hrir_file))
    {
        return {MR_HRTF_RESULT_HRIR_INVALID_SIZE, nullptr};
    }

    const auto hrir_file = static_cast<Hrir_file*>(raw_hrir.p);
    if (hrir_file->version != VERSION_HRIR)
    {
        return {MR_HRTF_RESULT_HRIR_INVALID_VERSION, nullptr};
    }

    // -------------------------------------------------------------------------

    const auto size_input = settings->size_input;
    const auto size_fft   =
        ((size_input + HRTF_OVERLAP) + 63) & ~63;
    // Round up the fft size to the nearest 64 bins so that it's more fft
    // than dft. That is, we want more 4, 2, 3 and 5 radix butterflys than
    // large prime number butterflys.

    auto result = static_cast<Mr_hrtf*>(config.malloc(sizeof(Mr_hrtf)));

    result->settings = config;
    result->size_fft = size_fft;

    result->channel_overlaps[MR_HRTF_SIDE_LEFT]  =
        malloc_span<Overlap_array>(config.channel_count, config.malloc);

    result->channel_overlaps[MR_HRTF_SIDE_RIGHT] =
        malloc_span<Overlap_array>(config.channel_count, config.malloc);

    {
        size_t fft_size_in_bytes = 0;

        kiss_fft_alloc(size_fft, 0, nullptr, &fft_size_in_bytes);

        result->fft  =
            kiss_fft_alloc
            (
                  size_fft
                , 0
                , config.malloc(fft_size_in_bytes)
                , &fft_size_in_bytes
            );

        result->ifft =
            kiss_fft_alloc
            (
                  size_fft
                , 1
                , config.malloc(fft_size_in_bytes)
                , &fft_size_in_bytes
            );
    }

    // -------------------------------------------------------------------------

    result->hrtfs =
        malloc_span<kiss_fft_cpx>
        (
              size_fft
            * CIPIC_COUNT_HRIR
            * MR_HRTF_SIDE_COUNT
            , config.malloc
        );

    result->temps =
        malloc_span<kiss_fft_cpx>
        (
              size_fft
            * MR_HRTF_TEMP_COUNT
            , config.malloc
        );

    result->mix =
        malloc_span<float>
        (
              size_input
            * MR_HRTF_SIDE_COUNT
            , config.malloc
        );

    result->last_indicies =
        malloc_span<int>
        (
              config.channel_count
            , config.malloc
        );

    for (auto& index: result->last_indicies)
    {
        index = -1;
    }

    // -------------------------------------------------------------------------

    const auto az_lookup = std::array<float, CIPIC_COUNT_HRIR_AZIMUTH>
    {
          -80, -65, -55, -45, -40, -35, -30, -25, -20, -15, -10, -5
        , 0
        , 5,    10,  15,  20,  25,  30,  35,  40,  45,  55,  65,  80
    };

    unsigned loc_index = 0;
    for (unsigned az = 0 ; az < CIPIC_COUNT_HRIR_AZIMUTH; ++az)
    {
        for (unsigned e = 0 ; e < CIPIC_COUNT_HRIR_ELEVATION; ++e)
        {
            result->locations[loc_index].elevation = e * 5.625f + -45.0f;
            result->locations[loc_index].azimuth   = az_lookup[az];

            ++loc_index;
        }
    }

    // -------------------------------------------------------------------------

    auto hrir_pre_fft = get_temp(*result, MR_HRTF_TEMP_INPUT_PRE_FFT);
    memset(hrir_pre_fft, 0, size_fft * sizeof(kiss_fft_cpx));

    static_assert(MR_HRTF_SIDE_LEFT  == 0, "MR_HRTF_SIDE_LEFT isn't 0.");
    static_assert(MR_HRTF_SIDE_RIGHT == 1, "MR_HRTF_SIDE_RIGHT isn't 1.");

    unsigned index = 0;
    for (unsigned az = 0 ; az < CIPIC_COUNT_HRIR_AZIMUTH; ++az)
    {
        for (unsigned e = 0 ; e < CIPIC_COUNT_HRIR_ELEVATION; ++e)
        {
            for (unsigned side = 0; side < 2; ++side)
            {
                const auto& hrir = hrir_file->hrirs[az][e][side];

                auto hrtf =
                    get_hrtf(*result, static_cast<Mr_hrtf_side>(side), index);

                for (unsigned i = 0; i < CIPIC_COUNT_HRIR_SAMPLES; ++i)
                {
                    hrir_pre_fft[i].r = hrir[i];
                }

                kiss_fft(result->fft, hrir_pre_fft, hrtf);
            }

            ++index;
        }
    }

    return {MR_HRTF_RESULT_ALL_OK, result};
}

void mr_hrtf_free_model(Mr_hrtf *model)
{
    if ((!model) || (!model->settings.free))
    {
        return;
    }

    const auto& settings = model->settings;

    free_span(model->mix,           settings.free);
    free_span(model->temps,         settings.free);
    free_span(model->hrtfs,         settings.free);
    free_span(model->last_indicies, settings.free);

    free_span(model->channel_overlaps[MR_HRTF_SIDE_LEFT],  settings.free);
    free_span(model->channel_overlaps[MR_HRTF_SIDE_RIGHT], settings.free);

    settings.free(model->fft);
    settings.free(model->ifft);

    settings.free(model);
}
