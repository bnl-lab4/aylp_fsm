#include "timing.h"

// The channel itself is private: only the two accessors below are exported, so
// plugins reach it through the executable's PLT rather than through a copy
// relocation of the object.
static struct aylp_loop_timing timing;
static bool timing_valid;


void aylp_timing_publish(double frame_period, double compute_time,
	size_t samples)
{
	if (!(frame_period > 0.0) || !(compute_time >= 0.0) || !samples) return;
	timing.frame_period = frame_period;
	timing.compute_time = compute_time;
	timing.samples = samples;
	timing_valid = true;
}


void aylp_timing_publish_exposure(double exposure)
{
	if (exposure >= 0.0) timing.exposure = exposure;
}


bool aylp_timing_get(struct aylp_loop_timing *out)
{
	if (!timing_valid) return false;
	*out = timing;
	return true;
}


static struct aylp_frame_geometry geometry;
static bool geometry_valid;


void aylp_frame_geometry_publish(size_t height, size_t width)
{
	if (!height || !width) return;
	geometry.height = height;
	geometry.width = width;
	geometry_valid = true;
}


bool aylp_frame_geometry_get(struct aylp_frame_geometry *out)
{
	if (!geometry_valid) return false;
	*out = geometry;
	return true;
}
