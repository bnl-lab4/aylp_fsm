#ifndef AYLP_TIMING_H_
#define AYLP_TIMING_H_

#include <stdbool.h>
#include <stddef.h>

/** Loop timing, published by the device that paces the loop and read by
 * devices that need to know how long the pipeline actually takes.
 *
 * Only the pacing device -- the source, which is the one that blocks waiting
 * for the next frame -- can tell WORK from WAITING. For any device further
 * down the pipeline the two are algebraically inseparable: the gap between its
 * own return and its next call is exactly (loop period - its own duration), an
 * identity that says nothing about how long anyone else took. The source, by
 * contrast, measures the whole downstream pipeline with no wait in it at all,
 * simply by timing from its own return to its next call.
 *
 * So the source publishes here, and controllers that need to know the transport
 * delay in frames (see fsp's "delay": "auto") read it.
 *
 * One process, one loop thread: publisher and readers are the same thread, so
 * this needs no locking. Devices loaded as plugins share the executable's copy
 * -- anyloop is linked with export_dynamic -- so an out-of-tree camera source
 * can publish to the same channel an in-tree controller reads.
 */
struct aylp_loop_timing {
	/** Seconds between successive frames delivered by the source. */
	double frame_period;
	/** Seconds the rest of the pipeline takes per iteration: centroiding,
	* control, DAC write and sinks, measured from the source's return to
	* its next call. This is the part of the loop delay that is WORK, and
	* it is the part that scales with what the pipeline is asked to do. */
	double compute_time;
	/** Iterations the two averages above were taken over. */
	size_t samples;
	/** Sensor integration time per frame, seconds; 0 when the source does
	* not know or does not expose one. Half of it is the lag from the
	* midpoint of a frame's integration -- where its information sits -- to
	* the end of that integration, which is the only part of the sensor's
	* contribution to the loop delay that software can know. */
	double exposure;
};


/** Publish a fresh timing window. Called by the pacing source device.
 * @param frame_period: mean seconds per frame over the window.
 * @param compute_time: mean seconds of downstream work over the window.
 * @param samples: iterations averaged.
 */
void aylp_timing_publish(double frame_period, double compute_time,
	size_t samples);


/** Publish the source's integration time, in seconds. Called once at init;
 * kept separate from the per-window publish because it does not change. */
void aylp_timing_publish_exposure(double exposure);


/** Read the most recently published window.
 * @param out: filled in on success; untouched otherwise.
 * @return true if a source has published at least once, false if not (no
 * source in the pipeline publishes timing, or the loop is too young).
 */
bool aylp_timing_get(struct aylp_loop_timing *out);


/** Frame geometry the source is delivering, in pixels.
 *
 * Same problem as the timing above: a centroider normalizes its output over the
 * whole frame, so everything after it works in units of "half a frame" and the
 * pixel scale is (dim - 1) / 2 -- but the frame's dimensions are gone from the
 * pipeline by then, since the header now describes a 2-element vector. A source
 * that sizes its own ROI (asi_source's auto ROI) makes this worse: the number
 * cannot even be written into the config ahead of time.
 */
struct aylp_frame_geometry {
	size_t height, width;	// px
};


/** Publish the frame size. Called once by the source, after its ROI is
 * resolved. */
void aylp_frame_geometry_publish(size_t height, size_t width);


/** Read the published frame size.
 * @param out: filled in on success; untouched otherwise.
 * @return true if a source has published a frame size, false otherwise.
 */
bool aylp_frame_geometry_get(struct aylp_frame_geometry *out);


#endif	// include guard
