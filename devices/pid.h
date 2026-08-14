#ifndef AYLP_DEVICES_PID_H_
#define AYLP_DEVICES_PID_H_

#include <stdbool.h>
#include <time.h>
#include "anyloop.h"

// narrowband internal-model line rejection (vector [y,x] type only): one
// adaptive quadrature oscillator per disturbance line. The two weights
// integrate the error demodulated at the line frequency (filtered-x LMS),
// which is the numerically robust realization of a resonator with poles ON
// the unit circle: infinite loop gain exactly at f, so a stable sinusoidal
// disturbance is rejected even above the loop's crossover, at the cost of a
// notch in S only ~g*|P|/pi Hz wide instead of a broadband waterbed hit.
struct aylp_pid_line {
	double f;	// line frequency (Hz)
	double g;	// adaptation gain (1/s); notch width ~ g*|P·S|/pi Hz
	double phi;	// phase of the command->error path at f (rad); the
			// demodulator is rotated by this so the weight update
			// descends the true gradient. Defaults to the delay lag
			// -2pi*f*line_delay; for lines below ~1.5x crossover the
			// sensitivity function's phase advance matters too, so
			// give an explicit line_phase there.
	double a, b;	// quadrature weights (command units)
	double th;	// running oscillator phase (rad)
};
#define AYLP_PID_MAX_LINES 8

struct aylp_pid_data {
	// param: consecutive iterations carrying AYLP_NO_SIGNAL after which the
	// output is parked at park_value and the controller's state is reset.
	// 0 disables parking, which is the behaviour this device has always had.
	//
	// Holding the last command is the right answer to a brief dropout -- it
	// leaves the loop where the beam last was -- but it is the wrong answer
	// to a beam that is actually gone: the integrator has no error it can
	// believe, and whatever DC correction it was holding is now arbitrary.
	// Past some number of frames the honest move is to stop pretending and
	// go to a known output. Set it well above the length of an ordinary
	// dropout so it can only fire on a real loss; measured on this rig,
	// ordinary dropouts run 1-3 frames while real losses run thousands.
	size_t park_after;
	// param: the output every element is driven to while parked (default 0,
	// which is 0 V through a DAC configured with no offset)
	double park_value;
	// consecutive no-signal iterations seen so far, and whether the output
	// is currently parked (so the transitions log once, not every frame)
	size_t no_signal;
	bool parked;

	// param type in ["vector", "matrix"]
	aylp_type type;
	// units to output
	aylp_units units;
	// accumulated error
	union {
		gsl_vector *acc_v;
		gsl_matrix *acc_m;
	};
	// previous error
	union {
		gsl_vector *pre_v;
		gsl_matrix *pre_m;
	};
	// correction result
	union {
		gsl_vector *res_v;
		gsl_matrix *res_m;
	};
	// previous (filtered) derivative, for the derivative-filter recursion
	union {
		gsl_vector *dfd_v;
		gsl_matrix *dfd_m;
	};
	// previous timestamp
	struct timespec tp;
	// pid params; clamp bounds the total correction, and the accumulator is
	// bounded to clamp/|i| so the integral term can reach that limit but never
	// wind past it (see acc_limit in pid.c)
	double p, i, d, clamp, g;
	// param: seconds to hold the command at zero, with the integrator kept
	// empty, before the loop is allowed to act
	double start_delay;
	// CLOCK_MONOTONIC time of the first proc (s); 0 until seeded
	double t0;
	// true once the startup hold has been released
	bool started;
	double p_y, i_y, d_y, g_y;
	// derivative filter: low-pass the D term at this cutoff (Hz) to limit
	// noise amplification. <=0 disables it (pure derivative). Per-axis, with
	// the same convention as d: dfilt_y < 0 inherits the base dfilt.
	double dfilt, dfilt_y;
	// deadband: errors with |error| < deadband are treated as zero, to
	// suppress quantization/relay limit cycles around the null. Per-axis.
	double deadband, deadband_y;
	// optional series lead compensator on the PID output, per axis:
	//   C(s) = lead_g * (s + 2pi*lead_fz) / (s + 2pi*lead_fp)
	// DC gain = lead_g*fz/fp, HF gain = lead_g. The pole (lead_fp > 0) enables
	// the section; the zero may sit at 0 Hz (lead_fz = 0) for a band-limited
	// differentiator s/(s+wp). fp > fz gives phase lead. Base fields are the x
	// axis; _y fields are the y axis (independent, no inheritance, so x-only /
	// y-only / both-different are all expressible).
	double lead_fz, lead_fp, lead_g;
	double lead_fz_y, lead_fp_y, lead_g_y;
	// lead compensator state: previous raw input and filtered output (vector
	// type only). Bilinear coeffs are recomputed each step from the live dt.
	gsl_vector *lead_in_v;
	gsl_vector *lead_out_v;
	// narrowband line rejection, per axis: [0] = y (element 0), [1] = x
	// (element 1); see struct aylp_pid_line above
	struct aylp_pid_line line[2][AYLP_PID_MAX_LINES];
	size_t n_lines[2];
	// total command->error loop delay (s) used for the default demod phase
	double line_delay;
};

// initialize pid device
int pid_init(struct aylp_device *self);

// process pid device once per loop
int pid_proc(struct aylp_device *self, struct aylp_state *state);

// close pid device when loop exits
int pid_fini(struct aylp_device *self);

#endif


