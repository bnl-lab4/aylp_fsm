#ifndef AYLP_DEVICES_TEST_SOURCE_H_
#define AYLP_DEVICES_TEST_SOURCE_H_

#include "anyloop.h"

struct aylp_test_source_data {
	// one of ["vector", "matrix", "matrix_uchar"]
	aylp_type type;
	// one of ["constant", "sine", "noise"]
	unsigned kind;
	// size of vector or height of matrix
	size_t size1;
	// width of matrix, if applicable
	size_t size2;
	// frequency of sine oscillation, if applicable; units are radians per
	// proc() call
	double frequency;
	// optional MULTITONE for kind "sine": if nfreqs>0, the output is the sum
	// of nfreqs equal-amplitude sines at freqs[] (radians per proc() call),
	// each scaled by amplitude/nfreqs. Lets one in-loop run integrate every
	// tone over the whole record (max SNR per frequency) so a phase-vs-freq
	// slope gives the transport delay. "frequency" is ignored when nfreqs>0.
	size_t nfreqs;
	double *freqs;
	// amplitude of sine wave, if applicable
	double amplitude;
	// offset of sine wave or value of constant
	double offset;
	// PRNG seed for kind "noise" (fixed default so the white-noise stimulus
	// is reproducible; the emitted sequence is recorded anyway)
	unsigned seed;
	// to put in pipeline
	union {
		gsl_vector *vector;
		gsl_matrix *matrix;
		gsl_matrix_uchar *matrix_uchar;
	};
	// accumulator
	size_t acc;
};

// initialize test_source device
int test_source_init(struct aylp_device *self);

// process test_source device once per loop
int test_source_proc(struct aylp_device *self, struct aylp_state *state);

// close test_source device when loop exits
int test_source_fini(struct aylp_device *self);

#endif

