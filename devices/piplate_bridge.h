#ifndef AYLP_DEVICES_PIPLATE_BRIDGE_H_
#define AYLP_DEVICES_PIPLATE_BRIDGE_H_

#include "anyloop.h"

struct aylp_piplate_bridge_data {
	int serial_fd;		// /dev/ttyACM0 (or configured port) file descriptor
	int board;		// DAQC board address (0–7)
	bool use_fg;		// use hardware function generator instead of setDAC
	// DAC mode -- one stage can drive several DAC channels in one proc call,
	// so each of these is an array of length n_outputs
	size_t n_outputs;	// number of DAC channels driven by this stage
	int *channels;		// DAC output channel per output; DAC: 0–3
	int *indices;		// which pipeline vector element to command
	double *scales;		// volts per pipeline unit
	double *offsets;	// volts at zero command
	int *last_codes;	// last DAC code written per channel (-1 = never)
	int *pending_codes;	// codes staged this iteration (-1 = not staged);
				// committed to last_codes only after a good write
	double *start_delays;	// per channel: seconds to hold at `offset` before
				// the pipeline command is allowed through
	bool *released;		// per channel: has the startup hold expired?
	bool has_start_delay;	// any channel has a nonzero start_delay
	double t0;		// CLOCK_MONOTONIC time of the first proc (s)
	bool wait_response;	// block for the BRIDGEplate ack after each write
	bool skip_unchanged;	// only write a channel when its voltage changes
	char *cmdbuf;		// scratch buffer for one iteration's commands
	size_t cmdbuf_sz;	// allocated size of cmdbuf
	// serial transmit model, so we never queue more command bytes than the
	// wire can carry: anything sitting in the tty buffer is transport delay
	long baud;		// configured line rate (bits/s)
	double byte_time;	// seconds per byte on the wire (8N1 => 10 bits)
	double cmd_time;	// seconds the board takes to service one setDAC:
				// >0 => tx_admit paces per command (measured at
				// init or set via cmd_time_us); <0 => pace on
				// byte_time; 0 => calibrate at init
	double max_backlog;	// seconds of queued tx we tolerate before dropping
	double tx_free;		// CLOCK_MONOTONIC time the tx queue drains (s)
	bool tx_primed;		// tx_free has been seeded
	// diagnostics, reported once per PIPLATE_DIAG_PERIOD
	long diag_writes;	// iterations that reached the wire
	long diag_drops;	// iterations dropped because the line was busy
	long diag_skips;	// channels skipped because the DAC code was unchanged
	long diag_cmds;		// setDAC commands that reached the wire
	double diag_t0;		// start of the current diagnostic window (s)
	// function generator mode
	int channel;		// FG channel: 1–2
	int fg_type;		// 1=sine 2=triangle 3=square 4=sawtooth 5=inv_saw 6=noise 7=sinc
	int fg_frequency;	// Hz, 10–10000
	int fg_level;		// 4=full 3=half 2=quarter 1=eighth
};

int piplate_bridge_init(struct aylp_device *self);
int piplate_bridge_proc(struct aylp_device *self, struct aylp_state *state);
int piplate_bridge_fini(struct aylp_device *self);

#endif
