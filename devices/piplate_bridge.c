// PiPlate BRIDGEplate device -- drives one DAQC2 DAC channel via USB serial.
//
// DAC mode (default): maps a pipeline vector element from MINMAX [-1, 1] to a
// voltage using offset + cmd*scale, and calls DAQC2.setDAC each iteration.
//
// A single DAC-mode stage can drive several DAC channels at once: pass arrays
// for channel/index/scale/offset (all the same length, or a scalar to apply
// the same value to every channel). Example (four channels):
//   "channel": [0, 2, 1, 3],
//   "index":   [0, 0, 0, 1],
//   "scale":   [-2, -2, 0, 0],
//   "offset":  [2, 4, 3.09, 1.924]
//
// Each setDAC on the BRIDGEplate normally costs a full serial round-trip
// (~0.5-0.7 ms) because we block for the board's ack, so writing N channels
// every iteration costs ~N round-trips. The fix that keeps that off the loop:
//   * wait_response=false (fire-and-forget): don't block for the ack. All the
//     channels' commands are written in one buffer and we return immediately,
//     flushing the board's stale acks from the input buffer next iteration.
//     This takes the serial round-trip off the loop entirely, so refreshing
//     every channel every iteration stays cheap. Leave it true (the default)
//     when you want the ack for flow control / debugging.
//
// THE BOARD IS THE BOTTLENECK. The BRIDGEplate is a USB CDC-ACM device, so the
// configured `baud` does not pace the USB wire: the bytes arrive at USB speed
// and the real cost is the firmware parsing the ASCII command and running the
// SPI transaction to the DAQC2 -- ~0.5-0.7 ms per setDAC as seen in the ack
// round-trip. A 2 kHz control loop wants to send one every 0.5 ms. The excess
// cannot vanish: without backpressure it piles up in the tty/firmware buffers,
// where it becomes pure transport delay in the feedback path -- and transport
// delay destroys phase margin regardless of the PID gains. So in fire-and-forget
// mode a transmit model (see tx_admit) refuses a write while the board is still
// busy, dropping that iteration's command instead of queueing it. The DAC holds
// its previous value one iteration longer and the next iteration sends a fresher
// target. The model paces on the board's per-command service time, measured at
// init by timing acked setDACs (override with cmd_time_us; a negative value
// falls back to the old baud-derived byte-time pacing). Watch the
// "piplate_bridge: ... writes/s | ...% dropped" line: a high drop fraction means
// the loop is running faster than the actuator can possibly be commanded, and
// you should slow the loop (longer camera exposure) or speed up the board's
// firmware, not push the gains.
//
// skip_unchanged (default false) additionally suppresses the setDAC for any
// channel whose DAC code has not changed since the last write -- e.g. a constant
// scale=0 bias channel gets written once and then skipped. The comparison is on
// the 12-bit code, not the floating-point target: the DAQC2 resolves 1 mV, so a
// slow integrator ramp that moves the target by microvolts per iteration would
// otherwise write every iteration while never moving the output. Only turn this
// on if the DAC outputs hold their value on their own: a channel driving a
// low-impedance load (say 10k) will droop between refreshes, so it must be
// re-sent every iteration and skip_unchanged must stay off for it. When in
// doubt, leave it off and refresh everything every iteration.
//
// FG mode (mode="fg"): configures the DAQC2 hardware function generator at
// init (fgTYPE/fgFREQ/fgLEVEL/fgON) and passes through each iteration with no
// USB traffic. The hardware updates its DAC at 200 kHz autonomously.
//
// Params (DAC mode):
//   port      -- serial device path (default "/dev/ttyACM0")
//   board     -- DAQC board address 0-7 (default 0)
//   channel   -- DAC output channel(s) 0-3, scalar or array (default 0)
//   index     -- pipeline vector element(s) to command, scalar or array (default 0)
//   scale     -- volts per pipeline unit, scalar or array (default 1.0)
//   offset    -- volts at zero command, scalar or array (default 0.0)
//   start_delay -- seconds to hold a channel at its `offset` before letting the
//                pipeline command through, scalar or array (default 0). Measured
//                from the first proc, which is when every channel first reaches
//                its bias. Use it to give a slow coarse stage time to walk the
//                beam near centre before the fine channels are allowed to move:
//                a fine channel that chases the large startup error meanwhile
//                just drives itself into its rail and has no travel left when
//                the coarse one arrives. This parks the *output* only -- give
//                the upstream `pid` a matching start_delay too, otherwise its
//                integrator winds up during the hold and the fine channel is at
//                its rail the instant the hold expires.
//   wait_response -- block for the BRIDGEplate ack each write (default true;
//                    set false for fire-and-forget in a fast loop)
//   skip_unchanged -- only send a channel when its DAC code changes (default
//                    false; keep false for loaded outputs that droop)
//   baud      -- nominal line rate: 9600/19200/38400/57600/115200/230400/
//                460800/921600 (default 115200). A USB CDC-ACM bridge ignores
//                it on the wire; it only sets the pacing when cmd_time_us < 0
//                (legacy byte-time model).
//   cmd_time_us -- the board's per-setDAC service time in microseconds, used
//                by the transmit model. 0 (default): measure it at init by
//                timing acked setDACs (each channel is parked at its `offset`
//                during the ~30 ms this takes). >0: use the given value.
//                <0: pace on baud-derived byte time instead (legacy).
//   max_backlog_ms -- how many ms of command time may sit queued at the board
//                before an iteration is dropped (default 0, i.e. only write
//                when the board is idle: lowest latency). Raise it to trade
//                command latency for a higher write rate. Ignored when
//                wait_response is true, which self-limits.
//
// Params (FG mode, mode="fg"):
//   port      -- serial device path (default "/dev/ttyACM0")
//   board     -- DAQC board address 0-7 (default 0)
//   channel   -- FG channel 1 or 2 (default 1)
//   frequency -- Hz, 10–10000 (default 1000)
//   type      -- waveform: 1=sine 2=triangle 3=square 4=sawtooth
//                          5=inv_sawtooth 6=noise 7=sinc (default 1)
//   level     -- amplitude: 4=full 3=half 2=quarter 1=eighth (default 4)

#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <termios.h>
#include <time.h>
#include <unistd.h>
#include <json-c/json.h>

#include "anyloop.h"
#include "logging.h"
#include "piplate_bridge.h"
#include "xalloc.h"

// DAQC2 analog outputs are 12-bit over 0-4.095 V: exactly 1 mV per code. There
// is no point expressing a target any finer than that.
#define DAC_LSB      0.001
#define DAC_VMAX     4.095
#define DAC_CODE_MAX 4095
// seconds between piplate_bridge throughput reports
#define PIPLATE_DIAG_PERIOD 5.0
// a setDAC command is at most 24 ASCII bytes ("DAQC2.setDAC(7,3,4.095)\n")
#define PIPLATE_CMD_BYTES 24
// acked round-trips used to measure the board's per-command service time
#define PIPLATE_CALIB_WARMUP  4
#define PIPLATE_CALIB_SAMPLES 24

static speed_t baud_to_speed(long b)
{
	switch (b) {
	case 9600:   return B9600;
	case 19200:  return B19200;
	case 38400:  return B38400;
	case 57600:  return B57600;
	case 115200: return B115200;
	case 230400: return B230400;
#ifdef B460800
	case 460800: return B460800;
#endif
#ifdef B921600
	case 921600: return B921600;
#endif
	default:     return 0;
	}
}

static int serial_open(const char *port, speed_t speed)
{
	int fd = open(port, O_RDWR | O_NOCTTY);
	if (fd < 0) return -1;

	struct termios tty;
	if (tcgetattr(fd, &tty)) { close(fd); return -1; }
	cfmakeraw(&tty);
	cfsetispeed(&tty, speed);
	cfsetospeed(&tty, speed);
	tty.c_cc[VMIN]  = 0;
	tty.c_cc[VTIME] = 5;  // 500ms read timeout in 0.1s units
	if (tcsetattr(fd, TCSANOW, &tty)) { close(fd); return -1; }
	return fd;
}

static double monotonic_s(void)
{
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return t.tv_sec + 1e-9 * t.tv_nsec;
}

/** Round a target voltage to the DAC code that will actually be produced.
* Printing more precision than the hardware can resolve is worse than useless:
* it costs a byte of wire time per command, and it makes skip_unchanged fire on
* changes the DAC cannot express, so a slow integrator ramp writes every single
* iteration while never moving the output. */
static int volts_to_code(double volts)
{
	if (volts < 0.0) volts = 0.0;
	if (volts > DAC_VMAX) volts = DAC_VMAX;
	long c = lround(volts / DAC_LSB);
	if (c < 0) c = 0;
	if (c > DAC_CODE_MAX) c = DAC_CODE_MAX;
	return (int)c;
}

/** Ask the transmit model whether `ncmds` more commands (`len` bytes) may go
* out now.
*
* The tty buffer will accept far more than the board can service -- one setDAC
* costs the firmware ~0.5-0.7 ms of parse + SPI time, while a 2 kHz loop offers
* a new command every 0.5 ms. Commands queued behind the board are pure
* transport delay in the feedback path, and transport delay destroys phase
* margin no matter what the gains are. So instead of queueing, we refuse: the
* DAC holds its previous value for one more iteration and the *next* iteration
* sends a fresher command. Same drop-stale-keep-newest discipline asi_source
* uses on the camera queue.
*
* Cost model: cmd_time > 0 charges per command (the measured firmware service
* time -- the real bottleneck on a USB CDC bridge, where baud is nominal);
* otherwise falls back to baud-derived byte time.
*
* A side effect worth having: because we never overfill the buffer, the blocking
* write() below cannot stall the loop. */
static bool tx_admit(struct aylp_piplate_bridge_data *data, size_t ncmds,
	size_t len)
{
	double now = monotonic_s();
	if (UNLIKELY(!data->tx_primed)) {
		data->tx_free = now;
		data->tx_primed = true;
	}
	// the board went idle while we weren't writing
	if (data->tx_free < now) data->tx_free = now;
	if (data->tx_free - now > data->max_backlog) return false;
	data->tx_free += data->cmd_time > 0.0 ? ncmds * data->cmd_time
		: len * data->byte_time;
	return true;
}

/** write(2) that copes with a short write, which would otherwise leave half a
* setDAC command on the wire and desynchronize the board's parser. */
static int write_all(int fd, const char *buf, size_t len)
{
	size_t off = 0;
	while (off < len) {
		ssize_t w = write(fd, buf + off, len - off);
		if (w < 0) {
			if (errno == EINTR) continue;
			return -1;
		}
		off += (size_t)w;
	}
	return 0;
}

static void drain_response(int fd)
{
	char c;
	while (read(fd, &c, 1) == 1 && c != '\n')
		;
}

static int send_cmd(int fd, const char *cmd)
{
	size_t n = strlen(cmd);
	if (write(fd, cmd, n) != (ssize_t)n) return -1;
	drain_response(fd);
	return 0;
}

/** Read one '\n'-terminated ack. Returns 1 on a full ack, 0 if the board went
* silent (VTIME read timeout, 500 ms). */
static int read_ack(int fd)
{
	char c;
	while (read(fd, &c, 1) == 1)
		if (c == '\n') return 1;
	return 0;
}

static int cmp_dbl(const void *a, const void *b)
{
	double x = *(const double *)a, y = *(const double *)b;
	return (x > y) - (x < y);
}

/** Measure the board's true per-setDAC service time by timing acked commands.
*
* Rotates through the configured channels writing each one's `offset` voltage
* (the same bias the startup hold parks them at, so nothing moves that wasn't
* about to). Warm-up round-trips are discarded, the median of the rest is
* returned. 0.0 on failure (write error or no ack), in which case the caller
* should fall back to byte-time pacing. */
static double calibrate_cmd_time(struct aylp_piplate_bridge_data *data)
{
	double samples[PIPLATE_CALIB_SAMPLES];
	char cmd[64];
	int total = PIPLATE_CALIB_WARMUP + PIPLATE_CALIB_SAMPLES;
	for (int k = 0; k < total; k++) {
		size_t i = (size_t)k % data->n_outputs;
		int code = volts_to_code(data->offsets[i]);
		int len = snprintf(cmd, sizeof cmd, "DAQC2.setDAC(%d,%d,%.3f)\n",
			data->board, data->channels[i], code * DAC_LSB);
		double t0 = monotonic_s();
		if (write_all(data->serial_fd, cmd, (size_t)len)) return 0.0;
		if (!read_ack(data->serial_fd)) return 0.0;
		if (k >= PIPLATE_CALIB_WARMUP)
			samples[k - PIPLATE_CALIB_WARMUP] = monotonic_s() - t0;
	}
	qsort(samples, PIPLATE_CALIB_SAMPLES, sizeof *samples, cmp_dbl);
	return samples[PIPLATE_CALIB_SAMPLES / 2];
}

// helpers for reading a param that may be either a scalar or an array, so a
// single stage can carry per-channel values. An absent value yields def; an
// array index past the end also yields def.
static size_t val_count(struct json_object *v)
{
	if (v && json_object_is_type(v, json_type_array))
		return json_object_array_length(v);
	return 1;
}

static int val_int_at(struct json_object *v, size_t i, int def)
{
	if (!v) return def;
	if (json_object_is_type(v, json_type_array)) {
		struct json_object *e = json_object_array_get_idx(v, i);
		return e ? (int)json_object_get_int64(e) : def;
	}
	return (int)json_object_get_int64(v);
}

static double val_dbl_at(struct json_object *v, size_t i, double def)
{
	if (!v) return def;
	if (json_object_is_type(v, json_type_array)) {
		struct json_object *e = json_object_array_get_idx(v, i);
		return e ? json_object_get_double(e) : def;
	}
	return json_object_get_double(v);
}

int piplate_bridge_init(struct aylp_device *self)
{
	self->proc = &piplate_bridge_proc;
	self->fini = &piplate_bridge_fini;
	struct aylp_piplate_bridge_data *data = xcalloc(1, sizeof *data);
	self->device_data = data;

	data->serial_fd   = -1;
	data->board       = 0;
	data->channel      = 1;
	data->wait_response = true;
	data->skip_unchanged = false;
	data->baud         = 115200;
	// 0 => admit a write only when the line is idle, which holds the queued
	// command latency down to a single command's wire time
	data->max_backlog  = 0.0;
	data->fg_type      = 1;
	data->fg_frequency = 1000;
	data->fg_level     = 4;

	const char *port = "/dev/ttyACM0";
	const char *mode = "dac";

	// DAC-mode params may be scalars or arrays; hold onto the json objects and
	// expand them into per-channel arrays once we know how many outputs there are
	struct json_object *ch_val = NULL, *idx_val = NULL;
	struct json_object *scale_val = NULL, *offset_val = NULL;
	struct json_object *delay_val = NULL;

	if (self->params) {
		json_object_object_foreach(self->params, key, val) {
			if (key[0] == '_') {
			} else if (!strcmp(key, "port")) {
				port = json_object_get_string(val);
			} else if (!strcmp(key, "mode")) {
				mode = json_object_get_string(val);
			} else if (!strcmp(key, "board")) {
				data->board = (int)json_object_get_int64(val);
			} else if (!strcmp(key, "channel")) {
				ch_val = val;
				data->channel = (int)json_object_get_int64(val);
			} else if (!strcmp(key, "index")) {
				idx_val = val;
			} else if (!strcmp(key, "scale")) {
				scale_val = val;
			} else if (!strcmp(key, "offset")) {
				offset_val = val;
			} else if (!strcmp(key, "start_delay")) {
				delay_val = val;
			} else if (!strcmp(key, "wait_response")) {
				data->wait_response = json_object_get_boolean(val);
			} else if (!strcmp(key, "skip_unchanged")) {
				data->skip_unchanged = json_object_get_boolean(val);
			} else if (!strcmp(key, "baud")) {
				data->baud = (long)json_object_get_int64(val);
			} else if (!strcmp(key, "cmd_time_us")) {
				data->cmd_time = json_object_get_double(val) * 1e-6;
			} else if (!strcmp(key, "max_backlog_ms")) {
				data->max_backlog = json_object_get_double(val) * 1e-3;
			} else if (!strcmp(key, "frequency")) {
				data->fg_frequency = (int)json_object_get_int64(val);
			} else if (!strcmp(key, "type")) {
				data->fg_type = (int)json_object_get_int64(val);
			} else if (!strcmp(key, "level")) {
				data->fg_level = (int)json_object_get_int64(val);
			} else {
				log_warn("piplate_bridge: unknown param \"%s\"", key);
			}
		}
	}

	data->use_fg = !strcmp(mode, "fg");

	speed_t speed = baud_to_speed(data->baud);
	if (!speed) {
		log_error("piplate_bridge: unsupported baud %ld (use 9600, 19200, "
			"38400, 57600, 115200, 230400, 460800 or 921600)",
			data->baud);
		return -1;
	}
	if (data->max_backlog < 0.0) {
		log_error("piplate_bridge: max_backlog_ms must be >= 0");
		return -1;
	}
	data->byte_time = 10.0 / (double)data->baud;	// 8N1 => 10 bits/byte

	// expand the DAC-mode params into per-channel arrays. The number of outputs
	// is the longest of the channel/index/scale/offset arrays; any param given
	// as a scalar is broadcast to every channel.
	if (!data->use_fg) {
		size_t n = val_count(ch_val);
		size_t counts[] = {
			val_count(idx_val), val_count(scale_val),
			val_count(offset_val), val_count(delay_val)
		};
		for (size_t k = 0; k < 4; k++)
			if (counts[k] > n) n = counts[k];

		data->n_outputs = n;
		data->channels   = xcalloc(n, sizeof *data->channels);
		data->indices    = xcalloc(n, sizeof *data->indices);
		data->scales     = xcalloc(n, sizeof *data->scales);
		data->offsets    = xcalloc(n, sizeof *data->offsets);
		data->last_codes = xcalloc(n, sizeof *data->last_codes);
		data->pending_codes = xcalloc(n, sizeof *data->pending_codes);
		data->start_delays = xcalloc(n, sizeof *data->start_delays);
		data->released     = xcalloc(n, sizeof *data->released);
		for (size_t i = 0; i < n; i++) {
			data->channels[i] = val_int_at(ch_val, i, 0);
			data->indices[i]  = val_int_at(idx_val, i, 0);
			data->scales[i]   = val_dbl_at(scale_val, i, 1.0);
			data->offsets[i]  = val_dbl_at(offset_val, i, 0.0);
			data->start_delays[i] = val_dbl_at(delay_val, i, 0.0);
			if (data->start_delays[i] < 0.0) {
				log_error("piplate_bridge: start_delay must be >= 0");
				return -1;
			}
			if (data->start_delays[i] > 0.0) data->has_start_delay = true;
			data->last_codes[i] = -1;  // force a write on iteration 0
		}
		// one setDAC command is at most PIPLATE_CMD_BYTES chars; 64 per
		// channel is plenty
		data->cmdbuf_sz = 64 * n;
		data->cmdbuf = xmalloc(data->cmdbuf_sz);
	}

	data->serial_fd = serial_open(port, speed);
	if (data->serial_fd < 0) {
		log_error("piplate_bridge: open %s: %s", port, strerror(errno));
		return -1;
	}

	if (data->use_fg) {
		char cmd[64];
		snprintf(cmd, sizeof cmd, "DAQC2.fgTYPE(%d, %d, %d)\n",
			data->board, data->channel, data->fg_type);
		if (send_cmd(data->serial_fd, cmd)) goto serial_err;
		snprintf(cmd, sizeof cmd, "DAQC2.fgFREQ(%d, %d, %d)\n",
			data->board, data->channel, data->fg_frequency);
		if (send_cmd(data->serial_fd, cmd)) goto serial_err;
		snprintf(cmd, sizeof cmd, "DAQC2.fgLEVEL(%d, %d, %d)\n",
			data->board, data->channel, data->fg_level);
		if (send_cmd(data->serial_fd, cmd)) goto serial_err;
		snprintf(cmd, sizeof cmd, "DAQC2.fgON(%d, %d)\n",
			data->board, data->channel);
		if (send_cmd(data->serial_fd, cmd)) goto serial_err;

		self->type_in   = AYLP_T_ANY;
		self->units_in  = AYLP_U_ANY;
		self->type_out  = AYLP_T_UNCHANGED;
		self->units_out = AYLP_U_UNCHANGED;
		log_info("piplate_bridge: %s board=%d fg_channel=%d type=%d freq=%dHz level=%d",
			port, data->board, data->channel,
			data->fg_type, data->fg_frequency, data->fg_level);
	} else {
		// stop any FG that may still be running on these channels from a
		// previous run
		for (size_t i = 0; i < data->n_outputs; i++) {
			char cmd[64];
			snprintf(cmd, sizeof cmd, "DAQC2.fgOFF(%d, %d)\n",
				data->board, data->channels[i]);
			send_cmd(data->serial_fd, cmd);
		}

		// pace on the board's real per-command service time, not the
		// nominal baud (fictional over USB CDC): measure it unless the
		// config supplied cmd_time_us (>0 fixed, <0 legacy byte pacing)
		if (data->cmd_time == 0.0) {
			double t = calibrate_cmd_time(data);
			if (t > 0.0) {
				data->cmd_time = t;
				log_info("piplate_bridge: measured %.0f us/command"
					" (ceiling %.0f commands/s)",
					1e6 * t, 1.0 / t);
			} else {
				data->cmd_time = -1.0;
				log_warn("piplate_bridge: cmd_time calibration "
					"got no ack; pacing on %ld baud byte "
					"time instead", data->baud);
			}
		}

		self->type_in   = AYLP_T_VECTOR;
		self->units_in  = AYLP_U_MINMAX;
		self->type_out  = AYLP_T_UNCHANGED;
		self->units_out = AYLP_U_UNCHANGED;
		log_info("piplate_bridge: %s board=%d %zu output(s) wait_response=%d",
			port, data->board, data->n_outputs, data->wait_response);
		for (size_t i = 0; i < data->n_outputs; i++)
			log_info("piplate_bridge:   channel=%d index=%d "
				"scale=%g offset=%g", data->channels[i],
				data->indices[i], data->scales[i], data->offsets[i]);
	}
	return 0;

serial_err:
	log_error("piplate_bridge: serial command failed: %s", strerror(errno));
	return -1;
}

int piplate_bridge_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_piplate_bridge_data *data = self->device_data;

	if (data->use_fg)
		return 0;

	// in fire-and-forget mode, discard last iteration's acks so the input
	// buffer stays bounded (we never read them otherwise)
	if (!data->wait_response)
		tcflush(data->serial_fd, TCIFLUSH);

	double now = monotonic_s();
	// Measure the startup hold from the first proc, not from init: this is the
	// iteration on which every channel first gets written (last_codes is -1), so
	// it is when the coarse channels actually reach their bias voltages.
	if (UNLIKELY(!data->t0)) data->t0 = now;
	double elapsed = now - data->t0;

	// build one buffer of setDAC commands. By default every channel is
	// refreshed each iteration (loaded outputs droop and must be re-sent); with
	// skip_unchanged, a channel whose DAC code is unchanged is omitted.
	size_t len = 0;
	size_t written = 0;
	for (size_t i = 0; i < data->n_outputs; i++) {
		data->pending_codes[i] = -1;
		int idx = data->indices[i];
		if (UNLIKELY((size_t)idx >= state->vector->size)) {
			log_error("piplate_bridge: index %d out of range "
				"(vector size %zu)", idx, state->vector->size);
			break;
		}
		double cmd = gsl_vector_get(state->vector, idx);
		// Startup hold: park this channel at its bias (cmd = 0) until its
		// start_delay has elapsed. A coarse stage needs time to walk the beam
		// near centre; a fine channel that chases the large startup error
		// meanwhile just drives itself into its rail and has no travel left by
		// the time the coarse one arrives. Note this only parks the *output* --
		// the upstream pid keeps integrating unless you give it a matching
		// start_delay, in which case it commands zero and this hold is belt and
		// braces for staggering channels.
		if (UNLIKELY(data->has_start_delay && !data->released[i])) {
			if (elapsed < data->start_delays[i]) {
				cmd = 0.0;
			} else {
				data->released[i] = true;
				if (data->start_delays[i] > 0.0)
					log_info("piplate_bridge: channel %d released "
						"after %.2f s hold", data->channels[i],
						elapsed);
			}
		}
		int code = volts_to_code(data->offsets[i] + cmd * data->scales[i]);

		if (data->skip_unchanged && code == data->last_codes[i]) {
			// the DAC already holds this code; a finer command change
			// than 1 mV cannot move the output, so don't spend wire
			// time saying so
			data->diag_skips++;
			continue;
		}
		// stage it: last_codes is only advanced once the bytes are away, so
		// a dropped iteration re-sends rather than silently forgetting
		data->pending_codes[i] = code;

		len += (size_t)snprintf(data->cmdbuf + len, data->cmdbuf_sz - len,
			"DAQC2.setDAC(%d,%d,%.3f)\n", data->board,
			data->channels[i], code * DAC_LSB);
		written++;
		log_trace("piplate_bridge: ch%d cmd=%g → code %d (%.3f V)",
			data->channels[i], cmd, code, code * DAC_LSB);
	}

	if (len) {
		// wait_response mode self-limits: drain_response blocks for the
		// board's ack, so nothing accumulates. Fire-and-forget has no such
		// backpressure, so the transmit model provides it.
		if (!data->wait_response && !tx_admit(data, written, len)) {
			data->diag_drops++;
			return 0;	// pending_codes discarded; retry next iteration
		}
		if (UNLIKELY(write_all(data->serial_fd, data->cmdbuf, len))) {
			log_warn("piplate_bridge: serial write failed: %s",
				strerror(errno));
			return 0;	// don't commit: the DAC didn't get these
		}
		// the bytes are away; now it's safe to remember them
		for (size_t i = 0; i < data->n_outputs; i++)
			if (data->pending_codes[i] >= 0)
				data->last_codes[i] = data->pending_codes[i];
		data->diag_writes++;
		data->diag_cmds += (long)written;
		// drain exactly one ack per command only when asked to; fire-and-
		// forget skips this and clears the acks via tcflush next iteration
		if (data->wait_response)
			for (size_t i = 0; i < written; i++)
				drain_response(data->serial_fd);
	}

	// report throughput: if drops dominate, the wire is the bottleneck and the
	// loop is running faster than the actuator can be commanded
	if (UNLIKELY(!data->diag_t0)) data->diag_t0 = now;
	double dt = now - data->diag_t0;
	if (UNLIKELY(dt >= PIPLATE_DIAG_PERIOD)) {
		long total = data->diag_writes + data->diag_drops;
		double cap = data->cmd_time > 0.0 ? 1.0 / data->cmd_time
			: data->baud / (10.0 * PIPLATE_CMD_BYTES);
		log_info("piplate_bridge: %.0f writes/s | %.0f%% dropped (board "
			"busy) | %.0f skips/s (code unchanged) | %.0f cmds/s "
			"of %.0f",
			data->diag_writes / dt,
			total ? 100.0 * data->diag_drops / total : 0.0,
			data->diag_skips / dt,
			data->diag_cmds / dt, cap);
		// the voltage each DAC channel is currently holding (last code
		// actually written; "--" means nothing has reached it yet)
		char vbuf[128];
		size_t voff = 0;
		for (size_t i = 0; i < data->n_outputs
				&& voff < sizeof vbuf; i++) {
			if (data->last_codes[i] < 0)
				voff += (size_t)snprintf(vbuf + voff,
					sizeof vbuf - voff, " ch%d=--",
					data->channels[i]);
			else
				voff += (size_t)snprintf(vbuf + voff,
					sizeof vbuf - voff, " ch%d=%.3fV",
					data->channels[i],
					data->last_codes[i] * DAC_LSB);
		}
		log_info("piplate_bridge: DAC%s", vbuf);
		data->diag_t0 = now;
		data->diag_writes = data->diag_drops = data->diag_skips = 0;
		data->diag_cmds = 0;
	}

	return 0;
}

int piplate_bridge_fini(struct aylp_device *self)
{
	struct aylp_piplate_bridge_data *data = self->device_data;
	if (data->serial_fd >= 0) {
		if (data->use_fg) {
			char cmd[64];
			snprintf(cmd, sizeof cmd, "DAQC2.fgOFF(%d, %d)\n",
				data->board, data->channel);
			send_cmd(data->serial_fd, cmd);
		}
		close(data->serial_fd);
	}
	if (data->channels) xfree(data->channels);
	if (data->indices) xfree(data->indices);
	if (data->scales) xfree(data->scales);
	if (data->offsets) xfree(data->offsets);
	if (data->last_codes) xfree(data->last_codes);
	if (data->pending_codes) xfree(data->pending_codes);
	if (data->start_delays) xfree(data->start_delays);
	if (data->released) xfree(data->released);
	if (data->cmdbuf) xfree(data->cmdbuf);
	xfree(data);
	return 0;
}
