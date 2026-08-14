// Thorlabs KDC101 K-Cube DC servo controller -- anyloop device.
// Drives one axis of the slow steering mirror via the Thorlabs APT serial
// protocol over FTDI USB-serial (115200 8N1, RTS/CTS hardware flow control).
//
// Params:
//   port    -- serial device path, e.g. "/dev/ttyUSB0"
//   index   -- pipeline vector element to command (default 0, i.e. pitch/y)
//   scale   -- encoder counts per pipeline unit (e.g. 21000 for a Z8 actuator
//              spanning ±1 MINMAX unit)
//   center  -- encoder counts at zero command (parked / centre position)

#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <termios.h>
#include <unistd.h>

#include "anyloop.h"
#include "kdc101.h"
#include "logging.h"
#include "xalloc.h"

// ── Thorlabs APT serial protocol ─────────────────────────────────────────────

#define MGMSG_MOD_SET_CHANENABLESTATE  0x0210u
#define MGMSG_MOT_MOVE_ABSOLUTE        0x0453u

#define APT_HOST  0x01u   // source: host PC
#define APT_DEST  0x50u   // dest:   generic USB K-Cube bay

// Pack a 6-byte short APT message (header only, no data payload).
static void apt_short(uint8_t buf[6], uint16_t id, uint8_t p1, uint8_t p2)
{
	buf[0] = id & 0xFFu;
	buf[1] = id >> 8;
	buf[2] = p1;
	buf[3] = p2;
	buf[4] = APT_DEST;
	buf[5] = APT_HOST;
}

// Pack the 6-byte header of a long APT message (data_len payload bytes follow).
static void apt_long_hdr(uint8_t buf[6], uint16_t id, uint16_t data_len)
{
	buf[0] = id & 0xFFu;
	buf[1] = id >> 8;
	buf[2] = data_len & 0xFFu;
	buf[3] = data_len >> 8;
	buf[4] = APT_DEST | 0x80u;  // MSB set signals long message
	buf[5] = APT_HOST;
}

// Blocking write that retries on EINTR.
static int apt_write(int fd, const void *buf, size_t n)
{
	const uint8_t *p = buf;
	while (n) {
		ssize_t r = write(fd, p, n);
		if (r < 0) {
			if (errno == EINTR) continue;
			return -1;
		}
		p += (size_t)r;
		n -= (size_t)r;
	}
	return 0;
}

// Send MGMSG_MOT_MOVE_ABSOLUTE for channel 1 to an absolute encoder position.
static int kdc101_move_absolute(int fd, int32_t counts)
{
	uint8_t msg[12];
	apt_long_hdr(msg, MGMSG_MOT_MOVE_ABSOLUTE, 6);
	msg[6]  = 0x01; msg[7]  = 0x00;       // ChanIdent = 1
	msg[8]  = (uint8_t)(counts);
	msg[9]  = (uint8_t)(counts >> 8);
	msg[10] = (uint8_t)(counts >> 16);
	msg[11] = (uint8_t)(counts >> 24);
	return apt_write(fd, msg, sizeof msg);
}

// ── anyloop interface ─────────────────────────────────────────────────────────

int kdc101_init(struct aylp_device *self)
{
	self->proc = &kdc101_proc;
	self->fini = &kdc101_fini;
	struct aylp_kdc101_data *data = xcalloc(1, sizeof *data);
	self->device_data = data;

	// defaults
	data->index  = 0;
	data->scale  = 1000.0;
	data->center = 0;
	const char *port = NULL;

	if (!self->params) {
		log_error("kdc101: no params object");
		return -1;
	}
	json_object_object_foreach(self->params, key, val) {
		if (key[0] == '_') {
		} else if (!strcmp(key, "port")) {
			port = json_object_get_string(val);
			log_trace("kdc101: port = %s", port);
		} else if (!strcmp(key, "index")) {
			data->index = (int)json_object_get_int64(val);
			log_trace("kdc101: index = %d", data->index);
		} else if (!strcmp(key, "scale")) {
			data->scale = json_object_get_double(val);
			log_trace("kdc101: scale = %g", data->scale);
		} else if (!strcmp(key, "center")) {
			data->center = (int32_t)json_object_get_int64(val);
			log_trace("kdc101: center = %d", (int)data->center);
		} else {
			log_warn("kdc101: unknown param \"%s\"", key);
		}
	}
	if (!port) {
		log_error("kdc101: \"port\" param is required");
		return -1;
	}

	// Open the FTDI serial port.
	// APT requires 115200 8N1 with RTS/CTS hardware flow control.
	data->fd = open(port, O_RDWR | O_NOCTTY);
	if (data->fd < 0) {
		log_error("kdc101: open(%s): %s", port, strerror(errno));
		return -1;
	}
	struct termios tty;
	if (tcgetattr(data->fd, &tty)) {
		log_error("kdc101: tcgetattr: %s", strerror(errno));
		return -1;
	}
	cfmakeraw(&tty);
	cfsetispeed(&tty, B115200);
	cfsetospeed(&tty, B115200);
	tty.c_cflag |= CRTSCTS | CLOCAL | CREAD;
	tty.c_cflag &= (tcflag_t)~CSTOPB;
	tty.c_cc[VMIN]  = 0;
	tty.c_cc[VTIME] = 0;  // non-blocking reads (we never read in proc)
	if (tcsetattr(data->fd, TCSANOW, &tty)) {
		log_error("kdc101: tcsetattr: %s", strerror(errno));
		return -1;
	}
	tcflush(data->fd, TCIOFLUSH);

	// Enable channel 1.
	uint8_t enable[6];
	apt_short(enable, MGMSG_MOD_SET_CHANENABLESTATE, 0x01, 0x01);
	if (apt_write(data->fd, enable, 6)) {
		log_error("kdc101: channel enable: %s", strerror(errno));
		return -1;
	}

	self->type_in   = AYLP_T_VECTOR;
	self->units_in  = AYLP_U_ANY;
	self->type_out  = AYLP_T_UNCHANGED;
	self->units_out = AYLP_U_UNCHANGED;

	log_info("kdc101: open %s  index=%d  scale=%g  center=%d",
		port, data->index, data->scale, (int)data->center);
	return 0;
}

int kdc101_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_kdc101_data *data = self->device_data;

	if (UNLIKELY((size_t)data->index >= state->vector->size)) {
		log_error("kdc101: index %d out of range (vector size %zu)",
			data->index, state->vector->size);
		return -1;
	}
	double cmd = gsl_vector_get(state->vector, data->index);
	int32_t counts = data->center + (int32_t)(cmd * data->scale);

	if (UNLIKELY(kdc101_move_absolute(data->fd, counts)))
		log_warn("kdc101: write failed: %s", strerror(errno));

	log_trace("kdc101: cmd=%g → %d counts", cmd, (int)counts);
	return 0;
}

int kdc101_fini(struct aylp_device *self)
{
	struct aylp_kdc101_data *data = self->device_data;
	if (data->fd >= 0) close(data->fd);
	xfree(data);
	return 0;
}
