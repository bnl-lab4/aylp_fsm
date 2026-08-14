#include <errno.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <unistd.h>
#include <json-c/json.h>
#include "anyloop.h"
#include "logging.h"
#include "udp_source.h"
#include "xalloc.h"


int udp_source_init(struct aylp_device *self)
{
	self->proc = &udp_source_proc;
	self->fini = &udp_source_fini;
	self->device_data = xcalloc(1, sizeof(struct aylp_udp_source_data));
	struct aylp_udp_source_data *data = self->device_data;

	unsigned short port = 0;
	size_t n = 0;

	if (!self->params) {
		log_error("No params object found.");
		return -1;
	}
	json_object_object_foreach(self->params, key, val) {
		if (key[0] == '_') {
			// keys starting with _ are comments
		} else if (!strcmp(key, "port")) {
			port = (unsigned short)strtoul(
				json_object_get_string(val), 0, 0
			);
			log_trace("port = %hu", port);
		} else if (!strcmp(key, "n")) {
			n = (size_t)json_object_get_uint64(val);
			log_trace("n = %zu", n);
		} else if (!strcmp(key, "units")) {
			const char *s = json_object_get_string(val);
			data->units = aylp_units_from_string(s);
			log_trace("units = %s (0x%hhX)", s, data->units);
		} else {
			log_warn("Unknown parameter \"%s\"", key);
		}
	}

	if (!port || !n) {
		log_error("You must provide all params: port, n.");
		return -1;
	}
	if (!data->units) {
		log_error("You must provide a valid units param.");
		return -1;
	}

	// zero-initialized output vector (holds last received command)
	data->vec = gsl_vector_calloc(n);

	// non-blocking UDP socket bound to loopback
	data->sock = socket(AF_INET, SOCK_DGRAM | SOCK_NONBLOCK, 0);
	if (data->sock < 0) {
		log_error("Couldn't create socket: %s", strerror(errno));
		return -1;
	}

	// keep receive buffer small so stale commands are dropped
	int rcvbuf = 4096;
	if (setsockopt(data->sock, SOL_SOCKET, SO_RCVBUF,
			&rcvbuf, sizeof(rcvbuf))) {
		log_warn("Couldn't set SO_RCVBUF: %s", strerror(errno));
	}

	memset(&data->sa, 0, sizeof(data->sa));
	data->sa.sin_family = AF_INET;
	data->sa.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
	data->sa.sin_port = htons(port);
	if (bind(data->sock, (struct sockaddr *)&data->sa, sizeof(data->sa))) {
		log_error("Couldn't bind to port %hu: %s", port, strerror(errno));
		return -1;
	}

	// accepts anything before it (or nothing, as a pipeline source)
	self->type_in = AYLP_T_ANY;
	self->units_in = AYLP_U_ANY;
	self->type_out = AYLP_T_VECTOR;
	self->units_out = data->units;
	return 0;
}


int udp_source_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_udp_source_data *data = self->device_data;
	size_t n = data->vec->size;
	size_t expected = n * sizeof(double);

	// non-blocking receive; if nothing is ready, keep the last value
	ssize_t r = recv(data->sock, data->vec->data, expected, 0);
	if (r < 0) {
		if (errno != EAGAIN && errno != EWOULDBLOCK) {
			log_error("recv error: %s", strerror(errno));
			return -1;
		}
		// no packet ready — fall through with last value
	} else if ((size_t)r != expected) {
		log_warn("udp_source: received %zd bytes, expected %zu; ignoring",
			r, expected);
	} else {
		log_trace("udp_source: received %zd bytes", r);
	}

	state->vector = data->vec;
	state->header.type = AYLP_T_VECTOR;
	state->header.units = data->units;
	state->header.log_dim.y = n;
	state->header.log_dim.x = 1;
	state->header.pitch.y = 0.0;
	state->header.pitch.x = 0.0;
	return 0;
}


int udp_source_fini(struct aylp_device *self)
{
	struct aylp_udp_source_data *data = self->device_data;
	close(data->sock);
	xfree_type(gsl_vector, data->vec);
	xfree(data);
	return 0;
}
