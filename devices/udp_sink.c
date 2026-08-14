#include <errno.h>
#include <limits.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <json-c/json.h>
#include <netinet/in.h>
#include <unistd.h>
#include <sys/socket.h>
#include <sys/uio.h>
#include "anyloop.h"
#include "block.h"
#include "logging.h"
#include "udp_sink.h"
#include "xalloc.h"


/** Reduce the pipeline data to [peak, mean, saturated_fraction].
* Written to out as three doubles. This exists so a viewer can follow beam
* intensity at the FULL loop rate: a 248x248 uchar frame is 61.5 kB, which is
* 233 MB/s at 3788 Hz and cannot be sent (or decoded) every iteration, while
* this is 24 bytes and costs one pass over data that is already hot in cache.
* Send the raw frame instead when you need the image itself -- the reduction is
* lossy by construction and there is no way back to the pixels from it.
*
* saturated_fraction is the share of pixels at the maximum the sample type can
* represent, which is meaningful only for uchar (camera) data; it is reported as
* 0 for double data, where there is no such ceiling. It matters because a peak
* pinned at 255 hides everything above it, so the mean is the only honest
* brightness reading left once the sensor clips.
*
* Walks the data in its native layout rather than calling
* get_contiguous_bytes(), which would copy the whole block for a non-contiguous
* matrix -- a malloc and a 61.5 kB memcpy per iteration is exactly the cost this
* is here to avoid. */
static int udp_sink_stats(struct aylp_state *state, double *out)
{
	size_t n = 0;
	// integer accumulators on the uchar path: comparing and summing into
	// doubles inside the loop stops the compiler from vectorizing it
	// 32 bits is enough for both and keeps the widening the vectorizer has
	// to do to a minimum: a 248x248 frame of 255s sums to 15.7e6, and the
	// saturated count can't exceed the pixel count
	unsigned char peak_u = 0;
	unsigned sum_u = 0, sat_u = 0;
	double peak_d = 0.0, sum_d = 0.0;
	switch (state->header.type) {
	case AYLP_T_MATRIX_UCHAR: {
		gsl_matrix_uchar *m = state->matrix_uchar;
		for (size_t i = 0; i < m->size1; i++) {
			const unsigned char *row = m->data + i*m->tda;
			for (size_t j = 0; j < m->size2; j++) {
				unsigned char v = row[j];
				if (v > peak_u) peak_u = v;
				sum_u += v;
				sat_u += (v == UCHAR_MAX);
			}
		}
		n = m->size1 * m->size2;
		out[0] = peak_u;
		out[1] = n ? (double)sum_u / n : 0.0;
		out[2] = n ? (double)sat_u / n : 0.0;
		return 0; }
	case AYLP_T_BLOCK_UCHAR: {
		gsl_block_uchar *b = state->block_uchar;
		for (size_t i = 0; i < b->size; i++) {
			unsigned char v = b->data[i];
			if (v > peak_u) peak_u = v;
			sum_u += v;
			sat_u += (v == UCHAR_MAX);
		}
		n = b->size;
		out[0] = peak_u;
		out[1] = n ? (double)sum_u / n : 0.0;
		out[2] = n ? (double)sat_u / n : 0.0;
		return 0; }
	case AYLP_T_MATRIX: {
		gsl_matrix *m = state->matrix;
		for (size_t i = 0; i < m->size1; i++) {
			const double *row = m->data + i*m->tda;
			for (size_t j = 0; j < m->size2; j++) {
				if (row[j] > peak_d) peak_d = row[j];
				sum_d += row[j];
			}
		}
		n = m->size1 * m->size2;
		break; }
	case AYLP_T_VECTOR: {
		gsl_vector *v = state->vector;
		for (size_t i = 0; i < v->size; i++) {
			double e = v->data[i*v->stride];
			if (e > peak_d) peak_d = e;
			sum_d += e;
		}
		n = v->size;
		break; }
	case AYLP_T_BLOCK: {
		gsl_block *b = state->block;
		for (size_t i = 0; i < b->size; i++) {
			if (b->data[i] > peak_d) peak_d = b->data[i];
			sum_d += b->data[i];
		}
		n = b->size;
		break; }
	default:
		log_error("Can't reduce data of type 0x%hhX", state->header.type);
		return -1;
	}
	out[0] = peak_d;
	out[1] = n ? sum_d / n : 0.0;
	out[2] = 0.0;	// no saturation ceiling for double data
	return 0;
}


int udp_sink_init(struct aylp_device *self)
{
	self->proc = &udp_sink_proc;
	self->fini = &udp_sink_fini;
	self->device_data = xcalloc(1, sizeof(struct aylp_udp_sink_data));
	struct aylp_udp_sink_data *data = self->device_data;

	data->sock = socket(AF_INET, SOCK_DGRAM, 0);
	if (data->sock == -1) {
		log_error("Couldn't initialize socket: %s", strerror(errno));
		return -1;
	}
	memset(&(data->dest_sa), 0, sizeof(data->dest_sa));
	data->dest_sa.sin_family = AF_INET;
	data->decimation = 1;
	data->downsample = 1;
	int got_params = 0;	// use this to check for params in case ip is 0
	if (!self->params) {
		log_error("No params object found.");
		return -1;
	}
	json_object_object_foreach(self->params, key, val) {
		if (key[0] == '_') {
			// keys starting with _ are comments
		} else if (!strcmp(key, "ip")) {
			inet_pton(AF_INET,
				json_object_get_string(val),
				&(data->dest_sa.sin_addr)
			);
			got_params |= 1;
			log_trace("ip = 0x%X", data->dest_sa.sin_addr);
		} else if (!strcmp(key, "port")) {
			data->dest_sa.sin_port = htons(
				(unsigned short)strtoul(
					json_object_get_string(val), 0, 0
				)
			);
			got_params |= 2;
			log_trace("port = 0x%X", data->dest_sa.sin_port);
		} else if (!strcmp(key, "decimation")) {
			data->decimation = json_object_get_uint64(val);
			log_trace("decimation = %zu", data->decimation);
		} else if (!strcmp(key, "downsample")) {
			data->downsample = json_object_get_uint64(val);
			log_trace("downsample = %zu", data->downsample);
		} else if (!strcmp(key, "reduce")) {
			const char *s = json_object_get_string(val);
			if (!s || !strcmp(s, "none")) {
				data->reduce = AYLP_UDP_SINK_RAW;
			} else if (!strcmp(s, "stats")) {
				data->reduce = AYLP_UDP_SINK_STATS;
			} else {
				log_error("Unknown reduce \"%s\"; expected "
					"\"none\" or \"stats\"", s);
				return -1;
			}
			log_trace("reduce = %d", data->reduce);
		} else {
			log_warn("Unknown parameter \"%s\"", key);
		}
	}
	// make sure we didn't miss any params
	if (got_params != (1|2)) {
		log_error("You must provide all params: ip, port.");
		return -1;
	}
	if (!data->decimation) {
		log_error("decimation must be nonzero.");
		return -1;
	}
	if (!data->downsample) {
		log_error("downsample must be nonzero.");
		return -1;
	}
	// we're not using sendto(), so we need to connect() the socket first
	int err = connect(
		data->sock,
		(struct sockaddr *)&data->dest_sa,
		sizeof(data->dest_sa)
	);
	if (err) {
		log_error("Couldn't connect: %s", strerror(errno));
		return -1;
	}
	// set types and units
	self->type_in = AYLP_T_ANY;
	self->units_in = AYLP_U_ANY;
	self->type_out = AYLP_T_UNCHANGED;
	self->units_out = AYLP_U_UNCHANGED;
	return 0;
}


int udp_sink_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_udp_sink_data *data = self->device_data;
	// bail before get_contiguous_bytes(), which may copy the whole block
	if (data->countdown) {
		data->countdown -= 1;
		return 0;
	}
	data->countdown = data->decimation - 1;
	if (data->reduce == AYLP_UDP_SINK_STATS) {
		if (udp_sink_stats(state, data->stats)) return -1;
		// keep the pipeline's status flags, but describe what we are
		// actually sending: a 3-vector in counts, not the original data
		data->stats_head = state->header;
		data->stats_head.type = AYLP_T_VECTOR;
		data->stats_head.units = AYLP_U_COUNTS;
		data->stats_head.log_dim.y = 3;
		data->stats_head.log_dim.x = 1;
		data->iovecs[0].iov_base = &data->stats_head;
		data->iovecs[0].iov_len = sizeof(data->stats_head);
		data->iovecs[1].iov_base = data->stats;
		data->iovecs[1].iov_len = sizeof(data->stats);
		ssize_t err = writev(data->sock, data->iovecs, 2);
		if (err < 0) {
			// A connected UDP socket reports ECONNREFUSED when no viewer is
			// listening. Monitoring is optional: its absence must neither fail
			// the pipeline nor synchronously spam stderr from the RT thread.
			if (errno == ECONNREFUSED)
				return 0;
			log_error("Couldn't send stats: %s", strerror(errno));
			return -1;
		}
		return 0;
	}
	if (data->downsample > 1 && state->header.type == AYLP_T_MATRIX_UCHAR) {
		gsl_matrix_uchar *m = state->matrix_uchar;
		size_t out_h = (m->size1 + data->downsample - 1) / data->downsample;
		size_t out_w = (m->size2 + data->downsample - 1) / data->downsample;
		size_t need = out_h * out_w;
		if (need > 65467) {
			log_error("Downsampled preview is still too large for UDP: %zu bytes", need);
			return -1;
		}
		if (need > data->preview_size) {
			data->preview = xrealloc(data->preview, need);
			data->preview_size = need;
		}
		for (size_t y = 0; y < out_h; y++)
			for (size_t x = 0; x < out_w; x++)
				data->preview[y*out_w + x] =
					gsl_matrix_uchar_get(m, y*data->downsample,
						x*data->downsample);
		data->preview_head = state->header;
		data->preview_head.log_dim.y = out_h;
		data->preview_head.log_dim.x = out_w;
		data->iovecs[0].iov_base = &data->preview_head;
		data->iovecs[0].iov_len = sizeof(data->preview_head);
		data->iovecs[1].iov_base = data->preview;
		data->iovecs[1].iov_len = need;
		ssize_t err = writev(data->sock, data->iovecs, 2);
		if (err < 0 && errno == ECONNREFUSED)
			return 0;
		if (err < 0)
			log_error("Couldn't send preview: %s", strerror(errno));
		return err < 0 ? -1 : 0;
	}
	// make data contiguous
	int needs_free = get_contiguous_bytes(&data->bytes, state);
	if (needs_free < 0) return needs_free;
	// first thing we send is the aylp header
	data->iovecs[0].iov_base = &state->header;
	data->iovecs[0].iov_len = sizeof(state->header);
	// second thing we send is the block data
	data->iovecs[1].iov_base = data->bytes.data;
	data->iovecs[1].iov_len = data->bytes.size;
	// write all the data in one go
	size_t n = data->iovecs[0].iov_len + data->iovecs[1].iov_len;
	log_trace("Writing %zu bytes to UDP", n);
	ssize_t err = writev(data->sock, data->iovecs, 2);
	if (err < 0 && errno == ECONNREFUSED) {
		if (needs_free)
			xfree(data->bytes.data);
		return 0;
	}
	if (err < 0) {
		// if n > SSIZE_MAX, this will fire, so we don't need to check
		// the sign of (ssize_t)n ourselves
		log_error("Couldn't send data: %s", strerror(errno));
	} else if (err != (ssize_t)n) {
		if (err == -1) {
			log_error("Couldn't send data: %s", strerror(errno));
		} else {
			log_error("Short write: %zu of %zu", err, n);
		}
	}
	if (needs_free) {
		xfree(data->bytes.data);
	}
	return (err < 0) ? -1 : 0;
}


int udp_sink_fini(struct aylp_device *self)
{
	struct aylp_udp_sink_data *data = self->device_data;
	xfree(data->preview);
	xfree(self->device_data);
	return 0;
}
