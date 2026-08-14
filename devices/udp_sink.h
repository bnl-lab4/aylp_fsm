#ifndef AYLP_DEVICES_UDP_SINK_H_
#define AYLP_DEVICES_UDP_SINK_H_

#include <netinet/in.h>
#include <sys/socket.h>
#include <sys/uio.h>
#include "anyloop.h"

// what to send: the data itself, or a reduction of it
enum aylp_udp_sink_reduce {
	// send the whole block, unchanged (the default)
	AYLP_UDP_SINK_RAW = 0,
	// send [peak, mean, saturated_fraction] as three doubles instead
	AYLP_UDP_SINK_STATS,
};

struct aylp_udp_sink_data {
	int sock;
	// destination
	struct sockaddr_in dest_sa;
	// iovec for writev so we can write the header and the data in one go
	struct iovec iovecs[2];
	// gsl_block_uchar that we will copy pointer to data to
	gsl_block_uchar bytes;
	// send only every decimation'th iteration (1 = send every one)
	size_t decimation;
	// counts iterations since the last send
	size_t countdown;
	// param: whether to reduce the data before sending it
	enum aylp_udp_sink_reduce reduce;
	// reduction output and its header, reused every iteration
	double stats[3];
	struct aylp_header stats_head;
	// Optional nearest-neighbour image preview. A 384x384 RAW8 frame cannot
	// fit in one UDP datagram; factor 2 produces 192x192 while preserving the
	// complete field of view and normalized centroid coordinates.
	size_t downsample;
	unsigned char *preview;
	size_t preview_size;
	struct aylp_header preview_head;
};

// initialize udp_sink device
int udp_sink_init(struct aylp_device *self);

// process udp_sink device once per loop
int udp_sink_proc(struct aylp_device *self, struct aylp_state *state);

// close udp_sink device when loop exits
int udp_sink_fini(struct aylp_device *self);

#endif
