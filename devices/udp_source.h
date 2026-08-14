#ifndef AYLP_DEVICES_UDP_SOURCE_H_
#define AYLP_DEVICES_UDP_SOURCE_H_

#include <netinet/in.h>
#include "anyloop.h"

struct aylp_udp_source_data {
	int sock;
	struct sockaddr_in sa;
	// last received (or zero-initialized) command vector
	gsl_vector *vec;
	// units to output
	aylp_units units;
};

// initialize udp_source device
int udp_source_init(struct aylp_device *self);

// process udp_source device once per loop
int udp_source_proc(struct aylp_device *self, struct aylp_state *state);

// close udp_source device when loop exits
int udp_source_fini(struct aylp_device *self);

#endif
