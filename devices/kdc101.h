#ifndef AYLP_DEVICES_KDC101_H_
#define AYLP_DEVICES_KDC101_H_

#include "anyloop.h"

// Private data for the kdc101 device.
struct aylp_kdc101_data {
	int fd;			// FTDI serial port file descriptor
	int index;		// which pipeline vector element to command
	double scale;		// encoder counts per pipeline unit
	int32_t center;		// encoder counts at zero command
};

int kdc101_init(struct aylp_device *self);
int kdc101_proc(struct aylp_device *self, struct aylp_state *state);
int kdc101_fini(struct aylp_device *self);

#endif
