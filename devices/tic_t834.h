#ifndef AYLP_DEVICES_TIC_T834_H_
#define AYLP_DEVICES_TIC_T834_H_

#include <libusb-1.0/libusb.h>
#include "anyloop.h"

// Private data for the tic_t834 device.
struct aylp_tic_t834_data {
	libusb_context *ctx;
	libusb_device_handle *devh;
	int index;		// which pipeline vector element to command
	double scale;		// steps per pipeline unit
	int32_t center;		// step position at zero command
};

int tic_t834_init(struct aylp_device *self);
int tic_t834_proc(struct aylp_device *self, struct aylp_state *state);
int tic_t834_fini(struct aylp_device *self);

#endif
