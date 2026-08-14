// Pololu Tic T834 stepper controller -- anyloop device.
// Drives one axis of the slow steering mirror via USB control transfers
// (libusb-1.0, no pololu-tic library required).
//
// The Tic must be configured for USB control (the default).  If the command
// timeout safety feature is enabled in the Tic settings, anyloop's loop rate
// must exceed 1/timeout; otherwise disable the timeout in the Tic settings.
//
// Params:
//   serial_number  -- (optional) USB serial number string to select a specific
//                     Tic when multiple are connected; first found if omitted
//   index          -- pipeline vector element to command (default 1, i.e. az/x)
//   scale          -- steps per pipeline unit
//   center         -- step position at zero command (parked / centre position)

#include <string.h>
#include <libusb-1.0/libusb.h>

#include "anyloop.h"
#include "tic_t834.h"
#include "logging.h"
#include "xalloc.h"

// ── Pololu Tic USB protocol ───────────────────────────────────────────────────

#define TIC_VID               0x1FFBu
#define TIC_T834_PID          0x00B5u

// bRequest command codes for USB vendor control transfers.
#define TIC_CMD_SET_TARGET_POSITION  0xE0u
#define TIC_CMD_ENERGIZE             0x85u
#define TIC_CMD_EXIT_SAFE_START      0x83u
#define TIC_CMD_DEENERGIZE           0x86u

// bmRequestType: vendor class, host-to-device, device recipient.
#define TIC_BMREQTYPE  0x40u

#define TIC_TIMEOUT_INIT_MS  2000
#define TIC_TIMEOUT_PROC_MS   200

// Send a command with a 32-bit signed integer parameter.
static int tic_cmd32(libusb_device_handle *devh, uint8_t cmd,
                     int32_t val, unsigned int timeout_ms)
{
	uint32_t u = (uint32_t)val;
	return libusb_control_transfer(devh,
		TIC_BMREQTYPE, cmd,
		(uint16_t)(u & 0xFFFFu),
		(uint16_t)(u >> 16),
		NULL, 0, timeout_ms);
}

// Send a command with no parameter.
static int tic_cmd0(libusb_device_handle *devh, uint8_t cmd,
                    unsigned int timeout_ms)
{
	return libusb_control_transfer(devh,
		TIC_BMREQTYPE, cmd,
		0, 0, NULL, 0, timeout_ms);
}

// Open the first Tic T834 whose serial number matches `serial` (NULL = any).
static libusb_device_handle *tic_open(libusb_context *ctx, const char *serial)
{
	libusb_device **list;
	ssize_t n = libusb_get_device_list(ctx, &list);
	if (n < 0) {
		log_error("tic_t834: get_device_list: %s",
			libusb_strerror((int)n));
		return NULL;
	}

	libusb_device_handle *result = NULL;
	for (ssize_t i = 0; i < n && !result; i++) {
		struct libusb_device_descriptor desc;
		if (libusb_get_device_descriptor(list[i], &desc)) continue;
		if (desc.idVendor != TIC_VID || desc.idProduct != TIC_T834_PID) continue;

		libusb_device_handle *h;
		if (libusb_open(list[i], &h)) continue;

		if (serial) {
			char sn[64] = {0};
			libusb_get_string_descriptor_ascii(h, desc.iSerialNumber,
				(unsigned char *)sn, (int)sizeof sn);
			if (strcmp(sn, serial)) {
				libusb_close(h);
				continue;
			}
		}
		result = h;
	}

	libusb_free_device_list(list, 1);
	return result;
}

// ── anyloop interface ─────────────────────────────────────────────────────────

int tic_t834_init(struct aylp_device *self)
{
	self->proc = &tic_t834_proc;
	self->fini = &tic_t834_fini;
	struct aylp_tic_t834_data *data = xcalloc(1, sizeof *data);
	self->device_data = data;

	// defaults
	data->index  = 1;
	data->scale  = 1000.0;
	data->center = 0;
	const char *serial = NULL;

	if (!self->params) {
		log_error("tic_t834: no params object");
		return -1;
	}
	json_object_object_foreach(self->params, key, val) {
		if (key[0] == '_') {
		} else if (!strcmp(key, "serial_number")) {
			serial = json_object_get_string(val);
			log_trace("tic_t834: serial_number = %s", serial);
		} else if (!strcmp(key, "index")) {
			data->index = (int)json_object_get_int64(val);
			log_trace("tic_t834: index = %d", data->index);
		} else if (!strcmp(key, "scale")) {
			data->scale = json_object_get_double(val);
			log_trace("tic_t834: scale = %g", data->scale);
		} else if (!strcmp(key, "center")) {
			data->center = (int32_t)json_object_get_int64(val);
			log_trace("tic_t834: center = %d", (int)data->center);
		} else {
			log_warn("tic_t834: unknown param \"%s\"", key);
		}
	}

	int r;
	if ((r = libusb_init(&data->ctx))) {
		log_error("tic_t834: libusb_init: %s", libusb_strerror(r));
		return -1;
	}

	data->devh = tic_open(data->ctx, serial);
	if (!data->devh) {
		log_error("tic_t834: no Tic T834 found%s%s",
			serial ? " with serial_number=" : "",
			serial ? serial : "");
		return -1;
	}

	// Detach any kernel driver (e.g. cdc_acm) and claim the vendor interface.
	libusb_set_auto_detach_kernel_driver(data->devh, 1);
	if ((r = libusb_claim_interface(data->devh, 0))) {
		log_error("tic_t834: claim interface 0: %s", libusb_strerror(r));
		return -1;
	}

	// Energize and clear safe-start so the motor accepts position commands.
	if ((r = tic_cmd0(data->devh, TIC_CMD_ENERGIZE, TIC_TIMEOUT_INIT_MS)) < 0) {
		log_error("tic_t834: energize: %s", libusb_strerror(r));
		return -1;
	}
	if ((r = tic_cmd0(data->devh, TIC_CMD_EXIT_SAFE_START,
	                  TIC_TIMEOUT_INIT_MS)) < 0) {
		log_error("tic_t834: exit safe start: %s", libusb_strerror(r));
		return -1;
	}

	self->type_in   = AYLP_T_VECTOR;
	self->units_in  = AYLP_U_ANY;
	self->type_out  = AYLP_T_UNCHANGED;
	self->units_out = AYLP_U_UNCHANGED;

	log_info("tic_t834: open%s%s  index=%d  scale=%g  center=%d",
		serial ? " sn=" : "", serial ? serial : "",
		data->index, data->scale, (int)data->center);
	return 0;
}

int tic_t834_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_tic_t834_data *data = self->device_data;

	if (UNLIKELY((size_t)data->index >= state->vector->size)) {
		log_error("tic_t834: index %d out of range (vector size %zu)",
			data->index, state->vector->size);
		return -1;
	}
	double cmd = gsl_vector_get(state->vector, data->index);
	int32_t pos = data->center + (int32_t)(cmd * data->scale);

	int r = tic_cmd32(data->devh, TIC_CMD_SET_TARGET_POSITION, pos,
	                  TIC_TIMEOUT_PROC_MS);
	if (UNLIKELY(r < 0))
		log_warn("tic_t834: set_target_position: %s", libusb_strerror(r));

	log_trace("tic_t834: cmd=%g → %d steps", cmd, (int)pos);
	return 0;
}

int tic_t834_fini(struct aylp_device *self)
{
	struct aylp_tic_t834_data *data = self->device_data;
	if (data->devh) {
		tic_cmd0(data->devh, TIC_CMD_DEENERGIZE, TIC_TIMEOUT_INIT_MS);
		libusb_release_interface(data->devh, 0);
		libusb_close(data->devh);
	}
	if (data->ctx) libusb_exit(data->ctx);
	xfree(data);
	return 0;
}
