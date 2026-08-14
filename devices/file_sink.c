#include <stdio.h>
#include <sys/stat.h>
#include <json-c/json.h>
#include <gsl/gsl_matrix.h>
#include "anyloop.h"
#include "block.h"
#include "logging.h"
#include "file_sink.h"
#include "xalloc.h"

/** Highest sequence number tried before giving up. */
#define FILE_SINK_MAX_SEQ 10000

/** Move an existing recording aside so this run cannot overwrite it.
*
* This device used to open with "wb", which truncates: every start silently
* destroyed the previous run's data, and there is no way to notice until you go
* looking for a run that is no longer there.
*
* The old file is RENAMED to `base.NNNN.ext` (lowest free NNNN) and the
* configured name is then created fresh, rather than this run being diverted to
* a new name. That ordering matters: the configured filename always holds the
* current run, so everything that reads a recording back by the name in the
* config -- contrib/calibration-scripts/tools/find_roi.py, the analysis scripts -- keeps working and keeps
* seeing the newest data, while every earlier run survives beside it.
*
* An existing but empty file is left alone and simply reused; rotating those
* would litter the directory with nothing every time a run failed at startup. */
static void rotate_existing(const char *fn)
{
	struct stat st;
	if (stat(fn, &st) != 0 || !S_ISREG(st.st_mode) || st.st_size == 0)
		return;

	// split off the extension so the sequence lands before it, keeping the
	// suffix meaningful to anything that dispatches on it
	const char *dot = strrchr(fn, '.');
	const char *slash = strrchr(fn, '/');
	if (dot && slash && dot < slash)
		dot = 0;	// the dot was in a directory name, not the file
	size_t stem = dot ? (size_t)(dot - fn) : strlen(fn);
	const char *ext = dot ? dot : "";

	size_t len = strlen(fn) + 16;
	char *cand = xmalloc(len);
	for (unsigned i = 1; i < FILE_SINK_MAX_SEQ; i++) {
		snprintf(cand, len, "%.*s.%04u%s", (int)stem, fn, i, ext);
		if (stat(cand, &st) == 0)
			continue;	// taken
		if (rename(fn, cand) == 0) {
			log_info("file_sink: kept the previous %s as %s", fn,
				cand
			);
		} else {
			log_warn("file_sink: could not move %s aside (%s); it "
				"is about to be overwritten", fn,
				strerror(errno)
			);
		}
		xfree(cand);
		return;
	}
	log_warn("file_sink: %d sequence slots for %s are all taken; the "
		"existing file is about to be overwritten", FILE_SINK_MAX_SEQ,
		fn
	);
	xfree(cand);
}


int file_sink_init(struct aylp_device *self)
{
	self->proc = &file_sink_proc;
	self->fini = &file_sink_fini;
	self->device_data = xcalloc(1, sizeof(struct aylp_file_sink_data));
	struct aylp_file_sink_data *data = self->device_data;
	if (!self->params) {
		log_error("No params object found.");
		return -1;
	}
	const char *fn = 0;
	json_object_object_foreach(self->params, key, val) {
		// parse parameters
		if (key[0] == '_') {
			// keys starting with _ are comments
		} else if (!strcmp(key, "filename")) {
			fn = json_object_get_string(val);
			log_trace("filename = %s", fn);
		} else if (!strcmp(key, "flush")) {
			data->flush = json_object_get_boolean(val);
			log_trace("flush = %d", data->flush);
		} else {
			log_warn("Unknown parameter \"%s\"", key);
		}
	}
	if (!fn) {
		log_error("You must provide the filename parameter.");
		return -1;
	}
	rotate_existing(fn);
	data->fp = fopen(fn, "wb");
	if (!data->fp) {
		log_error("Couldn't open file: %s", strerror(errno));
		return -1;
	}
	// set types and units
	self->type_in = AYLP_T_ANY;
	self->units_in = AYLP_U_ANY;
	self->type_out = AYLP_T_UNCHANGED;
	self->units_out = AYLP_U_UNCHANGED;
	return 0;
}


int file_sink_proc(struct aylp_device *self, struct aylp_state *state)
{
	struct aylp_file_sink_data *data = self->device_data;
	// write the header
	size_t n;
	n = fwrite(&state->header, 1, sizeof(struct aylp_header), data->fp);
	// write the data
	int needs_free = get_contiguous_bytes(&data->bytes, state);
	if (needs_free < 0) return needs_free;
	n += fwrite(data->bytes.data, 1, data->bytes.size, data->fp);
	size_t n_expect = sizeof(struct aylp_header) + data->bytes.size;
	if (n < n_expect) {
		log_error("Short write: %zu of %zu", n, n_expect);
		return -1;
	}
	if (needs_free) {
		xfree(data->bytes.data);
	}
	if (data->flush) {
		fflush(data->fp);
	}
	return 0;
}


int file_sink_fini(struct aylp_device *self)
{
	struct aylp_file_sink_data *data = self->device_data;
	fflush(data->fp);
	fclose(data->fp);
	xfree(data);
	return 0;
}

