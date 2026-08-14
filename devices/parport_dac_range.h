#ifndef AYLP_DEVICES_PARPORT_DAC_RANGE_H_
#define AYLP_DEVICES_PARPORT_DAC_RANGE_H_

#include <stdbool.h>
#include <stddef.h>

struct dac_range {
	const char *name;
	int code;		// DACRANGE nibble
	double vmin, vmax;
};

extern const struct dac_range dac_ranges[];
extern const size_t dac_range_count;

const struct dac_range *dac_find_range(const char *name);
int dac_find_range_index(const char *name);
bool dac_range_contains(int outer, int inner);
double dac_range_span(int index);

// Return the narrowest supported range inside [base, ceiling] that contains
// volts. If volts lies outside the ceiling, return the ceiling so conversion
// clips there rather than silently selecting a range the supply cannot drive.
int dac_pick_range(int base, int ceiling, double volts);

// Apply the runtime widening/narrowing policy without touching hardware.
// Returns true and stores the requested new range in `target` when a range
// change is due. Widening is immediate; narrowing requires `dwell`
// consecutive samples inside `frac` of the narrower candidate.
bool dac_range_update(int base, int ceiling, int current, double volts,
	double frac, long dwell, long *count, int *target);

#endif
