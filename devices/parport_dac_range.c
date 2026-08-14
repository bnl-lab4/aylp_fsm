#include <string.h>

#include "parport_dac_range.h"

// DAC81404 datasheet SLASEH2A, table 8-16.
const struct dac_range dac_ranges[] = {
	{ "0-5",   0x0,   0.0,  5.0 },
	{ "0-10",  0x1,   0.0, 10.0 },
	{ "0-20",  0x2,   0.0, 20.0 },
	{ "0-40",  0x3,   0.0, 40.0 },
	{ "+-5",   0x5,  -5.0,  5.0 },
	{ "+-10",  0x6, -10.0, 10.0 },
	{ "+-20",  0x7, -20.0, 20.0 },
	{ "0-6",   0x8,   0.0,  6.0 },
	{ "0-12",  0x9,   0.0, 12.0 },
	{ "0-24",  0xA,   0.0, 24.0 },
	{ "+-6",   0xD,  -6.0,  6.0 },
	{ "+-12",  0xE, -12.0, 12.0 },
};

const size_t dac_range_count = sizeof dac_ranges / sizeof dac_ranges[0];

const struct dac_range *dac_find_range(const char *name)
{
	// Accept ASCII spellings of a bipolar range, e.g. "+-10" or "-10".
	if (name[0] == '+' && name[1] == '-') name += 2;
	else if (name[0] == '-') name += 1;
	else goto unipolar;
	for (size_t i = 0; i < dac_range_count; i++)
		if (dac_ranges[i].vmin < 0.0
		&& !strcmp(dac_ranges[i].name + 2, name))
			return &dac_ranges[i];
	return NULL;
unipolar:
	for (size_t i = 0; i < dac_range_count; i++)
		if (!strcmp(dac_ranges[i].name, name)) return &dac_ranges[i];
	return NULL;
}

int dac_find_range_index(const char *name)
{
	const struct dac_range *range = dac_find_range(name);
	if (!range) return -1;
	return (int)(range - dac_ranges);
}

bool dac_range_contains(int outer, int inner)
{
	return dac_ranges[outer].vmin <= dac_ranges[inner].vmin
		&& dac_ranges[outer].vmax >= dac_ranges[inner].vmax;
}

double dac_range_span(int index)
{
	return dac_ranges[index].vmax - dac_ranges[index].vmin;
}

int dac_pick_range(int base, int ceiling, double volts)
{
	int best = -1;
	for (size_t r = 0; r < dac_range_count; r++) {
		if (!dac_range_contains((int)r, base)) continue;
		if (!dac_range_contains(ceiling, (int)r)) continue;
		if (volts < dac_ranges[r].vmin || volts > dac_ranges[r].vmax)
			continue;
		if (best < 0 || dac_range_span((int)r) < dac_range_span(best))
			best = (int)r;
	}
	return best < 0 ? ceiling : best;
}

bool dac_range_update(int base, int ceiling, int current, double volts,
	double frac, long dwell, long *count, int *target)
{
	int want = dac_pick_range(base, ceiling, volts);
	if (volts < dac_ranges[current].vmin
	|| volts > dac_ranges[current].vmax) {
		*count = 0;
		if (want == current) return false;
		*target = want;
		return true;
	}
	if (want == current || dac_range_span(want) >= dac_range_span(current)) {
		*count = 0;
		return false;
	}
	if (volts < dac_ranges[want].vmin * frac
	|| volts > dac_ranges[want].vmax * frac) {
		*count = 0;
		return false;
	}
	if (++*count < dwell) return false;
	*target = want;
	return true;
}
