#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "devices/parport_dac_range.h"

static void expect(const char *base_name, const char *max_name, double volts,
	const char *want_name)
{
	int base = dac_find_range_index(base_name);
	int ceiling = dac_find_range_index(max_name);
	assert(base >= 0 && ceiling >= 0);
	int got = dac_pick_range(base, ceiling, volts);
	if (strcmp(dac_ranges[got].name, want_name)) {
		fprintf(stderr, "%s..%s at %.9g V: got %s, expected %s\n",
			base_name, max_name, volts, dac_ranges[got].name,
			want_name);
		assert(0);
	}
}

int main(void)
{
	// The steering ladder: retain +-5 through its boundary, use the finer
	// +-6 intermediate range for modest excursions, and reach +-10 whenever
	// the command exceeds +-6 V.
	expect("+-5", "+-10",  0.0,    "+-5");
	expect("+-5", "+-10",  5.0,    "+-5");
	expect("+-5", "+-10",  5.0001, "+-6");
	expect("+-5", "+-10",  6.0,    "+-6");
	expect("+-5", "+-10",  6.0001, "+-10");
	expect("+-5", "+-10", 10.0,    "+-10");
	expect("+-5", "+-10", -5.0,    "+-5");
	expect("+-5", "+-10", -5.0001, "+-6");
	expect("+-5", "+-10", -6.0,    "+-6");
	expect("+-5", "+-10", -6.0001, "+-10");
	expect("+-5", "+-10", -10.0,   "+-10");
	// Outside the configured ceiling still selects +-10, where conversion
	// clamps; it must never silently select +-12 or +-20 on +-12 V rails.
	expect("+-5", "+-10", 11.0, "+-10");
	expect("+-5", "+-10", -11.0, "+-10");

	// A unipolar base retains the finer positive rungs but crosses to bipolar
	// immediately when a negative output is requested.
	expect("0-5", "+-10",  5.5, "0-6");
	expect("0-5", "+-10",  7.0, "0-10");
	expect("0-5", "+-10", -0.1, "+-5");
	expect("0-5", "+-10", -5.5, "+-6");
	expect("0-5", "+-10", -7.0, "+-10");

	assert(dac_find_range_index("-10") ==
		dac_find_range_index("+-10"));
	assert(dac_find_range_index("not-a-range") < 0);

	int base = dac_find_range_index("+-5");
	int ceiling = dac_find_range_index("+-10");
	int current = base, target = -1;
	long count = 0;
	assert(dac_range_update(base, ceiling, current, 5.01, 0.8, 3,
		&count, &target));
	assert(!strcmp(dac_ranges[target].name, "+-6"));
	current = target;
	assert(dac_range_update(base, ceiling, current, 6.01, 0.8, 3,
		&count, &target));
	assert(!strcmp(dac_ranges[target].name, "+-10"));
	current = target;
	// A sample outside the narrower range's hysteresis margin resets dwell.
	assert(!dac_range_update(base, ceiling, current, 3.0, 0.8, 3,
		&count, &target));
	assert(count == 1);
	assert(!dac_range_update(base, ceiling, current, 4.1, 0.8, 3,
		&count, &target));
	assert(count == 0);
	// Exactly three consecutive quiet samples narrow directly to the base.
	for (int i = 0; i < 2; i++)
		assert(!dac_range_update(base, ceiling, current, 3.0, 0.8, 3,
			&count, &target));
	assert(dac_range_update(base, ceiling, current, 3.0, 0.8, 3,
		&count, &target));
	assert(target == base);
	puts("parport_dac range-selection tests passed");
	return 0;
}
