// Type and units name <-> enum conversion.
//
// These live apart from anyloop.c purely so that linking a device does not drag
// in main(): every device that takes a `type` or `units` param calls into here,
// and a unit test that links one device would otherwise pull the whole program
// entry point in with it.

#include <string.h>

#include "anyloop.h"
#include "logging.h"


aylp_type aylp_type_from_string(const char *type_name)
{
	#define AYLP_TYPE_FROM_STRING_MATCH_TYPE(TYPE, type) \
	if (!strcasecmp(type_name, #type)) return TYPE;
	FOR_AYLP_TYPES(AYLP_TYPE_FROM_STRING_MATCH_TYPE)

	log_error("Couldn't parse type: %s", type_name);
	return AYLP_T_NONE;
}

const char *aylp_type_to_string(aylp_type type)
{
	switch (type) {
	#define AYLP_TYPE_TO_STRING_MATCH_TYPE(TYPE, type) \
	case TYPE: \
		return #type; \
		break;
	FOR_AYLP_TYPES(AYLP_TYPE_TO_STRING_MATCH_TYPE)
	default:
		log_error("Unknown type 0x%hhX", type);
		return "NONE";
	}
}

aylp_units aylp_units_from_string(const char *units_name)
{
	#define AYLP_UNITS_FROM_STRING_MATCH_UNITS(UNITS, units) \
	if (!strcasecmp(units_name, #units)) return UNITS;
	FOR_AYLP_UNITS(AYLP_UNITS_FROM_STRING_MATCH_UNITS)

	log_error("Couldn't parse units: %s", units_name);
	return AYLP_U_NONE;
}

const char *aylp_units_to_string(aylp_units units)
{
	switch (units) {
	#define AYLP_UNITS_TO_STRING_MATCH_UNITS(UNITS, units) \
	case UNITS: \
		return #units; \
		break;
	FOR_AYLP_UNITS(AYLP_UNITS_TO_STRING_MATCH_UNITS)
	default:
		log_error("Unknown units 0x%hhX", units);
		return "NONE";
	}
}
