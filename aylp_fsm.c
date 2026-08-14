// Entry point under this repo's own name.
//
// A host derives the symbol it looks for from the file name it was given, so
// `file:/usr/local/lib/anyloop/aylp_fsm.so` makes it dlsym("aylp_fsm_init").
// That is how every other aylp_* plug-in in this lab is addressed, and it is
// the same device as the `fsm.so` symlink beside it -- see README.md.

#include "anyloop.h"
#include "devices/fsm.h"

int aylp_fsm_init(struct aylp_device *self)
{
	return fsm_init(self);
}
