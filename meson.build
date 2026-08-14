project(
	'aylp_fsm',
	'c',
	default_options: ['c_std=gnu17', 'warning_level=3', 'optimization=3']
)

c = meson.get_compiler('c')
deps = [
	dependency('json-c'),
	dependency('gsl'),
	dependency('threads'),
	dependency('libusb-1.0'),
	c.find_library('dl')
]

incdir = include_directories('libaylp')

# libaylp, carried rather than symlinked from an anyloop checkout the way the
# other aylp_* plug-ins do it: these devices need a libaylp NEWER than the one
# upstream ships -- timing.c, and anyloop.h's measurement-validity flags -- so
# an upstream symlink would not compile. Only timing.c ends up in the plug-in;
# the rest is here so the unit tests can link without an anyloop checkout, and
# comes from the host at run time. See README.md.
libaylp_source_files = [
	'libaylp/logging.c',
	'libaylp/block.c',
	'libaylp/pretty.c',
	'libaylp/types.c',
	'libaylp/thread_pool.c',
	'libaylp/timing.c',
	'libaylp/xalloc.c',
]

device_source_files = [
	'devices/center_of_mass.c',
	'devices/wfs_com.c',
	'devices/fit_com.c',
	'devices/kdc101.c',
	'devices/tic_t834.c',
	'devices/clamp.c',
	'devices/biquad_filter.c',
	'devices/kalman_filter.c',
	'devices/device.c',
	'devices/delay.c',
	'devices/file_sink.c',
	'devices/logger.c',
	'devices/matmul.c',
	'devices/pid.c',
	'devices/fsp.c',
	'devices/fsm.c',
	'devices/poke.c',
	'devices/remove_piston.c',
	'devices/stop_after_count.c',
	'devices/test_source.c',
	'devices/udp_sink.c',
	'devices/udp_source.c',
	'devices/vonkarman_stream.c',
	'devices/piplate_bridge.c',
	'devices/parport_dac.c',
	'devices/parport_dac_range.c',
	'devices/bode_plot.c',
	'devices/latency_test.c',
	'devices/gain_test.c',
	'devices/prbs_test.c',
	'devices/attenuation_test.c',
	'devices/spiral_search.c',
]

# ---------------------------------------------------------------------------
# The plug-in. One library, loadable by a stock anyloop; see README.md for why
# libaylp is left undefined here and what the two linker flags are for.
plugin_map = meson.current_source_dir() / 'aylp_fsm.map'
aylp_fsm = shared_library('aylp_fsm',
	device_source_files + ['libaylp/timing.c', 'aylp_fsm.c'],
	name_prefix: '',
	install: true,
	install_dir: get_option('libdir') / 'anyloop',
	dependencies: deps,
	include_directories: [incdir, include_directories('.')],
	link_args: [
		'-Wl,--version-script=' + plugin_map,
		'-Wl,-Bsymbolic-functions',
	],
	link_depends: plugin_map,
	override_options: 'b_lundef=false'
)

# One name per device, all pointing at the one library above, because the host
# looks for an init symbol named after the file it was told to open.
plugin_devices = [
	'center_of_mass', 'wfs_com', 'fit_com', 'kdc101', 'tic_t834', 'clamp',
	'biquad_filter', 'kalman_filter', 'delay', 'file_sink', 'logger',
	'matmul', 'pid', 'fsp', 'fsm', 'poke', 'remove_piston',
	'stop_after_count', 'test_source', 'udp_sink', 'udp_source',
	'vonkarman_stream', 'piplate_bridge', 'parport_dac', 'bode_plot',
	'latency_test', 'gain_test', 'prbs_test', 'attenuation_test',
	'spiral_search',
]
make_link = find_program('make_link.sh')
foreach device : plugin_devices
	# in the build tree, so the plug-in can be exercised without installing
	custom_target(device + '_link',
		input: aylp_fsm,
		output: device + '.so',
		command: [make_link, '@INPUT@', '@OUTPUT@'],
		build_by_default: true,
		install: false
	)
	# ...and on install. This has to be install_symlink() rather than
	# installing the custom_target above: meson copies a custom_target's
	# output by dereferencing it, which would put thirty 1.3 MB copies of
	# the same library in libdir instead of thirty symlinks.
	install_symlink(device + '.so',
		pointing_to: 'aylp_fsm.so',
		install_dir: get_option('libdir') / 'anyloop'
	)
endforeach

# ---------------------------------------------------------------------------
# Unit tests. These link the device sources directly rather than going through
# the plug-in, so they need no anyloop, no camera and no config file --
# `meson test -C build` is safe to run on any machine.
test_center_of_mass = executable('test_center_of_mass',
	[
		'devices/center_of_mass_test.c',
		'devices/center_of_mass.c',
		'libaylp/logging.c',
		'libaylp/thread_pool.c',
		'libaylp/types.c',
		'libaylp/xalloc.c'
	],
	install: false,
	dependencies: deps,
	include_directories: [incdir, include_directories('devices')]
)
test('center_of_mass', test_center_of_mass)

test_pid_park = executable('test_pid_park',
	[
		'devices/pid_park_test.c',
		'devices/pid.c',
		'libaylp/logging.c',
		'libaylp/types.c',
		'libaylp/xalloc.c'
	],
	install: false,
	dependencies: deps,
	include_directories: [incdir, include_directories('devices')]
)
test('pid_park', test_pid_park)

test_gain_test = executable('test_gain_test',
	[
		'devices/gain_test_test.c',
		'devices/gain_test.c',
		'libaylp/logging.c',
		'libaylp/timing.c',
		'libaylp/xalloc.c'
	],
	install: false,
	dependencies: deps,
	include_directories: [incdir, include_directories('devices')]
)
test('gain_test', test_gain_test)

test_prbs_test = executable('test_prbs_test',
	[
		'devices/prbs_test_test.c',
		'devices/prbs_test.c',
		'libaylp/logging.c',
		'libaylp/xalloc.c'
	],
	install: false,
	dependencies: deps,
	include_directories: [incdir, include_directories('devices')]
)
test('prbs_test', test_prbs_test)

test_fsp_dark = executable('test_fsp_dark',
	[
		'devices/fsp_dark_test.c',
		'devices/fsp.c',
		'libaylp/logging.c',
		'libaylp/timing.c',
		'libaylp/types.c',
		'libaylp/xalloc.c'
	],
	install: false,
	dependencies: deps,
	include_directories: [incdir, include_directories('devices')]
)
test('fsp_dark', test_fsp_dark, timeout: 120)

# fsm composes real substages through init_device(), so it links the whole
# device registry. It still needs no hardware: the substages it drives are
# test_source, clamp and file_sink.
test_fsm = executable('test_fsm',
	['devices/fsm_test.c'] + libaylp_source_files + device_source_files,
	install: false,
	dependencies: deps,
	include_directories: [incdir, include_directories('devices'),
		include_directories('.')]
)
test('fsm', test_fsm)

test_fit_com = executable('test_fit_com',
	['devices/fit_com_test.c', 'devices/fit_com.c', 'libaylp/logging.c',
		'libaylp/xalloc.c'],
	install: false,
	dependencies: deps,
	include_directories: [incdir, include_directories('devices')])
test('fit_com', test_fit_com)

test_parport_dac_range = executable('test_parport_dac_range',
	['devices/parport_dac_range_test.c', 'devices/parport_dac_range.c'],
	install: false,
	include_directories: include_directories('.'))
test('parport_dac_range', test_parport_dac_range)
