anyloop:parport_dac
===================

Types and units: `[T_VECTOR, U_MINMAX|U_V] -> [T_UNCHANGED, U_UNCHANGED]`.

Drives a TI DAC81404 (quad 16-bit, high-voltage output) through a PCIe
parallel-port card — the fast replacement for the `piplate_bridge` /
DAQC2 serial DAC path, whose ~0.5–0.7 ms per `setDAC` was a large share of the
loop's transport delay.

Like `piplate_bridge`, one stage can drive several channels: pass arrays for
`channel`/`index`/`scale`/`offset`/`range`/`start_delay` (any scalar is
broadcast to every channel), and each channel's output voltage is
`offset + cmd*scale`, where `cmd` is the pipeline vector element named by
`index`.

**Status: working — the DAC drives real voltage.** As of 2026-07-29 the `spi`
link produces a measured output at the EVM's `DAC_VOUT_A` terminal, tracking
the commanded value (`contrib/calibration-scripts/configurations/conf_parport_dac_2v.json`: channel 0, range
`+-10`, `scale: 0.0`, so the output parks at `offset` — driven at 2 V, then
3 V). That is the end-to-end confirmation the "DAC power-up" note below was
missing: every earlier bench session had measured only SPI *bus* timing
(frames/s), which a device still in device-wide power-down would have produced
just as happily.

Hardware
--------

The transport is a StarTech PEX1P2 (ASIX AX99100) PCIe parallel card, seen on
this machine as PCI function `0000:05:00.2`. Its registers are available
through *memory* BARs, so a pin change is a plain store to mapped memory
(~100–300 ns, posted) rather than an `outb` (~1–2 µs, non-posted).

Two ways to get from the DB25 to the DAC, chosen with `link`:

  - `"spi"` (default): the DB25 data pins *are* the SPI bus — D0 = SCLK,
    D1 = SDIN, D2 = SYNC by default, GND on pins 18–25. SYNC (frame bracket)
    and LDAC (load strobe, unused here — see Update mode below) are both
    active-low pins, per the DAC81404's own datasheet naming (no "Z"
    suffix — see the pinout table's footnote). The 24-bit frames are
    bit-banged here: two stores per bit plus the frame
    brackets. **Measured
    2026-07-27: 156 k frames/s = 6.4 µs per channel update** (~123 ns per
    store, effective SCLK ≈ 3.75 MHz), so a two-channel update costs ~12.8 µs
    against the BRIDGEplate's ~0.5–1.4 ms. That is far below every SCLK
    ceiling the part specifies — 50 MHz for writes at IOVDD 2.7–5.5 V
    (SLASEH2A 7.7), but only 35 MHz for reads with `FSDO=1` and 20 MHz with
    `FSDO=0` (7.10, 7.11), which is worth knowing because `verify` reads.
    Needs nothing but wires between the card and the BP-DAC81404EVM.
  - `"pico"`: the DB25 carries a byte-parallel handshake to an RP2040/RP2350
    whose PIO runs the SPI at up to 50 MHz. Four stores per sample (~0.6 µs):
    low byte, clock high, high byte, clock low, with the channel number on
    spare control lines. The Pico owns the DAC's init and LDAC timing. Needs
    the bridge firmware.

**Levels — measured 2026-07-27, and the two pin groups are NOT the same.**

  - **Data pins D0–D7 (DB25 2–9): 3.3 V.** Driven high they measure 3.3 V,
    unbuffered, idle pins 0 V. Safe to wire straight to the EVM (IOVDD 3.3 V,
    V_IH = 0.7 × IOVDD = 2.31 V, so ~1 V of margin) or to an RP2040/RP2350.
  - **Control pins (DB25 1, 14, 16, 17): 4.94 V.** These *are* buffered to
    IEEE-1284 5 V levels. **Do not connect them to the DAC or to a Pico
    without level shifting** — 5 V exceeds the DAC81404's absolute maximum at
    3.3 V IOVDD, and RP2040/RP2350 GPIOs are not 5 V tolerant. This is the
    constraint on the `pico` link, whose clock and channel-select lines live
    on control pins: it needs a divider, series resistors, or a 74LVC245 on
    those four lines even though the data byte is fine as-is. The same applies
    to putting a hardware LDAC line on a control pin.

Keep leads under ~20 cm and put 33–100 Ω in series with SCLK for edge quality.

The DB25 breakout's screw terminals are fitted with **plastic covers** — pull
one to land or meter a wire, then put it back. Worth keeping on: the 4.94 V
control pins sit immediately alongside the 3.3 V data pins that run straight
to the DAC, so a slipped probe or a stray strand bridging the two groups is
the one mistake at this connector that reaches the DAC's absolute maximum.

Pinout (measure against ground on pins 18–25):

| signal | D0 | D1 | D2 | D3 | D4 | D5 | D6 | D7 |
|---|---|---|---|---|---|---|---|---|
| DB25 pin | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |

DB25 pin 1 is `/STROBE`, a *control* line, not D0 — the control lines are
pins 1 (STROBE), 14 (AUTOFD), 16 (INIT) and 17 (SELECT), and all but INIT are
inverted between the register bit and the pin. This device takes control
values in pin polarity and folds the inversion in for you.

Wiring to the BP-DAC81404EVM
---------------------------

From SLAU825 Figure 16 (schematic, p. 19) — note the signals are split across
**both** BoosterPack receptacles, and J4's physical pin numbers run opposite to
the BoosterPack numbering, which makes the Figure 7 pinout diagram easy to
misread:

| DB25 | signal | connector | physical pin | BoosterPack pin | net |
|---|---|---|---|---|---|
| 2 (D0) | SCLK | **J1** | **13** | 7 | `DAC_SCLK` |
| 3 (D1) | SDIN | **J4** | **12** | 15 | `DAC_SDIN` |
| 4 (D2) | SYNC | **J4** | **18** | 12 | `DAC_SYNCZ`(1) |
| 18–25 | GND | **J4** | **2** | 20 | `DGND` |
| 10 (ACK) | SDO, optional | **J4** | **14** | 14 | `DAC_SDO` |

(1) The "signal" column here is the chip's own pin name, straight from the
DAC81404 datasheet (SLASEH2A) — confirmed 2026-07-28 against a local copy,
`dac81404.pdf`: the pin table literally reads "SDO SCLK SDIN SYNC LDAC GND
IOVDD CLR", no "Z" on SYNC or LDAC anywhere. The "net" column is the
BP-DAC81404EVM board's *own* schematic net name (SLAU825), which is the EVM
designer's active-low notation and doesn't have to match the chip's pin
name — it doesn't, here. Don't take `DAC_SYNCZ` as evidence the pin itself
is called SYNCZ; the datasheet is the authority for that.

*As built, 2026-07-29:* this table is the wiring actually in place — D0 →
SCLK, D1 → SDIN, D2 → SYNC (EVM net `DAC_SYNCZ`). That matches the device's
default `pin_sclk`/`pin_sdin`/`pin_sync` of 0/1/2, so no pin overrides belong
in the config.

Count physical pins: pin 1 is the square pad, then double-row order (1,2 across
the top, 3,4 next), odd down one column, even down the other. Check J4 pin 2
reads continuity to ground before powering anything. Because SCLK sits on the
opposite edge from SDIN/SYNC, take a ground at **each** connector — J4 pin 2
and J1 pin 4 (`DGND`, BoosterPack 22) — so SCLK's return doesn't have to cross
the board.

**LDAC, CLR and RST already sit at IOVDD — no wires needed, and none fitted.**
*As built, 2026-07-29:* nothing from the DB25 is connected to any of the three.
The EVM holds all of them at 3.3 V through its own **10 kΩ pull-ups**, and
**J2 is removed**, so LDAC is *not* pulled to ground. That is exactly the state
this section had been asking for, reached through the board's pull-up network
instead of through added wires. (In EVM net names these are `LDACZ`, `CLRZ` and
`RSTZ`; the chip's own pin names drop the "Z" — see footnote (1) above.)

Leave them that way. This device drives none of them, and for LDAC and CLR the
DAC81404's pin table (SLASEH2A Table 6-1) asks for precisely this —
*"Connect to IOVDD if unused"*. All three are active-low inputs with no
documented **internal** pull-up, so those external 10 kΩ are what makes
"unused" mean "inactive" rather than "floating next to a bit-banged bus". They
are load-bearing, not decoration. 10 kΩ is a weak hold, so keep any stub short
and don't leave a probe lead hanging off one. Each fails differently if it does
get dragged low, worst last.

  - **LDAC** (chip pin 13, EVM jumper J2 — **removed**, pin held high by its
    pull-up). Grounding it (J2 1–2, the EVM default, and what the earlier
    bring-up notes here called for) is harmless *today* only because every
    channel is in asynchronous mode, where LDAC is ignored outright (see
    Update mode below). It stops being harmless the moment any channel is put
    in synchronous mode — which `sync_update` now does on request, so this is
    a live constraint rather than a hypothetical: synchronous channels update
    "when the LDAC pin is low", so a grounded LDAC would quietly turn the mode
    back into asynchronous. Don't refit J2 as a "restore the default" step
    during later rework — off is the wanted state.
  - **CLR** (chip pin 16, EVM jumper J3). This one matters more, because a
    single low glitch is destructive to loop state: "a clear command forces
    all DAC channels to clear the contents of their buffer **and active**
    registers to the clear code regardless of their synchronization setting"
    (SLASEH2A 8.3.3.3) — zero code on a unipolar range, midscale on a bipolar
    one, on all four channels at once. Worse, with `skip_unchanged: true` this
    device would not notice: it re-sends a channel only when its computed code
    *changes*, so the mirror would sit at the clear code until the command
    happened to move by ≥ 1 LSB. J3 is off, which by itself would leave the
    pin floating rather than deasserted — the 10 kΩ pull-up is what actually
    holds it high.

  - **RST** (chip pin 32) — **also to IOVDD**, though for a slightly different
    reason: the pin table gives no "connect to IOVDD if unused" note for this
    one, so that instruction is *not* in the datasheet. It doesn't need to be.
    RST is "an active-low reset input. Logic low on this pin causes the device
    to issue a power-on-reset event", so holding it high is simply what running
    the part means. §11.1 confirms TI expects it tied off rather than driven:
    "the RST and FAULT signals are static lines; therefore these lines can lie
    on the analog side of the ground plane."

    Of the three, a glitch here is the worst outcome. `tRSTW` is 20 ns, so a
    very short transient is enough, and a POR "causes all registers to
    initialize to default values" (§8.3.5) — `DEV-PWDWN` set, reference off,
    all four channels powered down and clamped to ground, ranges back to
    0–5 V, `FSDO` back to 0, `CRC-EN` off. This device configures the DAC once
    in `init` and never re-checks (a state now easy to recognise: the output stops
moving and stays put, while frames keep going out and `verify` keeps passing —
`verify` only proves the bus is intact, not that the part is configured). From
that moment every frame it sends is
    silently ignored and no output moves again until anyloop is restarted.
    There is no reason to keep a hardware reset line: `reset_at_init` already
    issues a `SOFT-RESET` through `TRIGGER`, which §8.3.5 defines as the same
    POR event. Holding it high is all that's wanted, and the EVM's 10 kΩ
    pull-up already does it — while still leaving the pin free for a button or
    testpoint if manual recovery is ever needed.

The EVM's own jumper wiring is in SLAU825, which is not saved locally, so don't
take a jumper position here on faith — confirm with a meter that chip pins 13,
16 and 32 all read ≈ IOVDD (3.3 V) against ground. On the current build they do,
with no jumper fitted on J2 or J3 and no wire on any of the three. Re-check
after any rework that disturbs those jumpers or the pull-up network.

Other jumpers: J16 none = internal reference, which init powers up; J11 2–3 =
unipolar. **J19 is the one
to think about**: its 2–3 default draws DAC IOVDD from the LaunchPad connector,
so with no LaunchPad either feed 3.3 V into J1 physical pin 1 (`AEVM_3V3`) or
move J19 to 1–2 and feed terminal block J20. It must be 3.3 V to match the
card's data-pin levels.

Supplies land on terminal block J17: pin 1 DAC_VDD (0–5.5 V), pins 2/3 GND,
pin 4 DAC_AVDD (≤41.5 V, needs headroom over the output range), pin 5 DAC_AVSS
(ground it in unipolar mode). AVDD − AVSS must not exceed 42 V. Outputs are on
the four 5-pin terminal blocks silkscreened `DAC_VOUT_A`…`DAC_VOUT_D`, each
with force and sense.

The EVM already carries 33 Ω in series on SCLK/SDIN/SYNC at its own end (R18,
R27, R28). That damps ringing into the DAC input but is not source termination
for a long cable, so it does not replace resistors at the DB25 end.

**Privileges and kernel lockdown.** The `mmio` backend needs `CAP_SYS_RAWIO`
(run anyloop as root) and needs `parport_pc` off the device, which `unbind`
does by default. Rebind afterwards with

```
echo 0000:05:00.2 | sudo tee /sys/bus/pci/drivers/parport_pc/bind
```

It also needs kernel lockdown off. With Secure Boot **on**, lockdown is
`integrity` and mapping the BAR fails with `EPERM` even as root (`Lockdown:
direct PCI access is restricted` in `dmesg`); lockdown blocks `ioperm`/`outb`
the same way, leaving no user-space fast path at all. Secure Boot was disabled
on this machine on 2026-07-27, so `cat /sys/kernel/security/lockdown` now
reads `[none] integrity confidentiality` — the *bracketed* word is the active
mode, which is easy to misread.

The `ppdev` backend works under lockdown and without root (given `lp` group
access to `/dev/parport0`), but every edge becomes an ioctl. Measured on this
machine: **~4 µs per edge, ~200 µs per 24-bit frame, ~4.9 k frames/s**, so a
two-channel update costs ~0.4 ms — about what the BRIDGEplate cost. Use it to
check wiring and frames, not to close a loop.

**Update mode: asynchronous by default, synchronous on request.** Which one a
channel uses is set by its `DACx-SYNC-EN` bit in `SYNCCONFIG` (bit N = channel
N, SLASEH2A 8.6.7), and `sync_update` picks between them for every channel this
stage drives.

  - **Asynchronous** (bit clear, the default): "a DAC data register write
    results in an immediate update of the DAC active register and DAC output
    on a SYNC rising edge" (SLASEH2A 8.3.3.1.2). Each channel therefore moves
    as its own frame closes, so two channels written in the same iteration
    step **~6.4 µs apart** on the `spi` link.
  - **Synchronous** (bit set, `sync_update: true`): the `DACx` write only
    loads that channel's *buffer* register and nothing moves. The outputs
    update on a trigger, which §8.3.3.1.1 says is "generated either through
    the SOFT-LDAC bit or by the LDAC pin". Since the LDAC pin is not wired
    here, this device uses the bit: it writes every channel, then issues one
    `TRIGGER` write with **SOFT-LDAC** (bit 4, table 8-17 — *"set this bit to
    1 to synchronously load the DACs that have been set in synchronous mode
    in the SYNCCONFIG register"*), and all of them step together on that
    frame's SYNC rising edge.

**LDAC the pin still does nothing**, in either mode — this device never drives
it, and on this build it is not connected at all (see "LDAC, CLR and RST"
under Wiring). That is a precondition for `sync_update`, not a detail:
channels in synchronous mode "are updated simultaneously when the LDAC pin is
low", so refitting **J2 1–2** to ground LDAC would hold the transfer
permanently open and push each buffer write straight through to the output —
silently degrading synchronous mode back into asynchronous, with the
simultaneity gone and nothing to indicate it. J2 stays off.

`dac_configure()` writes `SYNCCONFIG` explicitly either way rather than
trusting the reset default (all-clear). With `reset_at_init: false` against a
device a previous process left in the *other* mode, that write is the only
thing standing between you and outputs that accept and acknowledge every
`DACx` write while never moving — a miserable thing to debug from the outside.

**Output update rate.** SLASEH2A 8.3.3.1 (and `tDACWAIT` in tables 7-6/7-7)
requires **at least 2.4 µs between DAC output updates**, in both update modes.
A bit-banged 24-bit frame takes ~6.4 µs, so the `spi` link clears this on its
own — but the `pico` link hands a sample over in ~0.6 µs, so two channels in
one iteration would have been issued about four times too fast. The device now
enforces the spacing itself rather than leaving it to be true by accident.

**Codes and volts.** SLASEH2A equations 1 and 2 give
`VOUT = VREFIO × GAIN × CODE / 2^N` (unipolar) and the same minus
`GAIN × VREFIO / 2` (bipolar) — the range is divided into **2^16 steps, not
2^16 − 1**, so code `0xFFFF` sits one LSB *below* the nominal top of the range
(9.99985 V on `0-10`, not 10 V). `range`'s `vmax` is that nominal top, i.e. a
voltage the part cannot quite reach; asking for it clamps to `0xFFFF`.

**DAC power-up.** The DAC81404 comes out of reset *device-wide* powered down
(`SPICONFIG.DEV-PWDWN`), separately from its reference and all four output
amplifiers also being powered down (outputs clamped to ground through 10 kΩ).
Found 2026-07-28 reading a local copy of the datasheet: nothing here used to
clear `SPICONFIG.DEV-PWDWN` at all, which is a different, higher-level gate
than the reference/channel power bits below it — so every earlier bench
session that only measured SPI *bus* timing (frames/s) rather than an actual
DAC output voltage could plausibly have been talking to a device that never
left device-wide power-down. With that bit now cleared in `init`, the part
does drive its outputs (2026-07-29, see Status above) — which is the evidence
that this was the missing piece, and the reason to measure a voltage rather
than a frame rate when judging whether this device works. In the `spi` link, init writes `SPICONFIG`
(device active, and `FSDO=1` for `verify`'s benefit — see below), `GENCONFIG`
(reference on), `SYNCCONFIG` (asynchronous update), `DACRANGE` (the configured
ranges) and `DACPWDWN`, then parks every channel at its `offset` before the
loop starts. `DACPWDWN` powers up **only the channels this stage drives**; the
rest keep their reset power-down bit and stay clamped to ground through the
internal 10 kΩ, rather than being released into a range that was never
configured for them. When `verify` is set, each of those init writes is
checked over SDO/ACK the same way the per-channel writes are
(`dac_write_checked()`), logging once — so a clean startup rules out the
*whole* digital chain up front, not just whatever channel write happens to run
first. `TRIGGER`'s soft-reset is deliberately excluded: its fields are
momentary actions (SLASEH2A table 8-17 types them all "W", reset `0000h`), so
reading `SOFT-RESET` back afterward isn't expected to still show the `1010b`
that was sent — a "mismatch" there would be a false alarm, not a real one.

**Which registers can be read back at all.** Half the DAC81404's map is
write-only, per the TYPE column of SLASEH2A table 8-7: `NOP`, `DACRANGE`,
`BRDCAST` and `DACA`–`DACD` are typed plain `W`, and `TRIGGER`'s every field
is `W` in table 8-17. Only `DEVICEID`, `STATUS`, `SPICONFIG`, `GENCONFIG`,
`BRDCONFIG`, `SYNCCONFIG` and `DACPWDWN` have contents that come back. This
matters because `verify` used to issue a read command for *every* register it
wrote, including `DACRANGE` and every per-channel `DACx` write — the two it
writes most — and then compare the reply against what was sent. There is no
reply to compare, so that reliably produced a bogus "write not confirmed"
warning on exactly the writes the feature exists to check. See `verify` below
for what it does instead.

**Register safety.** SLASEH2A 8.6: "All register addresses not listed
should be considered as reserved locations and the register contents
should not be modified" — on parts like this, reserved space is commonly
factory trim/calibration, so a stray write there can do real, permanent
damage. Every write in this file funnels through one function
(`spi_frame_rw()`), which checks the target address against an allow-list
built from Table 8-7 (`dac_addr_writable()`) before anything is sent, and
refuses (logging an error, no frame transmitted) if it isn't one of NOP,
SPICONFIG, GENCONFIG, BRDCONFIG, SYNCCONFIG, DACPWDWN, DACRANGE, TRIGGER,
BRDCAST, or DACA–D. This is enforced in code, not just true by inspection of
the current callers — the only runtime-variable register address this file
ever constructs is `DACx = DAC0 + channel`, and `channel` is itself bounds
checked to 0–3 at init, so it can't drift outside 0x10–0x13 either.

Parameters
----------

### Transport

  - `backend` (string) (optional)
    - `"mmio"` (default) or `"ppdev"`.
  - `pci` (string) (optional)
    - PCI slot of the parallel function (default `"0000:05:00.2"`). Find it
      with `lspci -d 125b: -nn`.
  - `bar` (integer) (optional)
    - Memory BAR index holding the parallel port's registers (default 2).
  - `data_offset`, `status_offset`, `ctrl_offset` (integer) (optional)
    - Register offsets inside `bar`, defaults `0x280`, `0x284`, `0x288`.
      The AX99100 does **not** expose the byte-packed ISA layout through its
      memory BAR: the block is dword-spaced at +0x280, mirrored again at
      +0x2c0, with the whole 0x400 window repeating through the 4 kB page.
      Mapped and verified on this card on 2026-07-27 by writing through the
      memory window and reading the I/O window at 0x3010, and vice versa.
  - `ecr_bar`, `ecr_offset`, `ecr_value` (integer) (optional)
    - Where the ECP ECR lives — defaults 2 and `0x2a8`, i.e. the *same* BAR as
      the other registers (BAR5 holds the AX99100's own configuration
      registers and does not mirror the ECR) — and what to write to it
      (default `0x14`: SPP/compatibility mode, interrupts off) so the data
      pins are plain outputs. Set `ecr_bar` to -1 to leave the ECR alone.
  - `unbind` (boolean) (optional)
    - Unbind `parport_pc` from `pci` before mapping (default true).
  - `probe` (boolean) (optional)
    - Read the data register back at init and check it holds what was written
      (default true). Catches a wrong `bar`/`data_offset`, which otherwise
      looks like a DAC that silently ignores every frame. Only the idle
      pattern is written, so an attached DAC sees no edges.
  - `port` (string) (optional)
    - ppdev node (default `"/dev/parport0"`).

### Link

  - `link` (string) (optional)
    - `"spi"` (default) or `"pico"`.
  - `sclk_bit`, `sdin_bit`, `sync_bit` (integer) (optional)
    - Data pins (0–7) carrying the SPI signals in the `spi` link
      (defaults 0, 1, 2 = DB25 pins 2, 3, 4).
  - `clk_ctrl_bit` (integer) (optional)
    - Control line clocking bytes into the bridge in the `pico` link:
      0 = STROBE (pin 1), 1 = AUTOFD (14), 2 = INIT (16), 3 = SELECT (17).
      Default 0.
  - `chan_ctrl_shift` (integer) (optional)
    - First of the two control lines carrying the channel number in the
      `pico` link (default 1, i.e. AUTOFD + INIT).
  - `delay_ns` (integer) (optional)
    - Extra dwell after each edge, in nanoseconds (default 0). The store
      latency alone (~123 ns) is an order of magnitude past the DAC's
      `tSDIS`/`tSDIH` of 10 ns at IOVDD 2.7–5.5 V, 20 ns below that
      (SLASEH2A 7.6, 7.7); raise this only for long or unterminated wiring.
  - `verify` (boolean) (optional, `spi` link only)
    - After every write — including the one-shot init sequence
      (`SPICONFIG`/`GENCONFIG`/`SYNCCONFIG`/`DACRANGE`/`DACPWDWN`), not just
      the per-channel writes — check over SDO (DB25 pin 10 → EVM `DAC_SDO`,
      only meaningful if that wire is connected) that it landed, and warn if
      it didn't (default false).

      *How much* can be checked depends on the register, because half the map
      is write-only (see "Which registers can be read back at all" above):

      | register | what `verify` does | frames |
      |---|---|---|
      | `SPICONFIG`, `GENCONFIG`, `SYNCCONFIG`, `DACPWDWN` (R/W) | reads the value back and compares it | 3 |
      | `DACRANGE`, `DACA`–`DACD` (W) | checks the R/W+address echo the next cycle carries (Table 8-3) | 2 |
      | any of the above, with `crc` on | compares the value: a CRC-mode echo carries the previous cycle's *data* too (Table 8-5) | 2–3 |

      The log says which of the two actually happened — `... write confirmed
      by SDO/ACK readback: ...` versus `... write acknowledged by SDO/ACK
      (R/W+address echo) ...` — so a quiet run can't be mistaken for more
      assurance than it gives. Each init write logs once, since it only
      happens once; the first time a given channel checks out it likewise
      logs once at `info` level, then stays quiet on success and only logs
      again on a problem, so a clean run isn't spammed once per write.

      **This is meaningfully slower than the unverified path, not just "an
      extra frame or two" slower**: only the last frame's answer is ever
      used, but reading it back means a non-posted PCIe status read per bit —
      unlike every write in this file, which is posted (fire-and-forget), a
      read has to wait for an actual completion round-trip. Measured
      2026-07-28: verified writes run roughly 0.1–0.3 ms, against ~6 µs
      unverified — call it two orders of magnitude, not a small tax. Treat
      `verify` as a bring-up/diagnostic tool you turn on to confirm the
      hardware, not something to leave enabled across a real closed loop.

      **Confidence note, updated 2026-07-28 against a local copy of SLASEH2A
      (`dac81404.pdf`):** the frame format, R/W convention, 24-bit no-CRC
      framing, the register TYPE column, and the
      read-command-then-second-cycle readback protocol are all confirmed from
      the datasheet text directly (Tables 8-2, 8-3, 8-7). `init` forces
      `SPICONFIG.FSDO=1`, and that turns out to be *required* rather than
      merely convenient: SCLK idles high, so the first SCLK edge inside a
      cycle is a falling one, and only "SDO updates on SCLK falling edges"
      (8.6.4) lines the 24 samples up one-for-one with the 24 output bits.
      Under the reset default `FSDO=0` ("SDO updates on SCLK rising edges")
      the first rising edge doesn't arrive until after the first falling
      edge, so every sample would sit one bit position early and D0 would
      never be seen at all — that resolves the bit-alignment question this
      note used to leave open. One thing is still open: the *exact*
      status-register bit SDO/ACK lands on for this specific AX99100 card,
      which was never independently measured the way the data/control pin
      levels were, and isn't in the DAC's datasheet at all —
      `ack_status_bit`'s default is the generic SPP convention, not a
      bench-confirmed one. Run with `verify: true` and watch the log for
      `readback echo mismatch` before trusting a quiet run as proof anything
      is actually being checked.
  - `ack_status_bit` (integer) (optional)
    - Status-register bit carrying SDO (default 6, the standard SPP `/ACK`
      bit position). See the confidence note under `verify`.
  - `crc` (boolean) (optional, `spi` link only)
    - Use 32-bit frames with an appended CRC-8 instead of plain 24-bit
      ones, and turn on `SPICONFIG.CRC-EN` so the device itself rejects
      (rather than silently accepting) any write whose CRC doesn't check
      out (default false). Independent of `verify` — protects *every*
      write, not just ones being read back, so a marginal connection can't
      quietly land a corrupted voltage — and when `verify` *is* on it
      strictly improves it, since a CRC-mode echo carries the previous
      cycle's data (Table 8-5), which is the only way to confirm the value
      written to a write-only register like `DACx`. Requires `reset_at_init`
      (default true): the frame that turns CRC on has to be sent in the old,
      plain format the device is still expecting, which is only guaranteed
      if the soft-reset just put the device back to its `CRC-EN=0` default.
      Setting `crc` without `reset_at_init` now warns at init rather than
      failing silently later. That same enabling frame is deliberately sent
      *unverified*: `CRC-EN` takes effect on its own SYNC rising edge, so a
      verified write's trailing read-command and NOP frames would still be
      24 bits long while the device had already started demanding 32 — and
      the device ignores an access cycle with too few clock edges (8.5.1),
      producing a spurious "unconfirmed". It is re-issued immediately
      afterwards in the new format, where it *can* be verified.
      **Confidence note:** SLASEH2A 8.5.3 gives the polynomial
      (`x8+x2+x+1`, i.e. `0x07`) and no worked example — but it also says
      "if no error exists, the CRC remainder is zero", which pins down the
      initial remainder (`0x00`) and rules out a final XOR (textbook
      CRC-8-ATM/HEC adds `0x55`, which would leave a constant non-zero
      remainder instead). `crc8_atm()` matches that, verified by
      recomputation over content-plus-appended-CRC; bit order (MSB first,
      the natural one for this bus) is the one parameter still taken on
      faith. When `verify` is also on, the readback path recomputes the CRC
      over the echoed content and compares it to what the chip sent, *and*
      checks the `CRC-ERROR` bit the device returns in bit 30 of the
      following cycle (Tables 8-5, 8-6) — that bit is how the device reports
      that it rejected a frame, and it was previously masked away and never
      looked at, so a rejected write went by silently.

### Outputs

  - `channel` (integer or array) (optional)
    - DAC channel(s): 0 = DACA, 1 = DACB, 2 = DACC, 3 = DACD (default 0).
  - `index` (integer or array) (optional)
    - Pipeline vector element(s) to command (default 0).
  - `scale` (number or array) (optional)
    - Volts per pipeline unit (default 1.0).
  - `offset` (number or array) (optional)
    - Volts at zero command (default 0.0).
  - `range` (string or array) (optional)
    - Output range per channel (default `"0-10"`): `"0-5"`, `"0-6"`,
      `"0-10"`, `"0-12"`, `"0-20"`, `"0-24"`, `"0-40"`, `"+-5"`, `"+-6"`,
      `"+-10"`, `"+-12"`, `"+-20"` (the full set of valid `DACRANGE` codes,
      SLASEH2A table 8-16). Codes are MSB-aligned straight binary across the
      selected range, in 2^16 steps — see "Codes and volts" above for why the
      top of the range is one LSB out of reach.
  - `start_delay` (number or array) (optional)
    - Seconds to hold a channel at its `offset` before letting the pipeline
      command through (default 0). Same purpose as in `piplate_bridge`: give
      a coarse stage time to walk the beam near centre. Give the upstream
      `pid` a matching `start_delay`, or its integrator winds up during the
      hold.
  - `skip_unchanged` (boolean) (optional)
    - Only send a channel when its 16-bit code changes (default false).
  - `range_max` (string or array of strings) (optional)
    - Widest range auto-ranging may select, per channel (default: the same as
      `range`, i.e. auto-ranging **off**). When it names a wider range, the
      device switches a channel's range at runtime: it widens the moment the
      commanded voltage no longer fits, and narrows again once the command has
      sat inside `shrink_frac` of the narrower range for `shrink_dwell`
      consecutive iterations.
    - It walks a **ladder**, always choosing the narrowest range that both
      contains the voltage and lies inside the `[range, range_max]` window. So
      `range: "0-5"`, `range_max: "0-10"` goes 0-5 → 0-6 → 0-10 as needed;
      0-10 is a ceiling, not a jump.
    - The selection boundaries have a standalone regression test in
      `devices/parport_dac_range_test.c`. For an end-to-end hardware check, use
      `contrib/calibration-scripts/configurations/conf_parport_dac_autorange.json` only with the mirror/FSM input
      disconnected: its slow ±7 V command deliberately exercises
      ±5 → ±6 → ±10 and the hysteretic return path.
    - The window also fixes the ladder's **polarity**: `0-10` does not contain
      `-5..+5`, so a unipolar base with a unipolar ceiling never wanders into
      a bipolar range. If a channel does need to go below ground, name a
      bipolar `range_max` (`"+-5"` *does* contain `0-5`), and the ladder
      becomes 0-5 → ±5.
    - `range_max` must *contain* `range`, and init fails loudly if it doesn't.
      That window is also what keeps auto-ranging from selecting a range the
      **supply** cannot deliver — the part will happily accept `0-20` on a
      +12 V rail and simply saturate ~1.5 V short while the loop sees nothing
      wrong. SLASEH2A 7.13: unipolar needs `AVDD ≥ VMAX + 1.5 V`, bipolar also
      `AVSS ≤ VMIN − 1.5 V`. On ±12 V rails that permits 0-5, 0-6, 0-10, ±5,
      ±6 and ±10 — with the 10 V ones leaving only 0.5 V of margin.
    - **`AVSS` and a unipolar base range: an unavoidable trade.** SLASEH2A 7.5
      gives the unipolar accuracy figures (TUE ±0.07 %FSR, offset error,
      zero-code error 0.15 %FSR) at `AVSS = 0 V`; the `−21.5 V ≤ AVSS < 0`
      conditions belong to the *bipolar* rows. There is no specified
      unipolar-plus-negative-`AVSS` combination.
      - If the ladder never needs a bipolar range, **ground `AVSS`** — that is
        the characterized condition, and it halves analog dissipation
        (~96 mW vs ~192 mW on ±12 V, `AIDD`/`AISS` 8 mA typ).
      - If a bipolar range is in the ladder (`range_max: "+-10"` for a
        mirror that must be driven negative), `AVSS` **must** be negative, and
        the unipolar base range then runs outside its characterized condition.
        Grounding `AVSS` would make the bipolar rungs physically unreachable,
        so there is no way to have both. In mitigation, the bipolar TUE at
        negative `AVSS` is spec'd *tighter* (±0.05 %FSR), which suggests the
        part is comfortable with footroom below ground — but it is not
        guaranteed, and belongs in the run notes.
    - Reaching ±10 V needs `AVSS ≤ −11.25 V` (at ≤10 mA) or `≤ −11.5 V` (at
      ≤15 mA), from the footroom spec below. A −12 V rail leaves only
      0.5–0.75 V of margin, and a nominal 12 V rail at −5 % (11.4 V) **fails**
      the 15 mA case — meter the rails under load rather than trusting the
      label.
    - Runtime range changes require `link: "spi"`. The Pico sample protocol
      does not carry DACRANGE commands, so initialization rejects auto-ranging
      on that link instead of silently sending an invalid frame.
    - The 1.25–1.5 V output headroom/footroom figure is **not** a floor on the
      output: 7.5 conditions it on `−10 mA ≤ load ≤ 10 mA` (1.25 V) or ±15 mA
      (1.5 V). Unloaded — the EC table's own condition — zero code sits within
      the 0.15 %FSR zero-code error, ~7.5 mV of ground on a 5 V range. Only a
      load pulling ~10 mA that must still swing to 0 V justifies a negative
      `AVSS`, and it costs the characterized condition above.
    - **A range change shows one intermediate value, and its size depends on
      the write order.** `DACRANGE` is a resistor gain network (SLASEH2A
      8.3.2), so it takes effect the instant the write lands and re-interprets
      whatever code is *already* in the active register. There is no atomic
      range-plus-code update, and synchronous mode does not help — `DACRANGE`
      acts on the active register and ignores `SYNCCONFIG`. So for one frame
      (~6.4 µs) the output shows neither the old voltage nor the new one.
    - It cannot be reduced to *zero*: one code expresses the same voltage in
      two ranges only where they intersect, and solving
      `vmin_A + span_A·c = vmin_B + span_B·c` gives exactly one such voltage
      (0 V for two bipolar or two unipolar ranges). Away from it, something
      moves.
    - What it **is** reduced to is the smaller of the two orders, chosen per
      transition by comparing how far each candidate intermediate falls
      outside `[v_current, v_target]`:
      - *`DACRANGE` first* reinterprets the **old** code. Widening 0-5→0-10 at
        4.9 V puts **9.8 V** on the output — a 2× overshoot, in the direction
        the command was already heading, which is the worst way to throw a
        mirror.
      - *code first* writes the **new** code while the old range is still
        active, so it reads as the target scaled *down* by the range ratio —
        the output moves **toward zero** and then lands exactly on target when
        `DACRANGE` follows. Same case: **3.5 V** instead of 9.8 V, and it
        costs one frame less because no third write is needed.

      Measured over the implementation: 0-5→0-10 at 4.9→7.0 V gives 3.5 V
      rather than 9.8 V; the gentler 0-5→0-6 ladder step at 4.9→5.5 V gives
      4.58 V rather than 5.88 V, an excursion of only 0.32 V. Because
      `pick_range` always takes the narrowest range that fits, the 0-6 step is
      what actually happens for modest excursions.
    - **Holding LDAC does not avoid this**, which is the intuitive thing to
      try. SLASEH2A 8.3.3 gates only "the DAC **data** registers" (buffer →
      active); `DACRANGE` sets the programmable-gain output buffer (8.3.2),
      downstream of the ladder. Output is `gain(DACRANGE) × ladder(active)`,
      so freezing the active register freezes one factor and leaves the other
      free to change.
    - For the same reason `sync_update` forces the old, larger order: a `DACx`
      write only loads the buffer register there, so it cannot pre-position
      the active register at all. Expect the bigger intermediate if you
      enable both.
    - Widening is still immediate (clipping the command is worse) and
      narrowing still deliberately reluctant, because a transition is not free
      even at the smaller size.
    - A range change also invalidates `last_codes`, so `skip_unchanged` cannot
      conclude "unchanged" and leave the output sitting at a reinterpreted
      code.
  - `shrink_frac` (float) (optional)
    - Narrow only when the command fits inside this fraction of the candidate
      narrower range (default 0.8). Hysteresis against chatter at a boundary —
      every oscillation would be another one-frame glitch.
  - `shrink_dwell` (integer) (optional)
    - ...and only after it has fitted for this many consecutive iterations
      (default 4000, about 1 s at 3.8 kHz).
  - `sync_update` (boolean) (optional)
    - Put every driven channel in synchronous mode so they all step together
      (default false, i.e. asynchronous; `spi` link only). Each `DACx` write
      loads only that channel's buffer register, and one SOFT-LDAC per
      iteration moves all of them at once — see "Update mode" above.
    - The trade: one extra frame per iteration (~6.4 µs on the `spi` link),
      and the first channel's motion is deferred until the last one has been
      written, so the *group* moves slightly later than the first channel
      would have on its own. What you buy is that the channels stop stepping
      ~6.4 µs apart. Worth it when two channels drive one physical thing — an
      X/Y mirror, a coarse/fine pair — and the skew between them is a real
      disturbance rather than a rounding error against the loop period.
    - With `skip_unchanged`, the trigger is only issued when at least one
      channel was actually written; an iteration where every channel was
      skipped costs nothing.
    - Pointless with a single channel (it logs a warning) and rejected on the
      `pico` link, where the bridge firmware owns the DAC's init and LDAC
      timing — honouring it there would be a silent no-op.
  - `reset_at_init` (boolean) (optional)
    - Issue a DAC soft reset before configuring (default true; `spi` link
      only).

Example
-------

Two fine channels, ±2 V about a 2 V bias on a 0–10 V range, held at bias for
the first 30 s while a coarse stage settles:

```json
{
 "uri": "anyloop:parport_dac",
 "params": {
  "pci": "0000:05:00.2",
  "link": "spi",
  "channel": [0, 1],
  "index": [0, 1],
  "range": "0-10",
  "scale": [2, -2],
  "offset": [2, 2],
  "start_delay": [30, 30]
 }
}
```

An X/Y pair that must move together — same two channels, but with the ~6.4 µs
skew between them removed by a single SOFT-LDAC per iteration:

```json
{
 "uri": "anyloop:parport_dac",
 "params": {
  "pci": "0000:05:00.2",
  "link": "spi",
  "channel": [0, 1],
  "index": [0, 1],
  "range": "0-10",
  "scale": [2, -2],
  "offset": [5, 5],
  "sync_update": true
 }
}
```

`contrib/legacy/outdated_scripts/conf_parport_dac_smoke.json` is a retired
standalone bench check: a 0.2 Hz
sine into DACA/DACB with no camera and no mirror.

Measuring rise and fall time
----------------------------

`contrib/camera_pcie_hardware/dac_square.c` puts a square wave on one channel so the analog edge
can be measured on a scope. It is standalone — it maps the BAR and bit-bangs
the same frames this device does, without anyloop in the picture:

```
ninja -C build                           # or: gcc -O2 -o dac_square \
                                         #       contrib/camera_pcie_hardware/dac_square.c -lm
sudo build/dac_square                    # 1 kHz, 0 -> 5 V full-scale on DACA
sudo build/dac_square --rt               # SCHED_FIFO, less edge jitter
sudo build/dac_square --low 2.4 --high 2.6   # small-signal, the bode drive
sudo build/dac_square --channel 1 --freq 200 # DACB, slower, longer flat tops
sudo build/dac_square --help
```

`dac_square` is a meson target but is **not** installed — run it from the build
tree. `ninja -C build` also rebuilds `anyloop` itself, which is what you need
after any change to `fsp` or `parport_dac`.

It exists as a separate tool rather than a config because `test_source` has
kinds constant/sine/noise but no square, and because a pipeline-generated
edge would carry the loop's scheduling jitter into the very timing being
measured.

**What sets the rise time.** The output steps on *one* edge — the SYNC rising
edge that closes the frame, because every channel is in asynchronous mode (see
Update mode above). What follows is the DAC's own slew and settling, and it
does **not** depend on how fast the frame was clocked in. A slow frame delays
the analog edge; it does not smear it. Worth knowing before distrusting a
reading taken over the `ppdev` backend.

**Trigger.** D3 (DB25 pin 5, unused by the SPI link and never connected to the
EVM) is driven as a digital copy of the wave, and its edge is emitted on the
*same MMIO store* that raises SYNC — so it marks the update instant to within
one store (~123 ns). Probe `DAC_VOUT_x` on one channel and pin 5 on another,
trigger on pin 5, and t=0 is the update instant. That also makes the DAC's
dead time before it starts moving visible, which a scope triggered on the
analog edge itself cannot show. `--trig-bit` moves it to D4–D7 or `-1` turns
it off.

**Levels.** The defaults stay at 0 V and above deliberately. If the EVM is
strapped unipolar (J11 2–3, `DAC_AVSS` grounded) the part cannot drive below
ground whatever range is selected, so a negative `--low` measures the output
amplifier against its rail rather than a real settling edge. Check J11 and the
J17 supplies before widening the step. The tool prints the exact codes, the
achieved levels and the 10 %–90 % thresholds for the step it is about to
drive, so the numbers to read off the scope are on screen before it starts.

It refuses a frequency/duty combination whose shorter phase falls below the
2.4 µs `tDACWAIT` minimum, and reports at exit if any cycle started late —
i.e. if the wave was not actually at the requested frequency. Ctrl-C parks the
output at `--low` rather than leaving it wherever the last frame put it. It
unbinds `parport_pc` like the driver does, and prints the rebind command on
the way out.
