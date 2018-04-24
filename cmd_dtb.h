/* -------------------------------------------------------------
 *
 *  file:        cmd_dtb.h
 *
 *  description: command line interpreter
 *               DTB base functions
 *
 *  author:      Beat Meier
 *  modified:    21.4.2014
 *
 *  rev:
 *
 * -------------------------------------------------------------
 */

// =======================================================================
//  connection, communication, startup commands
// =======================================================================

HELP_CAT("dtb");

CMD_REG(scan, "", "Get infos of all connected DTBs")
CMD_REG(open, "[<name>]", "open a DTB (with name)")
CMD_REG(close, "", "close DTB connection")
CMD_REG(rpclink, "", "link all DTB functions")
CMD_REG(welcome, "", "blink with LEDs")
CMD_REG(setled, "<mask>", "set atb LEDs")
CMD_REG(log, "<text>", "writes text to log file")
CMD_REG(upgrade, "<filename>", "upgrade DTB")
CMD_REG(rpcinfo, "", "list all DTB functions")
CMD_REG(info, "", "show detailed DTB info")
CMD_REG(ver, "", "shows DTB software version number")
CMD_REG(version, "", "shows DTB software version")
CMD_REG(boardid, "", "get board id")
CMD_REG(init, "", "inits the testboard")
CMD_REG(flush, "", "flushes usb buffer")
CMD_REG(clear, "", "clears usb data buffer")


// =======================================================================
//  delay commands
// =======================================================================

CMD_REG(udelay, "<us>", "waits <us> microseconds")
CMD_REG(mdelay, "<ms>", "waits <ms> milliseconds")

// =======================================================================
//  test board commands
// =======================================================================

CMD_REG(clksrc, "<source>", "Select clock source (1=ext, 0=int)")
CMD_REG(clkok, "clkok", "Check if ext clock is present")

CMD_REG(clklvl, "<level>", "clk signal level")
CMD_REG(sdalvl, "<level>", "sda signel level")
CMD_REG(ctrlvl, "<level>", "ctr signel level")
CMD_REG(tinlvl, "<level>", "tin signel level")
CMD_REG(lvds, "", "LVDS inputs")
CMD_REG(lcds, "", "LCDS inputs")
// CMD_REG(tout, "", "")
// CMD_REG(trigout, "", "")
CMD_REG(pon, "", "switch ROC power on")
CMD_REG(poff, "", "switch ROC power off")
CMD_REG(va, "<mV>", "set VA in mV")
CMD_REG(vd, "<mV>", "set VD in mV")
CMD_REG(ia, "<mA>", "set IA in mA")
CMD_REG(id, "<mA>", "set ID in mA")
CMD_REG(getva, "", "get VA in V")
CMD_REG(getvd, "", "get VD in V")
CMD_REG(getia, "", "get IA in mA")
CMD_REG(getid, "", "get ID in mA")
CMD_REG(hvon, "", "switch HV on")
CMD_REG(hvoff, "", "switch HV off")
CMD_REG(status, "", "shows testboard status")
CMD_REG(d1, "<signal>", "assign signal to D1 output")
CMD_REG(d2, "<signal>", "assign signal to D2 outout")
CMD_REG(a1, "<signal>", "assign analog signal to A1 output")
CMD_REG(a2, "<signal>", "assign analog signal to A2 outout")
CMD_REG(probeadc, "<signal>", "assign analog signal to ADC")
CMD_REG(probeadc2, "<signal>", "assign analog signal to ADC gain 2")

// === DAQ ==================================================================

HELP_CAT("daq")

CMD_REG(dopen, "<buffer size>", "Open DAQ and allocate memory")
CMD_REG(dclose, "", "Close DAQ")
CMD_REG(dstart, "", "Enable DAQ")
CMD_REG(dstop, "", "Disable DAQ")
CMD_REG(dsize, "", "Show DAQ buffer fill state")
CMD_REG(dread, "", "Read Daq buffer and show as raw data")


// === ROC4sens =============================================================

HELP_CAT("r4s")

CMD_REG(show, "", "print ROC settings")
CMD_REG(cal, "<x> <y>", "Set calibrate to pixel (x.y)")
CMD_REG(cal2, "<x> <y>", "Set calibrate to pixel (x.y, x+1.y)")
CMD_REG(cal4, "<x> <y>", "Set calibrate to pixel (x.y, x+3.y)")
CMD_REG(hold, "<t>", "Set hold position in 6.25 ns steps")
CMD_REG(adcdel, "<t>", "Set ADC sampling point in 5ns steps")
CMD_REG(daqena, "<ext, slow, enable>", "Enable DAQ")

CMD_REG(go, "", "Start ROC4sens measurement sequence")
CMD_REG(goloop, "", "Start readout loop. Stop by pressing any key")
CMD_REG(running, "", "Check if measurement sequence running")

CMD_REG(vdig, "<mV>", "Set Roc4sens Vaux3")
CMD_REG(vana, "<mV>", "Set Roc4sens Vaux3")
CMD_REG(vcal, "<mV>", "Set Roc4sens Vcal")
CMD_REG(rgpr, "<mV>", "Set Roc4sens RgPr")
CMD_REG(rgsh, "<mV>", "Set Roc4sens RgSh")
CMD_REG(vref, "<mV>", "Set Roc4sens Vref")
CMD_REG(vaux1, "<mV>", "Set Roc4sens Vaux1")
CMD_REG(vaux2, "<mV>", "Set Roc4sens Vaux2")
