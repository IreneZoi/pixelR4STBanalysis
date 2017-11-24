/* -------------------------------------------------------------
 *
 *  author:      Beat Meier, PSI, 25.7.2017
 *  modified:    Aug 2017
 *
 * -------------------------------------------------------------
 */

CMD_REG( chip, "<id>", "define chip, load gain")
CMD_REG( getimg, "<Nev> <Npx> <prnt>", "Read an image")
CMD_REG( takeraw, "<prnt>", "take raw data, write to file")
CMD_REG( td, "Ntrg Nev plane", "take roi data, write to file")
CMD_REG( getped, "", "pedestal image")
CMD_REG( getcal, "", "pulse height image")
CMD_REG( scancal, "", "scan Vcal")
CMD_REG( scanhold, "", "scan hold delay")
CMD_REG( seqreadout, "", "Load measure -> readout sequence")
CMD_REG( seqreadcol, "", "Load measure -> readout columns") // FW 0.7
CMD_REG( seqcalscan, "", "Load calibrate scan sequence")
CMD_REG( scanva, "", "scan VA from r4s")

//CMD_REG( gui, "", "Start graphical user interface");
