/* -------------------------------------------------------------
 *
 *  file:        cmd_analyzer.cpp
 *
 *  description: command line interpreter
 *               experimental function
 *
 *  author:      Beat Meier
 *  modified:    25.7.2017
 *
 *  rev:
 *
 * -------------------------------------------------------------
 */

CMD_REG(getimg, "", "Read an image from r4s")
CMD_REG(getped, "", "pedestal image from r4s")
CMD_REG(getcal, "", "pulse height image from r4s")
CMD_REG(scancal, "", "scan Vcal from r4s")
CMD_REG(seqreadout, "", "Load measure -> readout sequence")
CMD_REG(seqcalscan, "", "Load calibrate scan sequence")

//CMD_REG(gui, "", "Start graphical user interface");
