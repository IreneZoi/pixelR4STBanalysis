/* -------------------------------------------------------------
 *
 *  file:        r4stest.h
 *
 *  description: globals for PSI46V2 Wafer tester
 *
 *  author:      Beat Meier
 *  modified:    24.7.2017
 *
 *  rev:
 *
 * -------------------------------------------------------------
 */

#ifndef R4STEST_H
#define R4STEST_H

#include "config.h"
#include "pixel_dtb.h"
#include "settings.h"
#include "protocol.h"
#include "roc.h"

#define VERSIONINFO TITLE " " VERSION " (" TIMESTAMP ")"

// global variables:

extern int nEntry; // counts the entries in the log file

extern CTestboard tb;
extern CSettings settings;  // global settings
extern CProtocol Log;  // log file
extern CRoc roc;

void cmd();

#endif
