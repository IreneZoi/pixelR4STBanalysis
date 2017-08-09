ROC4SENS Software
==============================
This repository contains _roc4senstest_ source code. _roc4senstest_ is a software to operate a _ROC4SENS_ chip using a DTB. The procedures to build and execute the software are explained below.

Build Requirements
------------------------------
The DTB uses a FTDI chip for USB communications. The FTDI D2XX drivers are needed to compile the software. The drivers are available at:
```
http://www.ftdichip.com/Drivers/D2XX.htm\
```
In addition, _ROOT_ libraries and _libusb_ headers are required.

Build
------------------------------
To build to software, change to the its directory and run:
```
make
```
If build is successful, the executable is created in the `bin` directory.

Running the Software
------------------------------
Run the following command to execute the software:
```
bin/r4stest [LogFileName]
```
Type `help` for a list of available commands.
