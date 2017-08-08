/* -------------------------------------------------------------
 *
 *  file:        r4stest.cpp
 *
 *  description: main program for ROC4sens test
 *
 *  author:      Beat Meier
 *  modified:    24.7.2017
 *
 *  rev:
 *
 * -------------------------------------------------------------
 */

#include <string>
#include <vector>
#include "r4stest.h"

#include "profiler.h"

using namespace std;


// --- globals ---------------------------------------
int nEntry; // counts the chip tests

CTestboard tb;
CSettings settings;  // global settings
CProtocol Log;


void help()
{
	printf("r4stest <log file name>\n");
}


string filename;


int main(int argc, char* argv[])
{
	string usbId;
	printf(VERSIONINFO "\n");

	if (argc != 2) { help(); return 1; }
	
	filename = argv[1];

	// --- load settings ----------------------------------
	if (!settings.Read("r4stest.ini"))
	{
		printf("error reading \"r4stest.ini\"\n");
		return 2;
	}

	// --- open log file ----------------------------------

	if (!Log.open(filename.c_str()))
	{
		printf("log: error creating file\n");
		return 3;
	}

	// --- open test board --------------------------------
	Log.section("DTB");

	try
	{
		if (!tb.FindDTB(usbId)) {}
		else if (tb.Open(usbId))
		{
			printf("\nDTB %s opened\n", usbId.c_str());
			string info;
			try
			{
				tb.GetInfo(info);
				printf("--- DTB info-------------------------------------\n"
					   "%s"
					   "-------------------------------------------------\n", info.c_str());
				Log.puts(info.c_str());
				tb.Welcome();
				tb.Flush();
			}
			catch(CRpcError &e)
			{
				e.What();
				printf("ERROR: DTB software version could not be identified, please update it!\n");
				tb.Close();
				printf("Connection to Board %s has been cancelled\n", usbId.c_str());
			}
		}
		else
		{
			printf("USB error: %s\n", tb.ConnectionError());
			printf("Connect testboard and try command 'scan' to find connected devices.\n");
			printf("Make sure you have permission to access USB devices.\n");
		}

		Log.flush();

		// --- call command interpreter -----------------------
		nEntry = 0;

		cmd();
		tb.Close();
	}
	catch (CRpcError &e)
	{
		e.What();
	}

	return 0;
}
