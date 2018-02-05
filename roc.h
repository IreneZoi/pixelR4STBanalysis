
#include <iostream>

extern CProtocol Log;  // log file

class CRoc
{
public:
  int chip;
  int VA;
  int VD;
  int Vana;
  int Vdig;
  int Vref;
  int VgPr;
  int VgSh;

  int Hold;
  int ADCdel;

  int Vcal;
  int CalX;
  int CalY;

  int DAQ; // daq enable flag
  bool ext; // external trigger flag

  bool haveGain;
  double p0[155][160];
  double p1[155][160];
  double p2[155][160];
  double p3[155][160];

  void print( bool wlog = 0 )
  {
    std::cout
      << "ROC settings:" << std::endl
      << "chip   " << chip << std::endl
      << "VA     " << VA << std::endl
      << "VD     " << VD << std::endl
      << "Vana   " << Vana << std::endl
      << "Vdig   " << Vdig << std::endl
      << "Vref   " << Vref << std::endl
      << "VgPr   " << VgPr << std::endl
      << "VgSh   " << VgSh << std::endl
      << "hold   " << Hold << std::endl
      << "ADCdel " << ADCdel << std::endl
      << "Vcal   " << Vcal << std::endl
      << "cal x  " << CalX << std::endl
      << "cal y  " << CalY << std::endl
      << "DAQ    " << DAQ << std::endl
      << "ext    " << ext << std::endl
      << "Gain   " << haveGain << std::endl
      << std::endl;
    if( wlog ) {
      Log.printf( "[ROC %i settings]\n", chip );
      //Log.printf( "  chip    %i\n", chip );
      Log.printf( "  VA      %i mV\n", VA );
      Log.printf( "  VD      %i mV\n", VD );
      Log.printf( "  Vana    %i mV\n", Vana );
      Log.printf( "  Vdig    %i mV\n", Vdig );
      Log.printf( "  Vref    %i mV\n", Vref );
      Log.printf( "  VgPr    %i mV\n", VgPr );
      Log.printf( "  VgSh    %i mV\n", VgSh );
      Log.printf( "  Hold    %i*6.25 ns\n", Hold );
      Log.printf( "  ADCdel  %i\n", ADCdel );
    }
  }

};
