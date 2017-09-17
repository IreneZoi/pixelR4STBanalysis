
#include <iostream>

class CRoc
{
public:
  int chip;
  int VA;
  int VD;
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

  void print()
  {
    std::cout
      << "ROC settings:" << std::endl
      << "chip   " << chip << std::endl
      << "VA     " << VA << std::endl
      << "VD     " << VD << std::endl
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
  }

};
