
chip 124 # new adapter, irradiated, chilled

id   200  mA # current limit
vd   2800 mV # DTB
vdig 2400 mV # Reg

ia   300  mA # current limit
va   2800 mV # DTB
vana 2500 mV # 166 mA against hole at right

#vref   0 mV # RMS  3, low gain at small pulses
vref 100 mV # RMS  6
#vref 250 mV # RMS 13
#vref 400 mV # RMS 13
#vref 600 mV # RMS 10.3
#vref 800 mV # RMS 5.8
#vref 999 mV # RMS 3.4 no gain: saturation at -1400 ADC
rgpr 600 mV # PSI: 600  James: 900
rgsh 700 mV # PSI: 600  James: 670
vcal 400 mV

vaux1 0 mV
vaux2 0 mV

adcdel 10
daqena 5 # 1 = 40 MHz, ext trg

cal 159 159 # outside pixel array = quiet mode

hold 0 # beam

 d1 3 # RBI = SEQ_START
#d1 4 # RBO = SEQ_END
#d1 5 # measure
 d1 11  # 11 = HOLD

#d2 10 # 10 = CAL_PULSE, for ext trg pulser
 d2 11 # 11 = HOLD

a1 1 # analog out
a2 6

#seqreadout 1 # row-wise, ext trg
seqreadcol 1 # column-wise, ext trg

pon
flush
mdelay 500
getvd
getva

hvon
getia
getid

go
