
chip 124 # new adapter, irradiated, chilled

id   200  mA # current limit
vd   2800 mV # DTB
vdig 2400 mV # Reg

ia   300  mA # current limit
va   2800 mV # DTB
vana 2500 mV # 166 mA against hole at right

#vref   0 mV # RMS 3, low gain at small pulses
#vref 100 mV # RMS 5
vref 250 mV # RMS 11
rgpr 600 mV # PSI: 600  James: 900
rgsh 700 mV # PSI: 600  James: 670
vcal 400 mV

vaux1 0 mV
vaux2 0 mV

adcdel 10
daqena 1 # 1 = 40 MHz
#daqena 3 # 3 = 20 MHz

cal 7 2

hold 20

 d1 3 # RBI = SEQ_START
#d1 4 # RBO = SEQ_END
#d1 5 # measure

#d2 10 # 10 = CAL_PULSE, for ext trg pulser
 d2 11 # 11 = HOLD

a1 1 # analog out
a2 6

#seqreadout 0 # row-wise, int trg
seqreadcol 0 # column-wise, int trg
#seqcalscan # test pulse

pon
flush
mdelay 500
getvd
getva

hvon
getia
getid

go