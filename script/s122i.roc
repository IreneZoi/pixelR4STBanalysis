
chip 122 # new adapter, irradiated

id   200  mA # current limit
vd   3000 mV # DTB
vdig 2600 mV # Reg

ia   300  mA # current limit
va   2800 mV # DTB
vana 2500 mV # 182 mA, against noise

vref 600 mV # PSI: 250
rgpr 1100 mV # PSI: 600  James: 900
rgsh 700 mV # PSI: 600  James: 670
vcal 400 mV

vaux1 0 mV
vaux2 0 mV

adcdel 10
#daqena 1 # 1 = 40 MHz, int trg
daqena 5 # 1 = 40 MHz, ext trg

#cal 7 2
cal 159 159 # outside pixel array = quiet mode

#hold 15 # test pulses
hold 0 # beam

#d1 3 # RBI = SEQ_START
#d1 4 # RBO = SEQ_END
#d1 5 # measure
 d1 11  # 11 = HOLD

#d2 10 # 10 = CAL_PULSE, for ext trg pulser
 d2 11 # 11 = HOLD

a1 1 # analog out
a2 6

#seqreadout 0 # row-wise, int trg
#seqreadcol 0 # column-wise, int trg
seqreadcol 1 # column-wise, ext trg
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
