
chip 119 # new adapter, irradiated, chilled

id   200  mA # current limit
vd   3200 mV # DTB, must be 350 mV above Vdig
vdig 2800 mV # Regulated

ia   300  mA # current limit
va   2900 mV # DTB, must be 350 mV above Vana
vana 2440 mV # 125 mA at -22 deg

vref 400 mV # baseline at 500 ADC
rgpr 800 mV # PSI: 600  James: 900
rgsh 600 mV # PSI: 600  James: 670
vcal 200 mV

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
