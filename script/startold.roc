
chip 107

vd   2400 mV # DTB
id   200  mA # current limit

va   2150 mV # DTB
ia   300  mA # current limit

vref 250 mV
rgpr 900 mV # PSI: 600
rgsh 670 mV # PSI: 600
vcal 400 mV

vaux1 0 mV
vaux2 0 mV

adcdel 10
daqena 1 # 1 = 40 MHz
#daqena 3 # 3 = 20 MHz

cal 7 2

hold 24

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
mdelay 500
getia
getid
go