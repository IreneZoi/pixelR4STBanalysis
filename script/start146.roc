
chip 146

id   200  mA # current limit
vd   2800 mV # DTB
vdig 2400 mV # Reg

ia   300  mA # current limit
va   2500 mV # DTB
vana 2100 mV # 125 mA

vref 250 mV

#rgpr 600 mV # peak 220 at 16
#rgsh 600 mV # underswing
#rgsh 630 mV # small underswing
#rgpr 700 mV # peak 50+220 at 20
#rgpr 650 mV # peak 30+210 at 17

rgpr 900 mV # James
rgsh 670 mV # James

vcal 400 mV

vaux1 0 mV
vaux2 0 mV

adcdel 10
daqena 1 # 1 = 40 MHz

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
getvd
getva

hvon
mdelay 500
getid
getia

go
