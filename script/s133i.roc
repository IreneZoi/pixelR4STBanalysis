
chip 133 # new adapter, irradiated, chilled: 42 uA at 120 V

id   200  mA # current limit
vd   2800 mV # DTB
vdig 2400 mV # Reg

ia   300  mA # current limit
va   2800 mV # DTB
vana 2320 mV # 125 mA for 133

vref 250 mV
rgpr 800 mV # was 600
rgsh 600 mV # was 400
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
