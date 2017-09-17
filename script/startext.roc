
chip 13

vd 2400 mV # less noise
#vdig 2400 mV # less noise
id 200  mA # current limit

va 2000 mV # 125 mA for chip 13
#vana 2000 mV # 125 mA for chip 13
ia 300  mA # current limit

vref 250 mV
rgpr 900 mV # PSI: 600
rgsh 670 mV # PSI: 600
vcal 400 mV

vaux1 0 mV
vaux2 0 mV

adcdel 10
daqena 5 # 5 = ext

cal 1 7
# cal 159 159 # outside pixel array = quiet mode
hold 24

d1 3 # RBI = SEQ_START
#d1 4 # RBO = SEQ_END
#d1 5 # trigger per row

d2 11 # 10 = CAL_PULSE, for ext trg pulser
#d2 11 # 11 = HOLD

a2 6
a1 1 # analog out

seqreadcol 1 # pedestal view ext
#seqcalscan # test pulse

pon
mdelay 500
getia
getid
go
