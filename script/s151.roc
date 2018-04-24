
chip 151

id   200  mA # current limit
vd   2800 mV # DTB
vdig 2400 mV # Reg

ia   300  mA # current limit
va   2600 mV # DTB, must be 350 above vana
vana 2120 mV # 125 mA

vref 250 mV
rgpr 900 mV # PSI: 600
rgsh 670 mV # PSI: 600
vcal 400 mV

vaux1 0 mV
vaux2 0 mV

adcdel 10
daqena 5 # 5 = ext 40 MHz

cal 159 159 # outside pixel array = quiet mode
hold 0  # TB21

#d1  3  # RBI = SEQ_START
#d1  4  # RBO = SEQ_END
 d1  5  # trigger per row
 d1 11  # 11 = HOLD

#d2 10  # 10 = CAL_PULSE, for ext trg pulser
 d2 11  # 11 = HOLD

a1 1 # analog out
a2 6

#seqreadout 1 # row-wise, ext trg
seqreadcol 1 # column-wise ext trg

pon
mdelay 500  # contains flush
getvd
getva

hvon
mdelay 500  # contains flush
getia
getid

go  # one readout cycle to initialize D1 = Hold
