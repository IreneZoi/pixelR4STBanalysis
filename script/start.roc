vd 2200 mV
id 200  mA

va 2050 mV
ia 300  mA

vref 250 mV
rgpr 600 mV
rgsh 600 mV
vcal 1000 mV

vaux1 0 mV
vaux2 0 mV
vaux3 0 mV

adcdel 10
daqena 2

cal 1 150
hold 11

d1 5 trigger
d2 0
a2 6
a1 1 roc4sens data out

seqreadcol

pon
mdelay 500
getia
getid
go
