## search combinations of tones in the Laskar solution

E5 = 4.248956  
E5p = 30.6427

E6 = 27.960603
E6p = 121.235372

I6 = -26.330003
I6p = 307.6644
I7 = -2.985415
I7p = 141.649322

test = 51.672231
testp = 42.373249

2*E6 - E5 - test
2*E6p - E5p - testp

test = 28.32799
testp = 28.724860

2*E6 - E5 + I6 - I7 - test
2*E6p - E5p + I6p - I7p - testp

test = -80.620542
testp = 128.898
-E6 + 2*I6
-E6p + 2*I6p - 360

# le degre de tolerance que semble s'autoriser laskar 88 (le laskar d'il y a quarante ans...) est de 1.e-4 arcsec/year sans trop se soucier des phases.  C'est a dire 10-7 rad/kyr


1.e-4/360/3600*2*pi*1000



