  module Parameters_for_INDO_C
    use vast_kind_param, ONLY:  double
    integer      :: isoki(80),nbfai(80)
    real(double) :: zetai(80),zetadi(2,80),zetawti(2,80), &
                    zcoreai(80),betaai(3,80), fgi(24,80)

! INDO parameter structure:
! isok      = 1 (full parameter set exists for element) or 0 (otherwise)
! nfba      = # basis functions on atom (corresponds to natorb)
! zcorea    = nuclear charge (corresponds to tore)
!
! zeta      = sp Slater orbital exponent
! zetad(1)  =  d Slater orbital exponent 1 (d orbitals are a linear combination
!              of two Slater functions)
! zetad(2)  =  d Slater orbital exponent 2
! zetawt(1) =  d Slater orbital coefficient 1
! zetawt(2) =  d Slater orbital coefficient 2
!
! betaa(1)  = Resonance integral for s orbitals
! betaa(2)  = Resonance integral for p orbitals
! betaa(3)  = Resonance integral for d orbitals
!
! fg(1)     = s ioniz pot from d(n-2)s2 or from sp  
! fg(2)     = s ioniz pot from d(n-1)s1
! fg(3)     = p ioniz pot from d(n-2)s2 or sp
! fg(4)     = p ioniz pot from d(n-1)s1
! fg(5)     = d ioniz pot from d(n-2)s2 or sp
! fg(6)     = d ioniz pot from d(n-1)s1
! fg(7)     = d ioniz pot from dns0
! fg(8)     = Probability of d(n-2)s2 or sp 
! fg(9)     = Probability of d(n-1)s1
! fg(10)    = Probability of dns0
!
! fg(11)    = F0(s,s) (=F0sp,F0pp) Slater-Condon one-center-two-electron integral
! fg(12)    = F0(s,d)
! fg(13)    = F0(d,d)
! fg(14)    = G1(s,p) [Higher-order Slater-Condon factors]
! fg(15)    = F2(p,p)
! fg(16)    = G2(s,d)
! fg(17)    = G1(p,d)
! fg(18)    = F2(p,d)
! fg(19)    = G3(p,d)
! fg(20)    = F2(d,d)
! fg(21)    = F4(d,d)
! fg(22)    = R1(sppd)
! fg(23)    = R3(sddd)
! fg(24)    = R2(sdpp)

!
!                    Data for Element   1     
!
      data        isoki(  1)/    1       /
      data        nbfai(  1)/    1       /
      data      zcoreai(  1)/    1.000000/
      data        zetai(  1)/    1.200000/
      data   betaai(  1,  1)/  -12.0000  /
      data      fgi(  1,  1)/  -13.0600  /
      data      fgi(  8,  1)/    1.000000/
      data      fgi( 11,  1)/   12.8500  /
      data      fgi( 12,  1)/   12.8500  /
!
!                    Data for Element   2     
!
      data        isoki(  2)/    0       /
      data        nbfai(  2)/    1       /
      data      zcoreai(  2)/    2.000000/
      data        zetai(  2)/    1.700000/
      data   betaai(  1,  2)/ -100.0000  /
      data      fgi(  8,  2)/    1.000000/
!
!                    Data for Element   3     
!
      data        isoki(  3)/    1       /
      data        nbfai(  3)/    4       /
      data      zcoreai(  3)/    1.000000/
      data        zetai(  3)/    0.650000/
      data   betaai(  1,  3)/   -1.0000  /
      data   betaai(  2,  3)/   -1.0000  /
      data      fgi(  1,  3)/   -5.4100  /
      data      fgi(  3,  3)/   -3.6100  /
      data      fgi(  8,  3)/    1.000000/
      data      fgi( 11,  3)/    4.5700  /
      data      fgi( 12,  3)/    4.5700  /
      data      fgi( 14,  3)/    2.503730/
      data      fgi( 15,  3)/    1.356879/
!
!                    Data for Element   4     
!
      data        isoki(  4)/    1       /
      data        nbfai(  4)/    4       /
      data      zcoreai(  4)/    2.000000/
      data        zetai(  4)/    0.975000/
      data   betaai(  1,  4)/  -13.0000  /
      data   betaai(  2,  4)/  -13.0000  /
      data      fgi(  1,  4)/   -9.3300  /
      data      fgi(  3,  4)/   -5.8800  /
      data      fgi(  8,  4)/    1.000000/
      data      fgi( 11,  4)/    6.7800  /
      data      fgi( 12,  4)/    6.7800  /
      data      fgi( 14,  4)/    3.828126/
      data      fgi( 15,  4)/    2.656354/
!
!                    Data for Element   5     
!
      data        isoki(  5)/    1       /
      data        nbfai(  5)/    4       /
      data      zcoreai(  5)/    3.000000/
      data        zetai(  5)/    1.300000/
      data   betaai(  1,  5)/   -8.0000  /
      data   betaai(  2,  5)/   -8.0000  /
      data      fgi(  1,  5)/  -14.0000  /
      data      fgi(  3,  5)/   -8.2400  /
      data      fgi(  8,  5)/    1.000000/
      data      fgi( 11,  5)/    8.6800  /
      data      fgi( 12,  5)/    8.6800  /
      data      fgi( 14,  5)/    5.401481/
      data      fgi( 15,  5)/    3.480847/
!
!                    Data for Element   6     
!
      data        isoki(  6)/    1       /
      data        nbfai(  6)/    4       /
      data      zcoreai(  6)/    4.000000/
      data        zetai(  6)/    1.625000/
      data   betaai(  1,  6)/  -17.0000  /
      data   betaai(  2,  6)/  -17.0000  /
      data      fgi(  1,  6)/  -19.4200  /
      data      fgi(  3,  6)/  -10.7000  /
      data      fgi(  8,  6)/    1.000000/
      data      fgi( 11,  6)/   11.1100  /
      data      fgi( 12,  6)/   11.1100  /
      data      fgi( 14,  6)/    6.897842/
      data      fgi( 15,  6)/    4.509913/
!
!                    Data for Element   7     
!
      data        isoki(  7)/    1       /
      data        nbfai(  7)/    4       /
      data      zcoreai(  7)/    5.000000/
      data        zetai(  7)/    1.950000/
      data   betaai(  1,  7)/  -26.0000  /
      data   betaai(  2,  7)/  -26.0000  /
      data      fgi(  1,  7)/  -25.5800  /
      data      fgi(  3,  7)/  -13.2500  /
      data      fgi(  8,  7)/    1.000000/
      data      fgi( 11,  7)/   12.0100  /
      data      fgi( 12,  7)/   12.0100  /
      data      fgi( 14,  7)/    8.958454/
      data      fgi( 15,  7)/    6.459559/
!
!                    Data for Element   8     
!
      data        isoki(  8)/    1       /
      data        nbfai(  8)/    4       /
      data      zcoreai(  8)/    6.000000/
      data        zetai(  8)/    2.275000/
      data   betaai(  1,  8)/  -34.0000  /
      data   betaai(  2,  8)/  -34.0000  /
      data      fgi(  1,  8)/  -32.4900  /
      data      fgi(  3,  8)/  -15.8800  /
      data      fgi(  8,  8)/    1.000000/
      data      fgi( 11,  8)/   13.0000  /
      data      fgi( 12,  8)/   13.0000  /
      data      fgi( 14,  8)/   11.815414/
      data      fgi( 15,  8)/    6.902802/
!
!                    Data for Element   9     
!
      data        isoki(  9)/    1       /
      data        nbfai(  9)/    4       /
      data      zcoreai(  9)/    7.000000/
      data        zetai(  9)/    2.600000/
      data   betaai(  1,  9)/  -44.0000  /
      data   betaai(  2,  9)/  -44.0000  /
      data      fgi(  1,  9)/  -40.1400  /
      data      fgi(  3,  9)/  -18.6100  /
      data      fgi(  8,  9)/    1.000000/
      data      fgi( 11,  9)/   14.0000  /
      data      fgi( 12,  9)/   14.0000  /
      data      fgi( 14,  9)/   14.484415/
      data      fgi( 15,  9)/    8.593198/
!
!                    Data for Element  10     
!
      data        isoki( 10)/    0       /
      data        nbfai( 10)/    4       /
      data      zcoreai( 10)/    8.000000/
      data        zetai( 10)/    2.925000/
      data   betaai(  1, 10)/ -100.0000  /
      data   betaai(  2, 10)/ -100.0000  /
      data      fgi(  8, 10)/    1.000000/
      data      fgi( 14, 10)/   12.781917/
      data      fgi( 15, 10)/    9.327345/
!
!                    Data for Element  11     
!
      data        isoki( 11)/    1       /
      data        nbfai( 11)/    4       /
      data      zcoreai( 11)/    1.000000/
      data        zetai( 11)/    0.836000/
      data   betaai(  1, 11)/   -5.0000  /
      data   betaai(  2, 11)/   -5.0000  /
      data      fgi(  1, 11)/   -4.8600  /
      data      fgi(  3, 11)/   -2.8600  /
      data      fgi(  8, 11)/    1.000000/
      data      fgi( 11, 11)/    3.3100  /
      data      fgi( 12, 11)/    3.3100  /
      data      fgi( 13, 11)/    1.6700  /
      data      fgi( 14, 11)/    1.667583/
      data      fgi( 15, 11)/    0.743903/
!
!                    Data for Element  12     
!
      data        isoki( 12)/    1       /
      data        nbfai( 12)/    4       /
      data      zcoreai( 12)/    2.000000/
      data        zetai( 12)/    1.103000/
      data   betaai(  1, 12)/   -6.0000  /
      data   betaai(  2, 12)/   -6.0000  /
      data      fgi(  1, 12)/   -8.1100  /
      data      fgi(  3, 12)/   -4.5500  /
      data      fgi(  8, 12)/    1.000000/
      data      fgi( 11, 12)/    4.7900  /
      data      fgi( 12, 12)/    4.7900  /
      data      fgi( 13, 12)/    2.4300  /
      data      fgi( 14, 12)/    2.476454/
      data      fgi( 15, 12)/    3.273174/
!
!                    Data for Element  13     
!
      data        isoki( 13)/    1       /
      data        nbfai( 13)/    4       /
      data      zcoreai( 13)/    3.000000/
      data        zetai( 13)/    1.370000/
      data   betaai(  1, 13)/   -7.0000  /
      data   betaai(  2, 13)/   -7.0000  /
      data      fgi(  1, 13)/  -11.4200  /
      data      fgi(  3, 13)/   -6.2900  /
      data      fgi(  8, 13)/    1.000000/
      data      fgi( 11, 13)/    6.2100  /
      data      fgi( 12, 13)/    6.2100  /
      data      fgi( 13, 13)/    3.4200  /
      data      fgi( 14, 13)/    3.359095/
      data      fgi( 15, 13)/    1.602491/
!
!                    Data for Element  14     
!
      data        isoki( 14)/    1       /
      data        nbfai( 14)/    4       /
      data      zcoreai( 14)/    4.000000/
      data        zetai( 14)/    1.520000/
      data   betaai(  1, 14)/   -9.0000  /
      data   betaai(  2, 14)/   -9.0000  /
      data      fgi(  1, 14)/  -14.7900  /
      data      fgi(  3, 14)/   -8.1000  /
      data      fgi(  8, 14)/    1.000000/
      data      fgi( 11, 14)/    7.5700  /
      data      fgi( 12, 14)/    7.5700  /
      data      fgi( 13, 14)/    4.6300  /
      data      fgi( 14, 14)/    4.812310/
      data      fgi( 15, 14)/    2.262706/
!
!                    Data for Element  15     
!
      data        isoki( 15)/    1       /
      data        nbfai( 15)/    4       /
      data      zcoreai( 15)/    5.000000/
      data        zetai( 15)/    1.730000/
      data   betaai(  1, 15)/  -15.0000  /
      data   betaai(  2, 15)/  -15.0000  /
      data      fgi(  1, 15)/  -18.2300  /
      data      fgi(  3, 15)/   -9.9800  /
      data      fgi(  8, 15)/    1.000000/
      data      fgi( 11, 15)/    8.8600  /
      data      fgi( 12, 15)/    8.8600  /
      data      fgi( 13, 15)/    6.0900  /
      data      fgi( 14, 15)/    1.047788/
      data      fgi( 15, 15)/    2.947716/
!
!                    Data for Element  16     
!
      data        isoki( 16)/    1       /
      data        nbfai( 16)/    4       /
      data      zcoreai( 16)/    6.000000/
      data        zetai( 16)/    1.925000/
      data   betaai(  1, 16)/  -15.0000  /
      data   betaai(  2, 16)/  -15.0000  /
      data      fgi(  1, 16)/  -21.7300  /
      data      fgi(  3, 16)/  -11.9200  /
      data      fgi(  8, 16)/    1.000000/
      data      fgi( 11, 16)/   10.0900  /
      data      fgi( 12, 16)/   10.0900  /
      data      fgi( 13, 16)/    7.7700  /
      data      fgi( 14, 16)/    3.075668/
      data      fgi( 15, 16)/    4.537809/
!
!                    Data for Element  17     
!
      data        isoki( 17)/    1       /
      data        nbfai( 17)/    4       /
      data      zcoreai( 17)/    7.000000/
      data        zetai( 17)/    2.130000/
      data   betaai(  1, 17)/  -11.0000  /
      data   betaai(  2, 17)/  -11.0000  /
      data      fgi(  1, 17)/  -25.2900  /
      data      fgi(  3, 17)/  -13.9300  /
      data      fgi(  8, 17)/    1.000000/
      data      fgi( 11, 17)/   11.2500  /
      data      fgi( 12, 17)/   11.2500  /
      data      fgi( 13, 17)/    9.6800  /
      data      fgi( 14, 17)/    8.802854/
      data      fgi( 15, 17)/    6.447161/
!
!                    Data for Element  18     
!
      data        isoki( 18)/    0       /
      data        nbfai( 18)/    4       /
      data      zcoreai( 18)/    8.000000/
      data        zetai( 18)/    2.365000/
      data   betaai(  1, 18)/ -100.0000  /
      data   betaai(  2, 18)/ -100.0000  /
      data      fgi(  8, 18)/    1.000000/
      data      fgi( 14, 18)/    7.762259/
      data      fgi( 15, 18)/    5.846134/
!
!                    Data for Element  19     
!
      data        isoki( 19)/    1       /
      data        nbfai( 19)/    4       /
      data      zcoreai( 19)/    1.000000/
      data        zetai( 19)/    1.180000/
      data   betaai(  1, 19)/   -1.0000  /
      data   betaai(  2, 19)/   -1.0000  /
      data      fgi(  1, 19)/   -4.3400  /
      data      fgi(  3, 19)/   -2.7300  /
      data      fgi(  8, 19)/    1.000000/
      data      fgi( 11, 19)/    3.1800  /
      data      fgi( 12, 19)/    3.1800  /
      data      fgi( 13, 19)/    5.0300  /
      data      fgi( 14, 19)/    1.110895/
      data      fgi( 15, 19)/    0.495935/
!
!                    Data for Element  20     
!
      data        isoki( 20)/    1       /
      data        nbfai( 20)/    9       /
      data      zcoreai( 20)/    2.000000/
      data        zetai( 20)/    1.210000/
      data   zetadi(  1, 20)/    1.850000/
      data  zetawti(  1, 20)/    1.000000/
      data   betaai(  1, 20)/    2.0000  /
      data   betaai(  2, 20)/    2.0000  /
      data   betaai(  3, 20)/  -11.4000  /
      data      fgi(  1, 20)/   -6.0300  /
      data      fgi(  2, 20)/   -5.1300  /
      data      fgi(  3, 20)/   -3.9600  /
      data      fgi(  4, 20)/   -2.9900  /
      data      fgi(  5, 20)/   -3.4400  /
      data      fgi(  6, 20)/   -3.4400  /
      data      fgi(  7, 20)/ -100.0000  /
      data      fgi(  8, 20)/    0.960800/
      data      fgi(  9, 20)/    0.039200/
      data      fgi( 11, 20)/    3.2500  /
      data      fgi( 12, 20)/    4.0000  /
      data      fgi( 13, 20)/    6.0300  /
      data      fgi( 14, 20)/    1.562197/
      data      fgi( 15, 20)/    0.288262/
      data      fgi( 16, 20)/    0.462460/
      data      fgi( 17, 20)/    0.730265/
      data      fgi( 18, 20)/    0.555448/
      data      fgi( 19, 20)/   -0.029508/
      data      fgi( 20, 20)/    2.343295/
      data      fgi( 21, 20)/    1.177847/
      data      fgi( 22, 20)/    2.146212/
      data      fgi( 23, 20)/    2.216599/
      data      fgi( 24, 20)/    1.575989/
!
!                    Data for Element  21     
!
      data        isoki( 21)/    1       /
      data        nbfai( 21)/    9       /
      data      zcoreai( 21)/    3.000000/
      data        zetai( 21)/    1.230000/
      data   zetadi(  1, 21)/    4.222400/
      data   zetadi(  2, 21)/    1.746500/
      data  zetawti(  1, 21)/    0.359220/
      data  zetawti(  2, 21)/    0.766010/
      data   betaai(  1, 21)/   -1.0000  /
      data   betaai(  2, 21)/   -1.0000  /
      data   betaai(  3, 21)/  -18.0000  /
      data      fgi(  1, 21)/   -6.7200  /
      data      fgi(  2, 21)/   -5.8300  /
      data      fgi(  3, 21)/   -4.2000  /
      data      fgi(  4, 21)/   -3.4300  /
      data      fgi(  5, 21)/   -8.1600  /
      data      fgi(  6, 21)/   -4.8500  /
      data      fgi(  7, 21)/ -100.0000  /
      data      fgi(  8, 21)/    0.939900/
      data      fgi(  9, 21)/    0.060100/
      data      fgi( 11, 21)/    3.8900  /
      data      fgi( 12, 21)/    4.7100  /
      data      fgi( 13, 21)/    7.0200  /
      data      fgi( 14, 21)/    1.500205/
      data      fgi( 15, 21)/    0.619919/
      data      fgi( 16, 21)/    0.727785/
      data      fgi( 17, 21)/    0.700509/
      data      fgi( 18, 21)/    1.363822/
      data      fgi( 19, 21)/    0.274004/
      data      fgi( 20, 21)/    3.657524/
      data      fgi( 21, 21)/    1.810164/
      data      fgi( 22, 21)/    1.961795/
      data      fgi( 23, 21)/    2.115305/
      data      fgi( 24, 21)/    1.424637/
!
!                    Data for Element  22     
!
      data        isoki( 22)/    1       /
      data        nbfai( 22)/    9       /
      data      zcoreai( 22)/    4.000000/
      data        zetai( 22)/    1.300000/
      data   zetadi(  1, 22)/    4.670000/
      data   zetadi(  2, 22)/    1.986000/
      data  zetawti(  1, 22)/    0.364610/
      data  zetawti(  2, 22)/    0.755610/
      data   betaai(  1, 22)/   -1.0000  /
      data   betaai(  2, 22)/   -1.0000  /
      data   betaai(  3, 22)/  -19.0000  /
      data      fgi(  1, 22)/   -7.2800  /
      data      fgi(  2, 22)/   -6.3400  /
      data      fgi(  3, 22)/   -4.4800  /
      data      fgi(  4, 22)/   -3.7500  /
      data      fgi(  5, 22)/   -9.0700  /
      data      fgi(  6, 22)/   -5.9300  /
      data      fgi(  7, 22)/ -100.0000  /
      data      fgi(  8, 22)/    0.906900/
      data      fgi(  9, 22)/    0.093100/
      data      fgi( 11, 22)/    4.5000  /
      data      fgi( 12, 22)/    5.3800  /
      data      fgi( 13, 22)/    7.9800  /
      data      fgi( 14, 22)/    1.624189/
      data      fgi( 15, 22)/    0.681911/
      data      fgi( 16, 22)/    0.768700/
      data      fgi( 17, 22)/    0.907562/
      data      fgi( 18, 22)/    1.698579/
      data      fgi( 19, 22)/    1.277034/
      data      fgi( 20, 22)/    5.566875/
      data      fgi( 21, 22)/    3.682321/
      data      fgi( 22, 22)/    1.571618/
      data      fgi( 23, 22)/    1.883618/
      data      fgi( 24, 22)/    1.107171/
!
!                    Data for Element  23     
!
      data        isoki( 23)/    1       /
      data        nbfai( 23)/    9       /
      data      zcoreai( 23)/    5.000000/
      data        zetai( 23)/    1.300000/
      data   zetadi(  1, 23)/    5.052000/
      data   zetadi(  2, 23)/    2.173000/
      data  zetawti(  1, 23)/    0.373780/
      data  zetawti(  2, 23)/    0.745640/
      data   betaai(  1, 23)/   -1.0000  /
      data   betaai(  2, 23)/   -1.0000  /
      data   betaai(  3, 23)/  -20.0000  /
      data      fgi(  1, 23)/   -7.7300  /
      data      fgi(  2, 23)/   -6.7100  /
      data      fgi(  3, 23)/   -4.7700  /
      data      fgi(  4, 23)/   -3.9500  /
      data      fgi(  5, 23)/   -9.8900  /
      data      fgi(  6, 23)/   -6.7700  /
      data      fgi(  7, 23)/ -100.0000  /
      data      fgi(  8, 23)/    0.839500/
      data      fgi(  9, 23)/    0.160500/
      data      fgi( 11, 23)/    5.0700  /
      data      fgi( 12, 23)/    6.0100  /
      data      fgi( 13, 23)/    8.9100  /
      data      fgi( 14, 23)/    1.872156/
      data      fgi( 15, 23)/    0.743903/
      data      fgi( 16, 23)/    0.773659/
      data      fgi( 17, 23)/    0.642236/
      data      fgi( 18, 23)/    1.388619/
      data      fgi( 19, 23)/    0.212012/
      data      fgi( 20, 23)/    6.298380/
      data      fgi( 21, 23)/    4.389029/
      data      fgi( 22, 23)/    1.298119/
      data      fgi( 23, 23)/    1.669320/
      data      fgi( 24, 23)/    0.894537/
!
!                    Data for Element  24     
!
      data        isoki( 24)/    1       /
      data        nbfai( 24)/    9       /
      data      zcoreai( 24)/    6.000000/
      data        zetai( 24)/    1.320000/
      data   zetadi(  1, 24)/    5.138000/
      data   zetadi(  2, 24)/    2.077000/
      data  zetawti(  1, 24)/    0.407140/
      data  zetawti(  2, 24)/    0.732420/
      data   betaai(  1, 24)/   -1.0000  /
      data   betaai(  2, 24)/   -1.0000  /
      data   betaai(  3, 24)/  -21.0000  /
      data      fgi(  1, 24)/   -8.0700  /
      data      fgi(  2, 24)/   -6.9700  /
      data      fgi(  3, 24)/   -5.0400  /
      data      fgi(  4, 24)/   -4.0600  /
      data      fgi(  5, 24)/  -10.6600  /
      data      fgi(  6, 24)/   -7.4300  /
      data      fgi(  7, 24)/ -100.0000  /
      data      fgi(  8, 24)/    0.705200/
      data      fgi(  9, 24)/    0.294800/
      data      fgi( 11, 24)/    5.6000  /
      data      fgi( 12, 24)/    6.6000  /
      data      fgi( 13, 24)/    9.8100  /
      data      fgi( 14, 24)/    1.785368/
      data      fgi( 15, 24)/    0.805895/
      data      fgi( 16, 24)/    0.647196/
      data      fgi( 17, 24)/    0.691830/
      data      fgi( 18, 24)/    1.413416/
      data      fgi( 19, 24)/    0.036823/
      data      fgi( 20, 24)/    7.872975/
      data      fgi( 21, 24)/    4.562606/
      data      fgi( 22, 24)/    1.138395/
      data      fgi( 23, 24)/    1.544197/
      data      fgi( 24, 24)/    0.770846/
!
!                    Data for Element  25     
!
      data        isoki( 25)/    1       /
      data        nbfai( 25)/    9       /
      data      zcoreai( 25)/    7.000000/
      data        zetai( 25)/    1.360000/
      data   zetadi(  1, 25)/    5.767000/
      data   zetadi(  2, 25)/    2.510000/
      data  zetawti(  1, 25)/    0.389840/
      data  zetawti(  2, 25)/    0.729650/
      data   betaai(  1, 25)/   -1.0000  /
      data   betaai(  2, 25)/   -1.0000  /
      data   betaai(  3, 25)/  -22.0000  /
      data      fgi(  1, 25)/   -8.3500  /
      data      fgi(  2, 25)/   -7.1500  /
      data      fgi(  3, 25)/   -5.2700  /
      data      fgi(  4, 25)/   -4.1000  /
      data      fgi(  5, 25)/  -11.4500  /
      data      fgi(  6, 25)/   -7.9900  /
      data      fgi(  7, 25)/ -100.0000  /
      data      fgi(  8, 25)/    0.665200/
      data      fgi(  9, 25)/    0.334800/
      data      fgi( 11, 25)/    6.0900  /
      data      fgi( 12, 25)/    7.1600  /
      data      fgi( 13, 25)/   10.6800  /
      data      fgi( 14, 25)/    2.343295/
      data      fgi( 15, 25)/    0.867887/
      data      fgi( 16, 25)/    0.757541/
      data      fgi( 17, 25)/    0.153740/
      data      fgi( 18, 25)/    0.993111/
      data      fgi( 19, 25)/    0.616200/
      data      fgi( 20, 25)/    8.182935/
      data      fgi( 21, 25)/    4.698988/
      data      fgi( 22, 25)/    1.054290/
      data      fgi( 23, 25)/    1.486488/
      data      fgi( 24, 25)/    0.704629/
!
!                    Data for Element  26     
!
      data        isoki( 26)/    1       /
      data        nbfai( 26)/    9       /
      data      zcoreai( 26)/    8.000000/
      data        zetai( 26)/    1.370000/
      data   zetadi(  1, 26)/    6.068000/
      data   zetadi(  2, 26)/    2.618000/
      data  zetawti(  1, 26)/    0.403790/
      data  zetawti(  2, 26)/    0.719840/
      data   betaai(  1, 26)/   -1.0000  /
      data   betaai(  2, 26)/   -1.0000  /
      data   betaai(  3, 26)/  -23.0000  /
      data      fgi(  1, 26)/   -8.5700  /
      data      fgi(  2, 26)/   -7.2700  /
      data      fgi(  3, 26)/   -5.4200  /
      data      fgi(  4, 26)/   -4.0800  /
      data      fgi(  5, 26)/  -12.3100  /
      data      fgi(  6, 26)/   -8.5300  /
      data      fgi(  7, 26)/ -100.0000  /
      data      fgi(  8, 26)/    0.314300/
      data      fgi(  9, 26)/    0.685700/
      data      fgi( 11, 26)/    6.5400  /
      data      fgi( 12, 26)/    7.6800  /
      data      fgi( 13, 26)/   11.5200  /
      data      fgi( 14, 26)/    2.020937/
      data      fgi( 15, 26)/    0.929879/
      data      fgi( 16, 26)/    0.823253/
      data      fgi( 17, 26)/    0.303760/
      data      fgi( 18, 26)/    0.622399/
      data      fgi( 19, 26)/    0.436423/
      data      fgi( 20, 26)/    7.563016/
      data      fgi( 21, 26)/    4.760980/
      data      fgi( 22, 26)/    0.937362/
      data      fgi( 23, 26)/    1.382743/
      data      fgi( 24, 26)/    0.616758/
!
!                    Data for Element  27     
!
      data        isoki( 27)/    1       /
      data        nbfai( 27)/    9       /
      data      zcoreai( 27)/    9.000000/
      data        zetai( 27)/    1.423000/
      data   zetadi(  1, 27)/    6.386000/
      data   zetadi(  2, 27)/    2.745000/
      data  zetawti(  1, 27)/    0.413330/
      data  zetawti(  2, 27)/    0.712620/
      data   betaai(  1, 27)/   -1.0000  /
      data   betaai(  2, 27)/   -1.0000  /
      data   betaai(  3, 27)/  -31.0000  /
      data      fgi(  1, 27)/   -8.7600  /
      data      fgi(  2, 27)/   -7.3800  /
      data      fgi(  3, 27)/   -5.4800  /
      data      fgi(  4, 27)/   -4.0200  /
      data      fgi(  5, 27)/  -13.3000  /
      data      fgi(  6, 27)/   -9.1000  /
      data      fgi(  7, 27)/ -100.0000  /
      data      fgi(  8, 27)/    0.206500/
      data      fgi(  9, 27)/    0.793500/
      data      fgi( 11, 27)/    6.9600  /
      data      fgi( 12, 27)/    8.1600  /
      data      fgi( 13, 27)/   12.3200  /
      data      fgi( 14, 27)/    2.814434/
      data      fgi( 15, 27)/    0.991871/
      data      fgi( 16, 27)/    0.786058/
      data      fgi( 17, 27)/    0.393029/
      data      fgi( 18, 27)/    0.779859/
      data      fgi( 19, 27)/    0.280204/
      data      fgi( 20, 27)/    7.996959/
      data      fgi( 21, 27)/    5.963624/
      data      fgi( 22, 27)/    0.927827/
      data      fgi( 23, 27)/    1.392741/
      data      fgi( 24, 27)/    0.606753/
!
!                    Data for Element  28     
!
      data        isoki( 28)/    1       /
      data        nbfai( 28)/    9       /
      data      zcoreai( 28)/   10.000000/
      data        zetai( 28)/    1.473000/
      data   zetadi(  1, 28)/    6.706000/
      data   zetadi(  2, 28)/    2.874000/
      data  zetawti(  1, 28)/    0.421200/
      data  zetawti(  2, 28)/    0.706580/
      data   betaai(  1, 28)/   -1.0000  /
      data   betaai(  2, 28)/   -1.0000  /
      data   betaai(  3, 28)/  -35.0000  /
      data      fgi(  1, 28)/   -8.9400  /
      data      fgi(  2, 28)/   -7.5100  /
      data      fgi(  3, 28)/   -5.4100  /
      data      fgi(  4, 28)/   -3.9300  /
      data      fgi(  5, 28)/  -14.4600  /
      data      fgi(  6, 28)/   -9.7900  /
      data      fgi(  7, 28)/ -100.0000  /
      data      fgi(  8, 28)/    0.142100/
      data      fgi(  9, 28)/    0.857900/
      data      fgi( 11, 28)/    7.3400  /
      data      fgi( 12, 28)/    8.6100  /
      data      fgi( 13, 28)/   13.1000  /
      data      fgi( 14, 28)/    2.405287/
      data      fgi( 15, 28)/    1.053863/
      data      fgi( 16, 28)/    0.830692/
      data      fgi( 17, 28)/    0.373191/
      data      fgi( 18, 28)/    0.750102/
      data      fgi( 19, 28)/    0.402948/
      data      fgi( 20, 28)/    9.893912/
      data      fgi( 21, 28)/    6.608340/
      data      fgi( 22, 28)/    0.914218/
      data      fgi( 23, 28)/    1.397009/
      data      fgi( 24, 28)/    0.594082/
!
!                    Data for Element  29     
!
      data        isoki( 29)/    1       /
      data        nbfai( 29)/    9       /
      data      zcoreai( 29)/   11.000000/
      data        zetai( 29)/    1.482000/
      data   zetadi(  1, 29)/    6.795000/
      data   zetadi(  2, 29)/    2.765000/
      data  zetawti(  1, 29)/    0.447290/
      data  zetawti(  2, 29)/    0.696830/
      data   betaai(  1, 29)/   -1.0000  /
      data   betaai(  2, 29)/   -1.0000  /
      data   betaai(  3, 29)/  -40.0000  /
      data      fgi(  1, 29)/   -9.1300  /
      data      fgi(  2, 29)/   -7.6900  /
      data      fgi(  3, 29)/   -5.1800  /
      data      fgi(  4, 29)/   -3.8400  /
      data      fgi(  5, 29)/  -15.8700  /
      data      fgi(  6, 29)/  -10.6700  /
      data      fgi(  7, 29)/ -100.0000  /
      data      fgi(  8, 29)/    0.095600/
      data      fgi(  9, 29)/    0.904400/
      data      fgi( 11, 29)/    7.6800  /
      data      fgi( 12, 29)/    9.0100  /
      data      fgi( 13, 29)/   13.8400  /
      data      fgi( 14, 29)/    2.566466/
      data      fgi( 15, 29)/    1.115855/
      data      fgi( 16, 29)/    0.552968/
      data      fgi( 17, 29)/    0.696789/
      data      fgi( 18, 29)/    1.326627/
      data      fgi( 19, 29)/    0.859208/
      data      fgi( 20, 29)/   10.660133/
      data      fgi( 21, 29)/    7.186725/
      data      fgi( 22, 29)/    0.815680/
      data      fgi( 23, 29)/    1.301797/
      data      fgi( 24, 29)/    0.521822/
!
!                    Data for Element  30     
!
      data        isoki( 30)/    1       /
      data        nbfai( 30)/    4       /
      data      zcoreai( 30)/    2.000000/
      data        zetai( 30)/    1.509000/
      data   betaai(  1, 30)/  -10.0000  /
      data   betaai(  2, 30)/  -10.0000  /
      data      fgi(  1, 30)/   -9.3600  /
      data      fgi(  2, 30)/   -9.3600  /
      data      fgi(  3, 30)/   -4.7700  /
      data      fgi(  4, 30)/   -4.7700  /
      data      fgi(  8, 30)/    1.000000/
      data      fgi( 11, 30)/    7.9800  /
      data      fgi( 12, 30)/    7.9800  /
      data      fgi( 13, 30)/   14.5500  /
      data      fgi( 14, 30)/    2.529271/
      data      fgi( 15, 30)/    1.177847/
!
!                    Data for Element  31     
!
      data        isoki( 31)/    0       /
      data        nbfai( 31)/    4       /
      data      zcoreai( 31)/    3.000000/
      data   betaai(  1, 31)/ -100.0000  /
      data   betaai(  2, 31)/ -100.0000  /
      data      fgi(  1, 31)/ -100.0000  /
      data      fgi(  3, 31)/ -100.0000  /
      data      fgi(  8, 31)/    1.000000/
!
!                    Data for Element  32     
!
      data        isoki( 32)/    0       /
      data        nbfai( 32)/    4       /
      data      zcoreai( 32)/    4.000000/
      data   betaai(  1, 32)/ -100.0000  /
      data   betaai(  2, 32)/ -100.0000  /
      data      fgi(  1, 32)/ -100.0000  /
      data      fgi(  3, 32)/ -100.0000  /
      data      fgi(  8, 32)/    1.000000/
!
!                    Data for Element  33     
!
      data        isoki( 33)/    0       /
      data        nbfai( 33)/    4       /
      data      zcoreai( 33)/    5.000000/
      data        zetai( 33)/    2.921000/
      data   betaai(  1, 33)/  -10.0000  /
      data   betaai(  2, 33)/  -10.0000  /
      data      fgi(  1, 33)/  -18.1800  /
      data      fgi(  3, 33)/   -9.1900  /
      data      fgi(  8, 33)/    1.000000/
      data      fgi( 14, 33)/    5.036471/
      data      fgi( 15, 33)/    3.802126/
!
!                    Data for Element  34     
!
      data        isoki( 34)/    0       /
      data        nbfai( 34)/    4       /
      data      zcoreai( 34)/    6.000000/
      data        zetai( 34)/    2.439000/
      data   betaai(  1, 34)/  -11.6700  /
      data   betaai(  2, 34)/  -11.6700  /
      data      fgi(  1, 34)/  -20.9500  /
      data      fgi(  3, 34)/  -10.3700  /
      data      fgi(  8, 34)/    1.000000/
      data      fgi( 14, 34)/    5.630697/
      data      fgi( 15, 34)/    4.230938/
!
!                    Data for Element  35     
!
      data        isoki( 35)/    1       /
      data        nbfai( 35)/    4       /
      data      zcoreai( 35)/    7.000000/
      data        zetai( 35)/    2.638000/
      data   betaai(  1, 35)/   -8.0000  /
      data   betaai(  2, 35)/   -8.0000  /
      data      fgi(  1, 35)/  -23.9400  /
      data      fgi(  3, 35)/  -12.4400  /
      data      fgi(  8, 35)/    1.000000/
      data      fgi( 11, 35)/    9.0800  /
      data      fgi( 12, 35)/    9.0800  /
      data      fgi( 14, 35)/    6.141120/
      data      fgi( 15, 35)/    4.608700/
!
!                    Data for Element  36     
!
      data        isoki( 36)/    0       /
      data        nbfai( 36)/    4       /
      data      zcoreai( 36)/    8.000000/
      data   betaai(  1, 36)/ -100.0000  /
      data   betaai(  2, 36)/ -100.0000  /
      data      fgi(  1, 36)/ -100.0000  /
      data      fgi(  3, 36)/ -100.0000  /
      data      fgi(  8, 36)/    1.000000/
      data      fgi( 13, 36)/    2.4800  /
!
!                    Data for Element  37     
!
      data        isoki( 37)/    0       /
      data        nbfai( 37)/    4       /
      data      zcoreai( 37)/    1.000000/
      data   betaai(  1, 37)/ -100.0000  /
      data   betaai(  2, 37)/ -100.0000  /
      data      fgi(  1, 37)/   -4.4300  /
      data      fgi(  3, 37)/   -2.6500  /
      data      fgi(  8, 37)/    1.000000/
      data      fgi( 11, 37)/    1.5000  /
      data      fgi( 12, 37)/    1.5000  /
      data      fgi( 13, 37)/    3.3800  /
!
!                    Data for Element  38     
!
      data        isoki( 38)/    0       /
      data        nbfai( 38)/    9       /
      data      zcoreai( 38)/    2.000000/
      data  zetawti(  1, 38)/    1.000000/
      data   betaai(  1, 38)/ -100.0000  /
      data   betaai(  2, 38)/ -100.0000  /
      data   betaai(  3, 38)/ -100.0000  /
      data      fgi(  1, 38)/   -5.8400  /
      data      fgi(  2, 38)/   -5.1900  /
      data      fgi(  3, 38)/   -3.7600  /
      data      fgi(  4, 38)/   -3.1600  /
      data      fgi(  5, 38)/   -3.6600  /
      data      fgi(  6, 38)/   -3.6600  /
      data      fgi(  7, 38)/   -2.4900  /
      data      fgi(  8, 38)/    0.939000/
      data      fgi(  9, 38)/    0.055000/
      data      fgi( 10, 38)/    0.007000/
      data      fgi( 11, 38)/    2.0500  /
      data      fgi( 12, 38)/    1.9700  /
      data      fgi( 13, 38)/    4.2800  /
!
!                    Data for Element  39     
!
      data        isoki( 39)/    1       /
      data        nbfai( 39)/    9       /
      data      zcoreai( 39)/    3.000000/
      data        zetai( 39)/    1.275000/
      data   zetadi(  1, 39)/    3.837000/
      data   zetadi(  2, 39)/    1.739000/
      data  zetawti(  1, 39)/    0.286550/
      data  zetawti(  2, 39)/    0.824680/
      data   betaai(  1, 39)/   -1.0000  /
      data   betaai(  2, 39)/   -1.0000  /
      data   betaai(  3, 39)/   -7.0000  /
      data      fgi(  1, 39)/   -6.5500  /
      data      fgi(  2, 39)/   -5.8500  /
      data      fgi(  3, 39)/   -4.1300  /
      data      fgi(  4, 39)/   -3.5800  /
      data      fgi(  5, 39)/   -6.6100  /
      data      fgi(  6, 39)/   -4.7400  /
      data      fgi(  7, 39)/   -3.2200  /
      data      fgi(  8, 39)/    0.920000/
      data      fgi(  9, 39)/    0.073000/
      data      fgi( 10, 39)/    0.007000/
      data      fgi( 11, 39)/    2.5700  /
      data      fgi( 12, 39)/    2.4100  /
      data      fgi( 13, 39)/    5.2200  /
      data      fgi( 14, 39)/    2.846438/
      data      fgi( 15, 39)/    2.231965/
      data      fgi( 16, 39)/    1.250107/
      data      fgi( 17, 39)/    1.607476/
      data      fgi( 18, 39)/    1.833888/
      data      fgi( 19, 39)/    1.011633/
      data      fgi( 20, 39)/    3.834798/
      data      fgi( 21, 39)/    2.562213/
      data      fgi( 22, 39)/    1.942611/
      data      fgi( 23, 39)/    1.968347/
      data      fgi( 24, 39)/    1.458400/
!
!                    Data for Element  40     
!
      data        isoki( 40)/    1       /
      data        nbfai( 40)/    9       /
      data      zcoreai( 40)/    4.000000/
      data        zetai( 40)/    1.337000/
      data   zetadi(  1, 40)/    3.639000/
      data   zetadi(  2, 40)/    1.804000/
      data  zetawti(  1, 40)/    0.399240/
      data  zetawti(  2, 40)/    0.713810/
      data   betaai(  1, 40)/   -1.0000  /
      data   betaai(  2, 40)/   -1.0000  /
      data   betaai(  3, 40)/  -10.0000  /
      data      fgi(  1, 40)/   -7.1600  /
      data      fgi(  2, 40)/   -6.4000  /
      data      fgi(  3, 40)/   -4.4300  /
      data      fgi(  4, 40)/   -3.9200  /
      data      fgi(  5, 40)/   -8.0800  /
      data      fgi(  6, 40)/   -5.7900  /
      data      fgi(  7, 40)/   -3.9500  /
      data      fgi(  8, 40)/    0.884000/
      data      fgi(  9, 40)/    0.108000/
      data      fgi( 10, 40)/    0.008000/
      data      fgi( 11, 40)/    3.0600  /
      data      fgi( 12, 40)/    2.9300  /
      data      fgi( 13, 40)/    6.1200  /
      data      fgi( 14, 40)/    2.984853/
      data      fgi( 15, 40)/    2.340500/
      data      fgi( 16, 40)/    1.064414/
      data      fgi( 17, 40)/    1.368699/
      data      fgi( 18, 40)/    1.759978/
      data      fgi( 19, 40)/    0.861363/
      data      fgi( 20, 40)/    4.369791/
      data      fgi( 21, 40)/    2.919668/
      data      fgi( 22, 40)/    1.777103/
      data      fgi( 23, 40)/    1.890505/
      data      fgi( 24, 40)/    1.315937/
!
!                    Data for Element  41     
!
      data        isoki( 41)/    1       /
      data        nbfai( 41)/    9       /
      data      zcoreai( 41)/    5.000000/
      data        zetai( 41)/    1.393000/
      data   zetadi(  1, 41)/    3.774100/
      data   zetadi(  2, 41)/    1.925100/
      data  zetawti(  1, 41)/    0.447750/
      data  zetawti(  2, 41)/    0.663010/
      data   betaai(  1, 41)/   -1.0000  /
      data   betaai(  2, 41)/   -1.0000  /
      data   betaai(  3, 41)/  -13.0000  /
      data      fgi(  1, 41)/   -7.6900  /
      data      fgi(  2, 41)/   -6.8400  /
      data      fgi(  3, 41)/   -4.6700  /
      data      fgi(  4, 41)/   -4.1500  /
      data      fgi(  5, 41)/   -9.5000  /
      data      fgi(  6, 41)/   -6.8200  /
      data      fgi(  7, 41)/   -4.7000  /
      data      fgi(  8, 41)/    0.806000/
      data      fgi(  9, 41)/    0.187000/
      data      fgi( 10, 41)/    0.007000/
      data      fgi( 11, 41)/    3.5200  /
      data      fgi( 12, 41)/    3.4500  /
      data      fgi( 13, 41)/    7.0100  /
      data      fgi( 14, 41)/    3.109873/
      data      fgi( 15, 41)/    2.438532/
      data      fgi( 16, 41)/    0.950646/
      data      fgi( 17, 41)/    1.222408/
      data      fgi( 18, 41)/    1.717580/
      data      fgi( 19, 41)/    0.769298/
      data      fgi( 20, 41)/    4.808812/
      data      fgi( 21, 41)/    3.212999/
      data      fgi( 22, 41)/    1.673072/
      data      fgi( 23, 41)/    1.842142/
      data      fgi( 24, 41)/    1.226002/
!
!                    Data for Element  42     
!
      data        isoki( 42)/    1       /
      data        nbfai( 42)/    9       /
      data      zcoreai( 42)/    6.000000/
      data        zetai( 42)/    1.440000/
      data   zetadi(  1, 42)/    3.954000/
      data   zetadi(  2, 42)/    2.047000/
      data  zetawti(  1, 42)/    0.480130/
      data  zetawti(  2, 42)/    0.628920/
      data   betaai(  1, 42)/   -1.0000  /
      data   betaai(  2, 42)/   -1.0000  /
      data   betaai(  3, 42)/  -15.0000  /
      data      fgi(  1, 42)/   -8.1300  /
      data      fgi(  2, 42)/   -7.1800  /
      data      fgi(  3, 42)/   -4.8700  /
      data      fgi(  4, 42)/   -4.2900  /
      data      fgi(  5, 42)/  -10.8500  /
      data      fgi(  6, 42)/   -7.8200  /
      data      fgi(  7, 42)/   -5.4600  /
      data      fgi(  8, 42)/    0.616000/
      data      fgi(  9, 42)/    0.382000/
      data      fgi( 10, 42)/    0.002000/
      data      fgi( 11, 42)/    3.9800  /
      data      fgi( 12, 42)/    3.9900  /
      data      fgi( 13, 42)/    7.8800  /
      data      fgi( 14, 42)/    3.214801/
      data      fgi( 15, 42)/    2.520808/
      data      fgi( 16, 42)/    0.872601/
      data      fgi( 17, 42)/    1.122052/
      data      fgi( 18, 42)/    1.688591/
      data      fgi( 19, 42)/    0.706140/
      data      fgi( 20, 42)/    5.170238/
      data      fgi( 21, 42)/    3.454485/
      data      fgi( 22, 42)/    1.599116/
      data      fgi( 23, 42)/    1.807038/
      data      fgi( 24, 42)/    1.162179/
!
!                    Data for Element  43     
!
      data        isoki( 43)/    1       /
      data        nbfai( 43)/    9       /
      data      zcoreai( 43)/    7.000000/
      data        zetai( 43)/    1.461000/
      data   zetadi(  1, 43)/    4.124000/
      data   zetadi(  2, 43)/    2.155000/
      data  zetawti(  1, 43)/    0.512410/
      data  zetawti(  2, 43)/    0.595420/
      data   betaai(  1, 43)/   -1.0000  /
      data   betaai(  2, 43)/   -1.0000  /
      data   betaai(  3, 43)/  -17.0000  /
      data      fgi(  1, 43)/   -8.4900  /
      data      fgi(  2, 43)/   -7.4200  /
      data      fgi(  3, 43)/   -5.0000  /
      data      fgi(  4, 43)/   -4.3400  /
      data      fgi(  5, 43)/  -12.1300  /
      data      fgi(  6, 43)/   -8.8000  /
      data      fgi(  7, 43)/   -6.2100  /
      data      fgi(  8, 43)/    0.242000/
      data      fgi(  9, 43)/    0.736000/
      data      fgi( 10, 43)/    0.022000/
      data      fgi( 11, 43)/    4.4000  /
      data      fgi( 12, 43)/    4.5400  /
      data      fgi( 13, 43)/    8.7300  /
      data      fgi( 14, 43)/    3.261683/
      data      fgi( 15, 43)/    2.557570/
      data      fgi( 16, 43)/    0.792260/
      data      fgi( 17, 43)/    1.018744/
      data      fgi( 18, 43)/    1.635130/
      data      fgi( 19, 43)/    0.641126/
      data      fgi( 20, 43)/    5.431609/
      data      fgi( 21, 43)/    3.629120/
      data      fgi( 22, 43)/    1.507582/
      data      fgi( 23, 43)/    1.744914/
      data      fgi( 24, 43)/    1.087087/
!
!                    Data for Element  44     
!
      data        isoki( 44)/    1       /
      data        nbfai( 44)/    9       /
      data      zcoreai( 44)/    8.000000/
      data        zetai( 44)/    1.470000/
      data   zetadi(  1, 44)/    4.259000/
      data   zetadi(  2, 44)/    2.094000/
      data  zetawti(  1, 44)/    0.534230/
      data  zetawti(  2, 44)/    0.592710/
      data   betaai(  1, 44)/   -2.0000  /
      data   betaai(  2, 44)/   -2.0000  /
      data   betaai(  3, 44)/  -20.0000  /
      data      fgi(  1, 44)/   -8.7500  /
      data      fgi(  2, 44)/   -7.5500  /
      data      fgi(  3, 44)/   -5.0700  /
      data      fgi(  4, 44)/   -4.3000  /
      data      fgi(  5, 44)/  -13.3800  /
      data      fgi(  6, 44)/   -9.7700  /
      data      fgi(  7, 44)/   -6.9800  /
      data      fgi(  8, 44)/    0.016000/
      data      fgi(  9, 44)/    0.717000/
      data      fgi( 10, 44)/    0.266000/
      data      fgi( 11, 44)/    4.8100  /
      data      fgi( 12, 44)/    5.1200  /
      data      fgi( 13, 44)/    9.5500  /
      data      fgi( 14, 44)/    2.479677/
      data      fgi( 15, 44)/    1.120070/
      data      fgi( 16, 44)/    0.990011/
      data      fgi( 17, 44)/    0.330045/
      data      fgi( 18, 44)/    1.440073/
      data      fgi( 19, 44)/    0.180025/
      data      fgi( 20, 44)/    6.170305/
      data      fgi( 21, 44)/    4.630177/
      data      fgi( 22, 44)/    1.431276/
      data      fgi( 23, 44)/    1.687971/
      data      fgi( 24, 44)/    1.025590/
!
!                    Data for Element  45     
!
      data        isoki( 45)/    1       /
      data        nbfai( 45)/    9       /
      data      zcoreai( 45)/    9.000000/
      data        zetai( 45)/    1.482000/
      data   zetadi(  1, 45)/    4.485000/
      data   zetadi(  2, 45)/    2.217000/
      data  zetawti(  1, 45)/    0.544090/
      data  zetawti(  2, 45)/    0.581370/
      data   betaai(  1, 45)/   -1.0000  /
      data   betaai(  2, 45)/   -1.0000  /
      data   betaai(  3, 45)/  -21.0000  /
      data      fgi(  1, 45)/   -8.9300  /
      data      fgi(  2, 45)/   -7.5800  /
      data      fgi(  3, 45)/   -5.0700  /
      data      fgi(  4, 45)/   -4.1700  /
      data      fgi(  5, 45)/  -14.5400  /
      data      fgi(  6, 45)/  -10.7000  /
      data      fgi(  7, 45)/   -7.7500  /
      data      fgi(  8, 45)/    0.003000/
      data      fgi(  9, 45)/    0.319000/
      data      fgi( 10, 45)/    0.678000/
      data      fgi( 11, 45)/    5.1900  /
      data      fgi( 12, 45)/    5.7100  /
      data      fgi( 13, 45)/   10.3600  /
      data      fgi( 14, 45)/    3.308566/
      data      fgi( 15, 45)/    2.594332/
      data      fgi( 16, 45)/    0.650041/
      data      fgi( 17, 45)/    0.835868/
      data      fgi( 18, 45)/    1.518541/
      data      fgi( 19, 45)/    0.526037/
      data      fgi( 20, 45)/    5.864504/
      data      fgi( 21, 45)/    3.918358/
      data      fgi( 22, 45)/    1.328619/
      data      fgi( 23, 45)/    1.608894/
      data      fgi( 24, 45)/    0.943435/
!
!                    Data for Element  46     
!
      data        isoki( 46)/    0       /
      data        nbfai( 46)/    9       /
      data      zcoreai( 46)/   10.000000/
      data        zetai( 46)/    1.497000/
      data   zetadi(  1, 46)/    4.723000/
      data   zetadi(  2, 46)/    2.342000/
      data  zetawti(  1, 46)/    0.549020/
      data  zetawti(  2, 46)/    0.575610/
      data   betaai(  1, 46)/   -1.0000  /
      data   betaai(  2, 46)/   -1.0000  /
      data   betaai(  3, 46)/  -21.5000  /
      data      fgi(  1, 46)/   -9.0200  /
      data      fgi(  2, 46)/   -7.5000  /
      data      fgi(  3, 46)/   -5.0200  /
      data      fgi(  4, 46)/   -3.9400  /
      data      fgi(  5, 46)/  -15.6500  /
      data      fgi(  6, 46)/  -11.6200  /
      data      fgi(  7, 46)/   -8.5300  /
      data      fgi(  8, 46)/    0.007000/
      data      fgi(  9, 46)/    0.121000/
      data      fgi( 10, 46)/    0.872000/
      data      fgi( 12, 46)/    6.3300  /
      data      fgi( 13, 46)/   11.1500  /
      data      fgi( 14, 46)/    3.342053/
      data      fgi( 15, 46)/    2.620590/
      data      fgi( 16, 46)/    0.583298/
      data      fgi( 17, 46)/    0.750046/
      data      fgi( 18, 46)/    1.460972/
      data      fgi( 19, 46)/    0.472026/
      data      fgi( 20, 46)/    6.121791/
      data      fgi( 21, 46)/    4.090263/
      data      fgi( 22, 46)/    1.240469/
      data      fgi( 23, 46)/    1.540221/
      data      fgi( 24, 46)/    0.873138/
!
!                    Data for Element  47     
!
      data        isoki( 47)/    1       /
      data        nbfai( 47)/    9       /
      data      zcoreai( 47)/   11.000000/
      data        zetai( 47)/    1.507000/
      data   zetadi(  1, 47)/    4.988960/
      data   zetadi(  2, 47)/    2.583740/
      data  zetawti(  1, 47)/    0.557630/
      data  zetawti(  2, 47)/    0.553590/
      data   betaai(  1, 47)/   -2.5000  /
      data   betaai(  2, 47)/   -2.5000  /
      data   betaai(  3, 47)/  -30.0000  /
      data      fgi(  2, 47)/   -7.3100  /
      data      fgi(  4, 47)/   -3.6200  /
      data      fgi(  6, 47)/  -12.5100  /
      data      fgi(  9, 47)/    1.000000/
      data      fgi( 11, 47)/    6.0000  /
      data      fgi( 12, 47)/    6.9600  /
      data      fgi( 13, 47)/    9.5500  /
      data      fgi( 14, 47)/    3.027900/
      data      fgi( 15, 47)/    2.374200/
      data      fgi( 16, 47)/    0.386800/
      data      fgi( 17, 47)/    0.508300/
      data      fgi( 18, 47)/    1.024000/
      data      fgi( 19, 47)/    0.308600/
      data      fgi( 20, 47)/    5.965000/
      data      fgi( 21, 47)/    3.828500/
!
!                    Data for Element  48     
!
      data        isoki( 48)/    0       /
      data        nbfai( 48)/    4       /
      data      zcoreai( 48)/    2.000000/
      data        zetai( 48)/    1.699000/
      data   betaai(  1, 48)/   -1.0000  /
      data   betaai(  2, 48)/   -1.0000  /
      data      fgi(  1, 48)/   -8.9400  /
      data      fgi(  2, 48)/   -8.9400  /
      data      fgi(  3, 48)/   -4.7500  /
      data      fgi(  4, 48)/   -4.7500  /
      data      fgi(  8, 48)/    1.000000/
      data      fgi( 14, 48)/    3.793018/
      data      fgi( 15, 48)/    2.974203/
!
!                    Data for Element  49     
!
      data        isoki( 49)/    0       /
      data        nbfai( 49)/    4       /
      data      zcoreai( 49)/    3.000000/
      data   betaai(  1, 49)/ -100.0000  /
      data   betaai(  2, 49)/ -100.0000  /
      data      fgi(  1, 49)/ -100.0000  /
      data      fgi(  3, 49)/ -100.0000  /
      data      fgi(  8, 49)/    1.000000/
!
!                    Data for Element  50     
!
      data        isoki( 50)/    0       /
      data        nbfai( 50)/    4       /
      data      zcoreai( 50)/    4.000000/
      data   betaai(  1, 50)/ -100.0000  /
      data   betaai(  2, 50)/ -100.0000  /
      data      fgi(  1, 50)/ -100.0000  /
      data      fgi(  3, 50)/ -100.0000  /
      data      fgi(  8, 50)/    1.000000/
!
!                    Data for Element  51     
!
      data        isoki( 51)/    0       /
      data        nbfai( 51)/    4       /
      data      zcoreai( 51)/    5.000000/
      data   betaai(  1, 51)/ -100.0000  /
      data   betaai(  2, 51)/ -100.0000  /
      data      fgi(  1, 51)/ -100.0000  /
      data      fgi(  3, 51)/ -100.0000  /
      data      fgi(  8, 51)/    1.000000/
!
!                    Data for Element  52     
!
      data        isoki( 52)/    0       /
      data        nbfai( 52)/    4       /
      data      zcoreai( 52)/    6.000000/
      data   betaai(  1, 52)/ -100.0000  /
      data   betaai(  2, 52)/ -100.0000  /
      data      fgi(  1, 52)/ -100.0000  /
      data      fgi(  3, 52)/ -100.0000  /
      data      fgi(  8, 52)/    1.000000/
!
!                    Data for Element  53     
!
      data        isoki( 53)/    1       /
      data        nbfai( 53)/    4       /
      data      zcoreai( 53)/    7.000000/
      data        zetai( 53)/    3.341000/
      data   betaai(  1, 53)/   -8.0000  /
      data   betaai(  2, 53)/   -8.0000  /
      data      fgi(  1, 53)/  -20.8400  /
      data      fgi(  3, 53)/  -11.2100  /
      data      fgi(  8, 53)/    1.000000/
      data      fgi( 11, 53)/    8.1500  /
      data      fgi( 12, 53)/    8.1500  /
      data      fgi( 14, 53)/    4.886521/
      data      fgi( 15, 53)/    3.842482/
!
!                    Data for Element  54     
!
      data        isoki( 54)/    0       /
      data        nbfai( 54)/    4       /
      data      zcoreai( 54)/    8.000000/
      data   betaai(  1, 54)/ -100.0000  /
      data   betaai(  2, 54)/ -100.0000  /
      data      fgi(  1, 54)/ -100.0000  /
      data      fgi(  3, 54)/ -100.0000  /
      data      fgi(  8, 54)/    1.000000/
!
!                    Data for Element  55     
!
      data        isoki( 55)/    0       /
      data        nbfai( 55)/    4       /
      data      zcoreai( 55)/    1.000000/
      data   betaai(  1, 55)/ -100.0000  /
      data   betaai(  2, 55)/ -100.0000  /
      data      fgi(  1, 55)/ -100.0000  /
      data      fgi(  3, 55)/ -100.0000  /
      data      fgi(  8, 55)/    1.000000/
!
!                    Data for Element  56     
!
      data        isoki( 56)/    0       /
      data        nbfai( 56)/    9       /
      data      zcoreai( 56)/    2.000000/
      data   betaai(  1, 56)/ -100.0000  /
      data   betaai(  2, 56)/ -100.0000  /
      data      fgi(  1, 56)/ -100.0000  /
      data      fgi(  3, 56)/ -100.0000  /
      data      fgi(  8, 56)/    1.000000/
!
!                    Data for Element  57     
!
      data        isoki( 57)/    0       /
      data        nbfai( 57)/    9       /
      data      zcoreai( 57)/    3.000000/
      data   betaai(  1, 57)/ -100.0000  /
      data   betaai(  2, 57)/ -100.0000  /
      data      fgi(  1, 57)/ -100.0000  /
      data      fgi(  3, 57)/ -100.0000  /
      data      fgi(  8, 57)/    1.000000/
!
!                    Data for Element  58     
!
      data        isoki( 58)/    0       /
      data        nbfai( 58)/    9       /
      data      zcoreai( 58)/    4.000000/
      data   betaai(  1, 58)/ -100.0000  /
      data   betaai(  2, 58)/ -100.0000  /
      data      fgi(  1, 58)/ -100.0000  /
      data      fgi(  3, 58)/ -100.0000  /
      data      fgi(  8, 58)/    1.000000/
!
!                    Data for Element  59     
!
      data        isoki( 59)/    0       /
      data        nbfai( 59)/    9       /
      data      zcoreai( 59)/    5.000000/
      data   betaai(  1, 59)/ -100.0000  /
      data   betaai(  2, 59)/ -100.0000  /
      data      fgi(  1, 59)/ -100.0000  /
      data      fgi(  3, 59)/ -100.0000  /
      data      fgi(  8, 59)/    1.000000/
!
!                    Data for Element  60     
!
      data        isoki( 60)/    0       /
      data        nbfai( 60)/    9       /
      data      zcoreai( 60)/    6.000000/
      data   betaai(  1, 60)/ -100.0000  /
      data   betaai(  2, 60)/ -100.0000  /
      data      fgi(  1, 60)/ -100.0000  /
      data      fgi(  3, 60)/ -100.0000  /
      data      fgi(  8, 60)/    1.000000/
!
!                    Data for Element  61     
!
      data        isoki( 61)/    0       /
      data        nbfai( 61)/    9       /
      data      zcoreai( 61)/    7.000000/
      data   betaai(  1, 61)/ -100.0000  /
      data   betaai(  2, 61)/ -100.0000  /
      data      fgi(  1, 61)/ -100.0000  /
      data      fgi(  3, 61)/ -100.0000  /
      data      fgi(  8, 61)/    1.000000/
!
!                    Data for Element  62     
!
      data        isoki( 62)/    0       /
      data        nbfai( 62)/    9       /
      data      zcoreai( 62)/    8.000000/
      data   betaai(  1, 62)/ -100.0000  /
      data   betaai(  2, 62)/ -100.0000  /
      data      fgi(  1, 62)/ -100.0000  /
      data      fgi(  3, 62)/ -100.0000  /
      data      fgi(  8, 62)/    1.000000/
!
!                    Data for Element  63     
!
      data        isoki( 63)/    0       /
      data        nbfai( 63)/    9       /
      data      zcoreai( 63)/    9.000000/
      data   betaai(  1, 63)/ -100.0000  /
      data   betaai(  2, 63)/ -100.0000  /
      data      fgi(  1, 63)/ -100.0000  /
      data      fgi(  3, 63)/ -100.0000  /
      data      fgi(  8, 63)/    1.000000/
!
!                    Data for Element  64     
!
      data        isoki( 64)/    0       /
      data        nbfai( 64)/    9       /
      data      zcoreai( 64)/   10.000000/
      data   betaai(  1, 64)/ -100.0000  /
      data   betaai(  2, 64)/ -100.0000  /
      data      fgi(  1, 64)/ -100.0000  /
      data      fgi(  3, 64)/ -100.0000  /
      data      fgi(  8, 64)/    1.000000/
!
!                    Data for Element  65     
!
      data        isoki( 65)/    0       /
      data        nbfai( 65)/    9       /
      data      zcoreai( 65)/   11.000000/
      data   betaai(  1, 65)/ -100.0000  /
      data   betaai(  2, 65)/ -100.0000  /
      data      fgi(  1, 65)/ -100.0000  /
      data      fgi(  3, 65)/ -100.0000  /
      data      fgi(  8, 65)/    1.000000/
!
!                    Data for Element  66     
!
      data        isoki( 66)/    0       /
      data        nbfai( 66)/    9       /
      data      zcoreai( 66)/   12.000000/
      data   betaai(  1, 66)/ -100.0000  /
      data   betaai(  2, 66)/ -100.0000  /
      data      fgi(  1, 66)/ -100.0000  /
      data      fgi(  3, 66)/ -100.0000  /
      data      fgi(  8, 66)/    1.000000/
!
!                    Data for Element  67     
!
      data        isoki( 67)/    0       /
      data        nbfai( 67)/    9       /
      data      zcoreai( 67)/   13.000000/
      data   betaai(  1, 67)/ -100.0000  /
      data   betaai(  2, 67)/ -100.0000  /
      data      fgi(  1, 67)/ -100.0000  /
      data      fgi(  3, 67)/ -100.0000  /
      data      fgi(  8, 67)/    1.000000/
!
!                    Data for Element  68     
!
      data        isoki( 68)/    0       /
      data        nbfai( 68)/    9       /
      data      zcoreai( 68)/   14.000000/
      data   betaai(  1, 68)/ -100.0000  /
      data   betaai(  2, 68)/ -100.0000  /
      data      fgi(  1, 68)/ -100.0000  /
      data      fgi(  3, 68)/ -100.0000  /
      data      fgi(  8, 68)/    1.000000/
!
!                    Data for Element  69     
!
      data        isoki( 69)/    0       /
      data        nbfai( 69)/    9       /
      data      zcoreai( 69)/   15.000000/
      data   betaai(  1, 69)/ -100.0000  /
      data   betaai(  2, 69)/ -100.0000  /
      data      fgi(  1, 69)/ -100.0000  /
      data      fgi(  3, 69)/ -100.0000  /
      data      fgi(  8, 69)/    1.000000/
!
!                    Data for Element  70     
!
      data        isoki( 70)/    0       /
      data        nbfai( 70)/    9       /
      data      zcoreai( 70)/   16.000000/
      data   betaai(  1, 70)/ -100.0000  /
      data   betaai(  2, 70)/ -100.0000  /
      data      fgi(  1, 70)/ -100.0000  /
      data      fgi(  3, 70)/ -100.0000  /
      data      fgi(  8, 70)/    1.000000/
!
!                    Data for Element  71     
!
      data        isoki( 71)/    0       /
      data        nbfai( 71)/    9       /
      data      zcoreai( 71)/   17.000000/
      data   betaai(  1, 71)/ -100.0000  /
      data   betaai(  2, 71)/ -100.0000  /
      data      fgi(  1, 71)/ -100.0000  /
      data      fgi(  3, 71)/ -100.0000  /
      data      fgi(  8, 71)/    1.000000/
!
!                    Data for Element  72     
!
      data        isoki( 72)/    0       /
      data        nbfai( 72)/    9       /
      data      zcoreai( 72)/    4.000000/
      data   betaai(  1, 72)/ -100.0000  /
      data   betaai(  2, 72)/ -100.0000  /
      data      fgi(  1, 72)/ -100.0000  /
      data      fgi(  3, 72)/ -100.0000  /
      data      fgi(  8, 72)/    1.000000/
!
!                    Data for Element  73     
!
      data        isoki( 73)/    0       /
      data        nbfai( 73)/    9       /
      data      zcoreai( 73)/    5.000000/
      data   betaai(  1, 73)/ -100.0000  /
      data   betaai(  2, 73)/ -100.0000  /
      data      fgi(  1, 73)/ -100.0000  /
      data      fgi(  3, 73)/ -100.0000  /
      data      fgi(  8, 73)/    1.000000/
!
!                    Data for Element  74     
!
      data        isoki( 74)/    0       /
      data        nbfai( 74)/    9       /
      data      zcoreai( 74)/    6.000000/
      data   betaai(  1, 74)/ -100.0000  /
      data   betaai(  2, 74)/ -100.0000  /
      data      fgi(  1, 74)/ -100.0000  /
      data      fgi(  3, 74)/ -100.0000  /
      data      fgi(  8, 74)/    1.000000/
!
!                    Data for Element  75     
!
      data        isoki( 75)/    0       /
      data        nbfai( 75)/    9       /
      data      zcoreai( 75)/    7.000000/
      data   betaai(  1, 75)/ -100.0000  /
      data   betaai(  2, 75)/ -100.0000  /
      data      fgi(  1, 75)/ -100.0000  /
      data      fgi(  3, 75)/ -100.0000  /
      data      fgi(  8, 75)/    1.000000/
!
!                    Data for Element  76     
!
      data        isoki( 76)/    0       /
      data        nbfai( 76)/    9       /
      data      zcoreai( 76)/    8.000000/
      data   betaai(  1, 76)/ -100.0000  /
      data   betaai(  2, 76)/ -100.0000  /
      data      fgi(  1, 76)/ -100.0000  /
      data      fgi(  3, 76)/ -100.0000  /
      data      fgi(  8, 76)/    1.000000/
!
!                    Data for Element  77     
!
      data        isoki( 77)/    0       /
      data        nbfai( 77)/    9       /
      data      zcoreai( 77)/    9.000000/
      data   betaai(  1, 77)/ -100.0000  /
      data   betaai(  2, 77)/ -100.0000  /
      data      fgi(  1, 77)/ -100.0000  /
      data      fgi(  3, 77)/ -100.0000  /
      data      fgi(  8, 77)/    1.000000/
!
!                    Data for Element  78     
!
      data        isoki( 78)/    0       /
      data        nbfai( 78)/    9       /
      data      zcoreai( 78)/   10.000000/
      data   betaai(  1, 78)/ -100.0000  /
      data   betaai(  2, 78)/ -100.0000  /
      data      fgi(  1, 78)/ -100.0000  /
      data      fgi(  3, 78)/ -100.0000  /
      data      fgi(  8, 78)/    1.000000/
!
!                    Data for Element  79     
!
      data        isoki( 79)/    1       /
      data        nbfai( 79)/    9       /
      data      zcoreai( 79)/   11.000000/
      data        zetai( 79)/    1.986000/
      data   zetadi(  1, 79)/    3.288000/
      data  zetawti(  1, 79)/    1.000000/
      data   betaai(  1, 79)/   -4.0000  /
      data   betaai(  2, 79)/   -4.0000  /
      data   betaai(  3, 79)/  -50.0000  /
      data      fgi(  1, 79)/  -10.1200  /
      data      fgi(  2, 79)/   -9.2260  /
      data      fgi(  3, 79)/   -5.4320  /
      data      fgi(  4, 79)/   -4.3580  /
      data      fgi(  5, 79)/  -14.6920  /
      data      fgi(  6, 79)/  -12.0170  /
      data      fgi(  9, 79)/    1.000000/
      data      fgi( 11, 79)/    7.3640  /
      data      fgi( 12, 79)/    7.3660  /
      data      fgi( 13, 79)/    9.1490  /
      data      fgi( 14, 79)/    3.456875/
      data      fgi( 15, 79)/    2.749111/
      data      fgi( 16, 79)/    1.011275/
      data      fgi( 17, 79)/    1.280173/
      data      fgi( 18, 79)/    1.942687/
      data      fgi( 19, 79)/    0.827134/
      data      fgi( 20, 79)/    3.569660/
      data      fgi( 21, 79)/    3.530346/
!
!                    Data for Element  80     
!
      data        isoki( 80)/    0       /
      data        nbfai( 80)/    4       /
      data      zcoreai( 80)/    2.000000/
      data   betaai(  1, 80)/ -100.0000  /
      data   betaai(  2, 80)/ -100.0000  /
      data      fgi(  1, 80)/ -100.0000  /
      data      fgi(  3, 80)/ -100.0000  /
      data      fgi(  8, 80)/    1.000000/

  end module Parameters_for_INDO_C

