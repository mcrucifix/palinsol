      subroutine berger75133f(
     +ia, ib, ic,id, ha, prma, pprma, tseta,
     +aa, a, c,
     +bb, b, d,
     +bf, pf, sa, ddr, 
     +baf, s, d11, 
     +paf, ss, d111)

      implicit double precision(a-h,p-z)                                p7500030
      dimension aa(80),a(80),c(80),cyc1(80),bb(80),b(80),d(80),cyc2(80),p7500040
     *il(5), ih(5), 
     *bf(65000),pf(10000),dfp(10000), sa(10000), ddr(10000), 
     *baf(9640),paf(9640),s(9640),d11(9640),
     *baff(9640),paff(9640),d111(9640),ss(9640)                         p7500060
c                                                                       p7500070
c   pal 7513n                                                           p7500080
c   sortie tableau astron astrophys 1976 paper                          p7500090
c   +0.5 pour arrondir                                                  p7500100
c                                                                       p7500110
c                                                                       p7500120
c   avec classement par ordre decroissant                               p7500130
c   papier pour inqua 1977                                              p7500140
c                                                                       p7500150
!      open(unit=5,file='ber907505_1.dat',status='old')
!      open(unit=7,file='ber90_5.dat',status='new')
 
      intent (in) :: ia, ib, ic,id, ha, prma, pprma, tseta
      intent (in) :: aa, a, c
      intent (in) :: bb, b, d
      intent (in) :: bf, pf, sa, ddr
      intent (inout) :: baf, s, d11, paf, ss, d111

      z=1.0d-08                                                         p7500160
      pi=dacos(-1.0d0)                                                  p7500170
      pir=pi/180.0d0                                                    p7500180
      pirr=pir/3600.0d0                                                 p7500190
      atp=1.0d0                                                         p7500200
      atb=0.1d0                                                         p7500210
      il(1)=1                                                           p7500260
      ih(1)=ib                                                          p7500270
      il(2)=ih(1)+1                                                     p7500280
      ih(2)=2*ib                                                        p7500290
      il(3)=ih(2)+1                                                     p7500300
      ih(3)=2*ib+(ib*(ib-1))/2                                          p7500310
      il(4)=ih(3)+1                                                     p7500320
      ih(4)=ic                                                          p7500330
      il(5)=ih(4)+1                                                     p7500340
      ih(5)=id                                                          p7500350
      do 50 i=1,ia                                                      p7500360
 5000 format(i4,f20.8,f20.7,f20.6)                                      p7500380
      aza=a(i)*pirr                                                     p7500390
      aza=dabs(aza)                                                     p7500400
      if (aza.le.z) go to 40                                            p7500410
      cyc1(i)=pi*2.0d0/aza                                              p7500420
      go to 50                                                          p7500430
   40 cyc1(i)=0.0d0                                                     p7500440
   50 continue                                                          p7500450
      do 51 i=1,ib                                                      p7500460
      aza=b(i)*pirr                                                     p7500490
      aza=dabs(aza)                                                     p7500500
      if (aza.le.z) go to 41                                            p7500510
      cyc2(i)=pi*2.0d0/aza                                              p7500520
      go to 51                                                          p7500530
   41 cyc2(i)=0.0d0                                                     p7500540
   51 continue                                                          p7500550
      call claser(aa,a,c,cyc1,ia)                                       p7500560
      call claser(bb,b,d,cyc2,ib)                                       p7500570
c                                                                       p7500900
c   obliquity                                                           p7500910
c       
      do 54 i=1,ic                                                      p7500930
      baf(i)=bf(i)/pirr                                                    p7502580
      paf(i)=pf(i)/pirr                                                    p7502590
      s(i)=sa(i)/pirr                                                      p7502610
      d11(i)=ddr(i)/pir                                                    p7502620
  54  continue
 5002 format(i4,f15.7,f15.7,f13.7,f13.6,f12.4)                          p7500950
      icc=ic+1                                                          p7500970
      do 55 i=icc,id                                                    p7500980
      paf(i)=pf(i)/pirr                                                    p7502680
      s(i)=sa(i)/pirr                                                      p7502700
      d11(i)=ddr(i)/pir                                                    p7502710
 5003 format(i4,15x,f15.7,f13.7,f13.6,f12.4)                            p7501000
   55 continue                                                          p7501010
c
c                                                                       p7501070
c   mean motion positif phase comprise entre 0 et 360 degres            p7501080
c                                                                       p7501090
      do 300 i=1,id                                                     p7501100
      as=s(i)                                                           p7501110
      ad=d11(i)                                                         p7501120
      if (as.ge.0.0d0) go to 301                                        p7501130
      as=-as                                                            p7501140
      ad=-ad                                                            p7501150
      paf(i)=-paf(i)                                                    p7501160
  301 call subad(ad)                                                    p7501170
      s(i)=as                                                           p7501180
      d11(i)=ad                                                         p7501190
  300 continue                                                          p7501200
      k=0                                                               p7501210
      iz=0                                                              p7501220
      do 200 jj=1,4                                                     p7501230
      icc=il(jj)                                                        p7501240
      iccc=ih(jj)                                                       p7501250
      do 201 i=icc,iccc                                                 p7501260
      if (i.eq.1) go to 203                                             p7501270
      ad=d11(i)                                                         p7501280
      as=s(i)                                                           p7501290
c                                                                       p7501300
c   je regarde si 1 valeur est egale a une valeur precedente d'un nouveap7501310
c       vecteur dont toutes les composantes sont differentes            p7501320
c                                                                       p7501330
      do 202 j=1,k                                                      p7501340
      if (as.ne.ss(j)) go to 202                                        p7501350
      if (ad.ne.d111(j)) go to 202                                      p7501360
      baff(j)=baff(j)+baf(i)                                            p7501370
      paff(j)=paff(j)+paf(i)                                            p7501380
 6101 format(1x,i4,2x,i4,4x,f15.7,4x,f15.7,4x,f13.6,4x,f13.5)           p7501410
      go to 201                                                         p7501420
  202 continue                                                          p7501430
  203 k=k+1                                                             p7501440
      ss(k)=s(i)                                                        p7501450
      paff(k)=paf(i)                                                    p7501460
      baff(k)=baf(i)                                                    p7501470
      d111(k)=d11(i)                                                    p7501480
!      write(6,6101) k,i,baff(k),paff(k),ss(k),d111(k)                   p7501490
  201 continue                                                          p7501500
!      write(6,6010)                                                     p7501510
!      write(6,6010)                                                     p7501520
      ih(jj)=k                                                          p7501530
  200 continue                                                          p7501540
      icc=il(5)                                                         p7501550
      iccc=ih(5)                                                        p7501560
      do 210 i=icc,iccc                                                 p7501570
      ad=d11(i)                                                         p7501580
      as=s(i)                                                           p7501590
      do 212 j=1,k                                                      p7501600
      if (as.ne.ss(j)) go to 212                                        p7501610
      if (ad.ne.d111(j)) go to 212                                      p7501620
      paff(j)=paff(j)+paf(i)                                            p7501630
      iz=0                                                              p7501640
!      write(6,6102) j,iz,paff(j),ss(j),d111(j)                          p7501650
!      write(6,6102) iz,i,paf(i),s(i),d11(i)                             p7501660
 6102 format(1x,i4,2x,i4,19x,4x,f15.7,4x,f13.6,4x,f13.5)                p7501670
      go to 210                                                         p7501680
  212 continue                                                          p7501690
  213 k=k+1                                                             p7501700
      ss(k)=s(i)                                                        p7501710
      paff(k)=paf(i)                                                    p7501720
      d111(k)=d11(i)                                                    p7501730
!      write(6,6102) k,i,paff(k),ss(k),d111(k)                           p7501740
  210 continue                                                          p7501750
      ih(5)=k                                                           p7501760
c                                                                       p7501770
      do 204 i=2,5                                                      p7501780
      il(i)=ih(i-1)+1                                                   p7501790
  204 continue                                                          p7501800
c                                                                       p7501810
      iab=ih(4)                                                         p7501820
      do 215 i=1,iab                                                    p7501830
      s(i)=ss(i)                                                        p7501840
      d11(i)=d111(i)                                                    p7501850
      paf(i)=paff(i)                                                    p7501860
      baf(i)=baff(i)                                                    p7501870
  215 continue                                                          p7501880
c                                                                       p7501890
      iab=ih(4)+1                                                       p7501900
      iabb=ih(5)                                                        p7501910
      do 216 i=iab,iabb                                                 p7501920
      s(i)=ss(i)                                                        p7501930
      d11(i)=d111(i)                                                    p7501940
      paf(i)=paff(i)                                                    p7501950
  216 continue                                                          p7501960
c                                                                       p7501970
c   obliquity and precession print                                      p7501980
c                                                                       p7501990
!      write(6,6005)                                                     p7502000
 6005 format(1h1,///,1x,'table 3  obliquity and precession relative to',p7502010
     */,1x,9x,'mean ecliptic and mean equinox of date',/)               p7502020
!      write(6,6006)                                                     p7502030
 6006 format(4x,'i','   obliquity',3x,'precession',3x,'mean  rate ',4x,2p7502040
     *x,'phase',3x,4x,2x,'period',/)                                    p7502050
!      write(6,6106)                                                     p7502060
 6106 format(11x,'('''')',9x,'('''')',6x,'(''''/year)',6x,'(degree)',7x,p7502070
     *'(years)',/)                                                      p7502080
      k=0                                                               p7502090
      do 57 j=1,4                                                       p7502100
      icc=il(j)                                                         p7502110
      iccc=ih(j)                                                        p7502120
      do 58 i=icc,iccc                                                  p7502130
      ax=baf(i)                                                         p7502140
      ay=paf(i)                                                         p7502150
      ad=d11(i)                                                         p7502160
      as=s(i)                                                           p7502170
   59 if (as.eq.0.0d0) go to 100                                        p7502180
      ac=360.0d0*3600.0d0/as                                            p7502190
      go to 101                                                         p7502200
  100 ac=0.0d0                                                          p7502210
  101 continue                                                          p7502220
      call subxy(k,ax,ay,as,ad,ac)                                      p7502230
   58 continue                                                          p7502240
!      write(6,6010)                                                     p7502250
 6010 format(1x)                                                        p7502260
   57 continue                                                          p7502270
c                                                                       p7502280
c   precession                                                          p7502290
c                                                                       p7502300
!      write(6,6012)                                                     p7502310
 6012 format(8x,'additional terms for precession',/)                    p7502320
      icc=il(5)                                                         p7502330
      iccc=ih(5)                                                        p7502340
      do 63 i=icc,iccc                                                  p7502350
      ay=paf(i)                                                         p7502360
      as=s(i)                                                           p7502370
      ad=d11(i)                                                         p7502380
   61 if (as.eq.0.0d0) go to 102                                        p7502390
      ac=360.0d0*3600.0d0/as                                            p7502400
      go to 103                                                         p7502410
  102 ac=0.0d0                                                          p7502420
  103 continue                                                          p7502430
      ayy=dabs(ay)                                                      p7502440
      if (ayy.ge.atp) go to 62                                          p7502450
      go to 63                                                          p7502460
   62 k=k+1                                                             p7502470
!      write(6,6011) k,ay,as,ad,ac                                       p7502480
 6011 format(1x,i4,16x,f8.2,4x,f10.6,4x,f10.5,4x,f10.0)                 p7502490
   63 continue                                                          p7502500
c                                                                       p7502510
c   ces cartes jusqu a la prochaine rouge sont a considerer c rouges    p7502520
c                                                                       p7502530
c   il est necessaire de mettre le classement ici                       p7502540
c    car le groupement des divers termes avait ete fait par les sommatiop7502550
c                                                                       p7502560
c   il faut sauver les  s et d11                                        p7502570
c                                                                       p7502580
      iac=ih(5)                                                         p7502590
      do 900 i=1,iac                                                    p7502600
      ss(i)=s(i)                                                        p7502610
      d111(i)=d11(i)                                                    p7502620
  900 continue                                                          p7502630
      iab=ih(4)                                                         p7502640
      call clase(baf,s,d11,iab)                                         p7502650
!      write(6,6500)                                                     p7502660
 6500 format(1h1,///,1x,'table    obliquity relative to mean ecliptic ofp7502670
     * date',/)                                                         p7502680
 6550 format(1h1,///,1x,'table   precession relative to mean ecliptic ofp7502690
     * date',/)                                                         p7502700
!      write(6,6501)                                                     p7502710
 6501 format(4x,'i','   obliquity',4x,'mean  rate',8x,'phase',9x,'periodp7502720
     *',/)                                                              p7502730
 6551 format(4x,'i','  precession',4x,'mean  rate',8x,'phase',9x,'periodp7502740
     *',/)                                                              p7502750
!      write(6,6502)                                                     p7502760
 6502 format(9x,'('''')',9x,'(''''/year)',6x,'(degree)',7x,'(years)',/) p7502770
      m=0                                                               p7502780
      mm=0                                                              p7502790
      do 500 i=1,iab                                                    p7502800
!      write (*,*) 'i=',iab
      wxa=dabs(baf(i))                                                  p7502810
      if(wxa.le.atb) go to 500                                          p7502820
      if(s(i).eq.0.0d0) go to 501                                       p7502830
      ac=360.0d0*3600.0d0/s(i)                                          p7502840
      go to 502                                                         p7502850
  501 ac=0.0d0                                                          p7502860
!  502 write(6,6505) i,baf(i),s(i),d11(i),ac                             p7502870
  502 continue
 6505 format(1x,i4,4x,f8.2,4x,f10.6,4x,f10.4,4x,f10.0)                  p7502880
 7000 format(i5,2x,f13.7,2x,f10.6,2x,f10.4,2x,f10.0)                    p7502900
      m=m+1                                                             p7502910
      mm=mm+1                                                           p7502920
      if(m.eq.10) go to 505                                             p7502930
      go to 500                                                         p7502940
!  505 write(6,6010)                                                     p7502950
  505 continue
      m=0                                                               p7502960
  500 continue                                                          p7502970
      act=atb                                                           p7502980
      xct=ha                                                            p7502990
!      write(*,*) 'calling control'
      call cntrol(baf,s,d11,iab,xct,act)                                p7503000
      iac=ih(5)                                                         p7503010
      call clase(paf,ss,d111,iac)                                       p7503020
!      write(6,6550)                                                     p7503030
!      write(6,6551)                                                     p7503040
!      write(6,6502)                                                     p7503050
      m=0                                                               p7503060
      mn=0                                                              p7503070
      do 550 i=1,iac                                                    p7503080
      wxa=dabs(paf(i))                                                  p7503090
      if(wxa.le.atp) go to 550                                          p7503100
      if(ss(i).eq.0.0d0) go to 551                                      p7503110
      ac=360.0d0*3600.0d0/ss(i)                                         p7503120
      go to 552                                                         p7503130
  551 ac=0.0d0                                                          p7503140
!  552 write(6,6505) i,paf(i),ss(i),d111(i),ac                           p7503150
  552 continue
!      write(*,*) 'm=',m 
      m=m+1                                                             p7503170
      mn=mn+1                                                           p7503180
      if(m.eq.10) go to 555                                             p7503190
      go to 550                                                         p7503200
  555 continue
!  555 write(6,6010)                                                     p7503210
      m=0                                                               p7503220
  550 continue                                                          p7503230
c                                                                       p7503240
!      write(6,6666) mm,mn                                               p7503250
 6666 format(//,1x,2i6)                                                 p7503260
 7777 format(2i6)                                                       p7503280
      end                                                               p7503300
      subroutine subad(ad)                                              p7503310
      implicit double precision(a-h,p-z)                                p7503320
      cennn=360.0d0                                                     p7503330
   59 if (ad.ge.cennn) go to 57                                         p7503340
      go to 58                                                          p7503350
   57 ad=ad-cennn                                                       p7503360
      go to 59                                                          p7503370
   58 if (ad.lt.0.0d0) go to 61                                         p7503380
      go to 62                                                          p7503390
   61 ad=ad+cennn                                                       p7503400
      go to 58                                                          p7503410
   62 return                                                            p7503420
      end                                                               p7503430
      subroutine subxy(k,ax,ay,as,ad,ac)                                p7503440
      implicit double precision(a-h,p-z)                                p7503450
      atb=0.1d0                                                         p7503460
      atp=1.0d0                                                         p7503470
      ic1=1                                                             p7503480
      ic2=1                                                             p7503490
      axx=dabs(ax)                                                      p7503500
      ayy=dabs(ay)                                                      p7503510
      if (axx.le.atb) ic1=0                                             p7503520
      if (ayy.le.atp) ic2=0                                             p7503530
      if ((ic1.eq.0).and.(ic2.eq.0)) go to 1                            p7503540
      if (ic1.eq.0) go to 2                                             p7503550
      if (ic2.eq.0) go to 3                                             p7503560
      k=k+1                                                             p7503570
!      write (6,6001) k,ax,ay,as,ad,ac                                   p7503580
 6001 format(1x,i4,4x,f8.2,4x,f8.2,4x,f10.6,4x,f10.5,4x,f10.0)          p7503590
      return                                                            p7503600
    1 return                                                            p7503610
    2 k=k+1                                                             p7503620
!      write(6,6002) k,ay,as,ad,ac                                       p7503630
 6002 format(1x,i4,12x,4x,f8.2,4x,f10.6,4x,f10.5,4x,f10.0)              p7503640
      return                                                            p7503650
    3 k=k+1                                                             p7503660
!      write(6,6003) k,ax,as,ad,ac                                       p7503670
 6003 format(1x,i4,4x,f8.2,12x,4x,f10.6,4x,f10.5,4x,f10.0)              p7503680
      return                                                            p7503690
      end                                                               p7503700
      subroutine claser(xx,x,y,z,n)                                     p7503710
      implicit double precision(a-h,p-z)                                p7503720
      dimension xx(80),x(80),y(80),z(80)                                p7503730
   50 continue                                                          p7503740
      m=0                                                               p7503750
      do 1 i=2,n                                                        p7503760
      wa=dabs(xx(i))                                                    p7503770
      wb=dabs(xx(i-1))                                                  p7503780
      if(wa.gt.wb) go to 2                                              p7503790
      go to 1                                                           p7503800
    2 aa=xx(i)                                                          p7503810
      ab=x(i)                                                           p7503820
      ac=y(i)                                                           p7503830
      ad=z(i)                                                           p7503840
      xx(i)=xx(i-1)                                                     p7503850
      x(i)=x(i-1)                                                       p7503860
      y(i)=y(i-1)                                                       p7503870
      z(i)=z(i-1)                                                       p7503880
      xx(i-1)=aa                                                        p7503890
      x(i-1)=ab                                                         p7503900
      y(i-1)=ac                                                         p7503910
      z(i-1)=ad                                                         p7503920
      m=m+1                                                             p7503930
    1 continue                                                          p7503940
      if(m.eq.0) go to 51                                               p7503950
      go to 50                                                          p7503960
   51 return                                                            p7503970
      end                                                               p7503980
      subroutine clase(xx,x,y,n)                                        p7503990
      implicit double precision(a-h,p-z)                                p7504000
      dimension xx(9640),x(9640),y(9640)                                p7504010
   50 continue                                                          p7504020
      m=0                                                               p7504030
      do 1 i=2,n                                                        p7504040
      wa=dabs(xx(i))                                                    p7504050
      wb=dabs(xx(i-1))                                                  p7504060
      if(wa.gt.wb) go to 2                                              p7504070
      go to 1                                                           p7504080
    2 aa=xx(i)                                                          p7504090
      ab=x(i)                                                           p7504100
      ac=y(i)                                                           p7504110
      xx(i)=xx(i-1)                                                     p7504120
      x(i)=x(i-1)                                                       p7504130
      y(i)=y(i-1)                                                       p7504140
      xx(i-1)=aa                                                        p7504150
      x(i-1)=ab                                                         p7504160
      y(i-1)=ac                                                         p7504170
      m=m+1                                                             p7504180
    1 continue                                                          p7504190
      if(m.eq.0) go to 51                                               p7504200
      go to 50                                                          p7504210
   51 return                                                            p7504220
      end                                                               p7504230
      subroutine cntrol(xx,x,y,n,xct,act)                               p7504240
      implicit double precision(a-h,p-z)                                p7504250
      dimension xx(9640),x(9640),y(9640)                                p7504260
      pir=dacos(-1.0d0)                                                 p7504270
!      write(6,6001)                                                     p7504280
 6001 format(1h1,///,' control',//)                                     p7504290
      do 2 j=1,11                                                       p7504300
      t=-(j-1)*1000.0d0*100.0d0                                         p7504310
      xctt=xct                                                          p7504320
      do 1 i=1,n                                                        p7504330
      xxa=dabs(xx(i))                                                   p7504340
      if(xxa.gt.act) go to 3                                            p7504350
      go to 1                                                           p7504360
    3 continue                                                          p7504370
      xa=xx(i)/3600.0d0                                                 p7504380
      xp=(x(i)*t/3600.0d0+y(i))*pir                                     p7504390
      xctt=xctt+xa*dcos(xp)                                             p7504400
    1 continue                                                          p7504410
!      write(6,6000) t,xctt                                              p7504420
 6000 format(1x,f10.0,f15.7)                                            p7504430
    2 continue                                                          p7504440
      return                                                            p7504450
      end                                                               p7504540
