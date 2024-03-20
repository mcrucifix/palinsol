      subroutine berger(
     +  ia, ib, ic, id, 
     +  bea, prea, prega, ala, apoa, 
     +  prma,pprma, ha, hha,  tseta,aa,a,c,bb,b,d,bf,pf,dpf,sa,ddr, 
     +  qaa,qa,qc,
     +  filename, icl)
        implicit double precision (a-h,p-z)                        
        integer, intent(in) :: icl
        character(len=icl), intent(in) :: filename
        dimension aa(80),a(80),c(80),
     *bbt(80),bt(80),dt(80),                                            p7500040
     1ar(80),cr(80),br(500),dr(500),bb(500),b(500),d(500),              p7500050
     *bf(65000),pf(10000),dpf(10000),sa(10000),ddr(10000),              p7500060
     2ca(9),cb(3),cbb(3),zeb(3)                                         p7500070
     *,qaa(12900),qa(12900),qc(12900),ind(12900)                        p7500080
        intent(in)  :: bea,prea, prega, ala, apoa
        intent(out) :: ia,ib,ic,id
        intent(inout) :: prma, ha, hha,  tseta
        intent(out) :: pprma,
     *               aa,a,c,bb,b,d,bf,pf,dpf,sa,ddr,qaa,qa,qc
        common kktilde
c                                                                       p7500090
c       pal7505ee modifie   data a1                                     p7500100
c                                                                       p7500110
c                                                                       p7500120
c  ----
c  comments my mc
c  summary of the constants that we provide
c  - ala is the "newcomb" contant, which berger (or sharaf budnikova?)
c    redined as 'l0'. see more info in R/ls_constants
c  - apoa is P0, and is defined as dl/de^2 
c    if my understading of  the code below is correct, both l0 and p0
c    are defined (and have to be defined for) zero earch eccentricity
c    both ala and apoa have to be defined consistently for a given
c    Earth-Moon distance, etc. 
c    prea is, I am reasnoably confident, the dpsi/dt given in his eq. 66
c             p. 114
c             which will be used to compute psi_bar (still need to
c             this).  This is an initial condition at time t-0. 
c              it is important to get it right for realistic solutions
c              over the latest 2 Myr. But for a speculative solution far
c              back in time this should be sampled 
c   - prega is I is the "psi_0", that is , the value of "psi" at the 
c              time=0 (the reference time). THis again should be sampled
c              for far-back-in-time solutions. 
c              Psi  est la longitude du noeud gamma (de la date), tel
c                   que ramene sur l'eclipitique de reference
c                   voir Fig. 3 de Berger et Loutre 1990 QRes. 
c              it could be sampled anywhere between -pi and pi
c    - bea is the obliquity it time zero that is used to compute prma
c           in degrees
c    - ha is the initial guess proposed for iterative search for
c      the 'h' that finally appears in eq. 24 of Berger and Loutre
c      and is also output
c    - prma is --- I guess - output only and contains the mean
c           precession rate 'k'. But wait, this is this different
c           from the "kbar" that, really, is the one that should be used
c           SO TODO: CHECK THE DIFFERENCE
c                                                                       p7500500
c     pal 7505ee   calcul precession 2d degre                           p7500510
c   calcul des coeff developpement sin i sin omega    obli  precession  p7500520
c   pour paleoclimatologie directement                                  p7500530
c   precession par an julien    trani=1.0d-05                           p7500540
c   obl prec  to  1sec    pour tableaux prog 7513                       p7500550
c                                                                       p7500560
!      write(6,6000)                                                     p7500570
 6000 format(1h1,///,10x,'solution berger-laskar-sharaf   prec degre
     * 2 excen   pal 7505ff',///)                                       p7500590
c                                                                       p7500600
c     constantes                                                        p7500610
c
c   prm .eq. pprm ---->kktilde=1
c   prm .ne. pprm ---->kktilde=0
      kktilde=0
c
      pi=dacos(-1.0d0)                                                  p7500620
      pir=pi/180.0d0                                                    p7500630
      pirr=pir/3600.0d0                                                 p7500640
      pip=2.0d0*pi                                                      p7500650
      un=1.0d0                                                          p7500660
      z=0.0d0                                                           p7500670
      zz=1.0d-08                                                        p7500680
!  open(unit=7,status='new',file='ber907505_1')
c                                                                       p7500690
c constantes du probleme                                                p7500700
c ala apoa constante                                                    p7500710
c bea prea valeurs en 1950.0       refer to mean elements of 1850       p7500720
c                                                                       p7500730
c lecture of the data for e,i,pi,omega
c
***************************************************************************
*c pour une solution bretagnon                                            
c       ia=19                                                           p7500740
c       ibt=15                                                          p7500750
c       open(unit=1,status='old',name='ioep2.data')
c	read(1,*)
c	read(1,*)
c	do 1 i=1,ia
c       read(1,1000)aa(i),a(i),c(i)
c1	continue
c	read(1,*)
c	read(1,*)
c	do 2 i=1,ibt
c       read(1,1000)bbt(i),bt(i),dt(i)
c2	continue
c       kz=5
c1000	format(f12.8,f11.6,f12.6)
c
c  constantes 1950.0 par rapport a 1850.0
c
c     bea=23.4458d0                                                     p7500840
c     prea=50.2686d0                                                    p7500860
c     prega=5025.74d0                                                   p7500880
c
***************************************************************************
***************************************************************************
*c pour une solution laskar
 	open(unit=10,file=filename,status='old')
        read(10,*)
 	read(10,1001)ia
 	do 1 i=1,ia
 	 read(10,1000)k, a(i),ampl,c(i)
 	 aa(i)=ampl*1.0d-08
 1	continue
        read(10,*)
 	read(10,1001)ibt
 	do 2 i=1,ibt
 	 read(10,1000)k, bt(i),ampl,dt(i)
 	 bbt(i)=ampl*1.0d-08
 2	continue
         kz=1                                                           p7500780
 1000	format(i4,f15.6,f15.0,f15.6)
 1001	format(i4)
c
c  constantes 1950.0 par rapport a 2000.0 (Andoyer)
c
c     bea=23.4458d0                                                     p7500840
c      prea=50.2686d0                                                    p7500860
c      prega=-2513.65d0                                                  p7500880
c
c
c  constantes 1950.0 par rapport a 2000.0 (Lieske)
c
!      bea=23.44579d0                                                    p7500840
!      prea=50.27985d0                                                   p7500860
!      prega=-2514.27d0                                                  p7500880
c
c   constantes utlisees dans ber90

!      bea=23.44579d0                                                    p7500840
!      prea=50.273147d0                                                  p7500860
!      prega=-2514.27d0                                                  p7500880

c
c  constantes pour le Devonien

!      bea=23.4458d0                                                     p7500840
!      prea=58.577730                                                    p7500860
!      prega=-2513.65d0                                                  p7500880
***************************************************************************
c   kz  5th term for incl on invar plane                                p7500760
c    ne  3 equations a resoudre                                         p7500770
      ne=3                                                              p7500790
      ks=0
!      ala=54.9066d0                                                     p7500800
      al=ala*pirr                                                       p7500810
!      apoa=17.3919d0                                                    p7500820
      apo=apoa*pirr                                                     p7500830
      be=bea*pir                                                        p7500850
      pre=prea*pirr                                                     p7500870
      preg=prega*pirr                                                   p7500890
c
c
c   origine est a present 1950                                          p7500910
c                                                                       p7500920
***************************************************************************
* bretagnon : t0=1850
c       ct=100.0d0/3600.0d0                                             p7500930
***************************************************************************
***************************************************************************
* laskar : t0=2000
  	ct=-50.0d0/3600.0d0
***************************************************************************
      do 98 iam=1,ia                                                    p7500940
      c(iam)=c(iam)+a(iam)*ct                                           p7500950
   98 continue                                                          p7500960
      do 99 iam=1,ibt                                                   p7500970
      dt(iam)=dt(iam)+bt(iam)*ct                                        p7500980
   99 continue                                                          p7500990
c                                                                       p7501000
c impression des data et transformation en radians                      p7501010
c data berger-bretagnon sin(i/2) bbt bt dt                              p7501020
c data berger-laskar sin(i/2) bbt bt dt                                 p7501020
c                                                                       p7501030
!      write(6,6001)                                                     p7501040
 6001 format(1h ,1x,'solution   laskar 1988',/)                         p7501050
      anumi=2.0d0                                                       p7501060
      call impd(ia,ibt,aa,a,c,bbt,bt,dt,kz,ar,cr,br,dr,pip,pirr,pir,anump7501070
     *i)                                                                p7501080
c                                                                       p7501090
c calcul developpement sin(i) from sin(i/2)                             p7501100
c ibs=nbre suppose de termes                                            p7501110
c ib=nbre termes retenus                                                p7501120
      ibs=500                                                           p7501130
      call trani(ibt,bbt,bt,dt,ibs,bb,b,d,ib)                           p7501140
c                                                                       p7501150
c impressions des donnees transformees et transf. en radians            p7501160
c data berger sin(i)   bb b d                                           p7501170
c                                                                       p7501180
!      write(6,6002)                                                     p7501190
 6002 format(1h1,///,1x,'solution   berger 1975',/)                     p7501200
      anumi=1.0d0                                                       p7501210
      call impd(ia,ib,aa,a,c,bb,b,d,kz,ar,cr,br,dr,pip,pirr,pir,anumi)  p7501220
c                                                                       p7501230
c resolution du systeme en h prm tset                                   p7501240
c                                                                       p7501250
!      write(6,6003)                                                     p7501260
! 6003 format(1h1,///,1x,'calcul iteration h  prm  tset',//)             p7501270
      ic=ib*(ib+1)                                                      p7501280
      id=ic+((ia-1)*ia)/2                                               p7501290
! first guesses for epsilonbar, psibar, et zeta_bar (that he spells
! tseta)

! commented out for being called in R.

      !ha=23.4d0            ! value of epsilon used in eq. 67
      !prma=50.44d0         
!tseta=1.964d0                                                     p7501320
      h=ha*pir                                                          p7501330
      prm=prma*pirr                                                     p7501340
      tset=tseta*pir                                                    p7501350
      icim=0                                                            p7501360
c                                                                       p7501370
c   debut de l'iteration                                                p7501380
  103 continue                                                          p7501390
c                                                                       p7501400
c calcul coeff bf pf dpf pour h prm tset donne                          p7501410
c                                                                       p7501420
      call coef(h,prm,tset,apo,al,icim,ia,aa,a,c,ar,cr,ib,bb,b,d,br,dr, p7501430
     *bf,pf,dpf,sa,ddr,ic,id,pi,pir,pirr,hh)                            p7501440
c                                                                       p7501450
c calcul des 3 fonctions                                                p7501460
c                                                                       p7501470
      t=0.0d0                                                           p7501480
      call fbpf(t,ic,id,bf,pf,sa,ddr,fbf,fpf)                           p7501490
      fbfo=hh+fbf-be                                                    p7501500
      fpfo=tset+fpf-preg                                                p7501510
      call fdpf1(t,id,dpf,sa,ddr,fdpf)                                  p7501520
c
c prn .eq. pprn
c
      if (kktilde.eq.1)then
       fdpfo=prm+fdpf-pre                                               p7501530
         endif
c
c prn .ne. pprn
c
      if (kktilde.eq.0)then
      so=0.0d0
      do 104 i=1,ib
       xw=br(i)+prm
       ai=prm/xw
       x=ai*ai
       bc=bb(i)*bb(i)
       so=so+bc*(1.5d0+0.75d0*x-2.5d0*ai-0.5*ai*(ai-1.0d0)
     $       *dtan(h)*dtan(h))
104   continue
c     tw=3.0d0*apo/al                                                   p7508680
c     do 105 i=1,ia-1                                                   p7508690
c     do 105 k=i+1,ia                                                   p7508710
c     so=so+tw*aa(i)*aa(k)*dcos(cr(i)-cr(k))                            p7508750
c 105 continue                                                          p7508780
      fdpfo=prm*(1.0d0-so)+fdpf-pre                                     p7501530
        endif
c
c                                                                       p7501540
c calcul des derivees des trois fonctions                               p7501550
c                                                                       p7501560
c                                                                       p7501570
c sbfh = derivee de obl mobi par rapport a h                            p7501580
c spfk = derivee de prec mobi par rapport a (k=prm)                     p7501590
c sdpfa = derivee de taux prec mobi par rapport a (a=tset)              p7501600
c                                                                       p7501610
      call dff(h,prm,tset,ib,bb,br,dr,sbfh,sbfk,sbfa,spfh,spfk,spfa,    p7501620
     *sdpfh,sdpfk,sdpfa,ia,aa,ar,cr,apo,al)                             p7501630
c                                                                       p7501640
c calcul du nouveau pas d'iteration                                     p7501650
c ne=3 (equations en h,prm,tset)                                        p7501660
c apres simq  sol. sont dans cb(1)=h cb(2)=prm cb(3)=tset               p7501670
c                                                                       p7501680
      ca(1)=sbfh                                                        p7501690
      ca(4)=sbfk                                                        p7501700
      ca(7)=sbfa                                                        p7501710
      ca(2)=spfh                                                        p7501720
      ca(5)=spfk                                                        p7501730
      ca(8)=spfa                                                        p7501740
      ca(3)=sdpfh                                                       p7501750
      ca(6)=sdpfk                                                       p7501760
      ca(9)=sdpfa                                                       p7501770
      cb(1)=-fbfo                                                       p7501780
      cb(2)=-fpfo                                                       p7501790
      cb(3)=-fdpfo                                                      p7501800
      ca1=sbfh/pirr                                                     p7501810
      ca4=sbfk/pirr                                                     p7501820
      ca7=sbfa/pirr                                                     p7501830
      ca2=spfh/pirr                                                     p7501840
      ca5=spfk/pirr                                                     p7501850
      ca8=spfa/pirr                                                     p7501860
      ca3=sdpfh/pirr                                                    p7501870
      ca6=sdpfk/pirr                                                    p7501880
      ca9=sdpfa/pirr                                                    p7501890
      cb1=-fbfo/pirr                                                    p7501900
      cb2=-fpfo/pirr                                                    p7501910
      cb3=-fdpfo/pirr                                                   p7501920
      hha=hh/pir                                                        p7501930
!      write(6,6004) ha,hha,prma,tseta                                   p7501940
! 6004 format(1h,4f14.6)                                                 p7501950
!      write(6,6005) ca1,ca4,ca7,ca2,ca5,ca8,ca3,ca6,ca9                 p7501960
! 6005 format(1h ,9d14.6)                                                p7501970
!      write(6,6006) cb1,cb2,cb3                                         p7501980
! 6006 format(1h ,3d14.6)                                                p7501990
      call simq(ca,cb,ne,ks)                                            p7502000
      cbb(1)=h                                                          p7502010
      cbb(2)=prm                                                        p7502020
      cbb(3)=tset                                                       p7502030
      cb1=cb(1)/pir                                                     p7502040
      cb2=cb(2)/pirr                                                    p7502050
      cb3=cb(3)/pir                                                     p7502060
!      write(6,6006) cb1,cb2,cb3                                         p7502070
      icon=0                                                            p7502080
      zzz=1.0d-05                                                       p7502090
c                                                                       p7502100
c   precision du calcul                                                 p7502110
c   ha  1.0d-04 degre   prma  1.0d-04 sec   tseta  1.0d-04 degre        p7502120
c                                                                       p7502130
      zeb(1)=zzz*pir                                                    p7502140
      zeb(2)=zzz*pirr                                                   p7502150
      zeb(3)=zzz*pir                                                    p7502160
      do 101 i=1,3                                                      p7502170
      cbi=dabs(cb(i))                                                   p7502180
      cbb(i)=cbb(i)+cb(i)                                               p7502190
      if (cbi.le.zeb(i)) go to 101                                      p7502200
      icon=icon+1                                                       p7502210
  101 continue                                                          p7502220
      h=cbb(1)                                                          p7502230
      ha=h/pir                                                          p7502240
      prm=cbb(2)                                                        p7502250
      prma=prm/pirr                                                     p7502260
      tset=cbb(3)                                                       p7502270
      tseta=tset/pir                                                    p7502280
      call conpr(pprm,h,al,apo,prm,ib,bb,br,ia,aa,cr)                   p7502290
      pprma=pprm/pirr                                                   p7502300
!      write(6,6007) ha,prma,pprma,tseta                                 p7502310
! 6007 format(1h ,4f14.6,//)                                             p7502320
      if (icon.eq.0) go to 102                                          p7502330
      go to 103                                                         p7502340
  102 icim=1                                                            p7502350
      call coef(h,prm,tset,apo,al,icim,ia,aa,a,c,ar,cr,ib,bb,b,d,br,dr, p7502360
     *bf,pf,dpf,sa,ddr,ic,id,pi,pir,pirr,hh)                            p7502370
      hha=hh/pir                                                        p7502380
!      write(6,6900) ha,hha,prma,pprma,tseta                             p7502390
!      write(7,6900) ha,hha,prma,pprma,tseta                             p7502390
! 6900 format(///,1x,5f14.6)                                             p7502400
      ha=hha                                                            p7502410
c   car pour la suite je n utiliserai plus que h*                       p7502420
c                                                                       p7502430
c   impression sur cartes                                               p7502440
c                                                                       p7502450
      ipr=7505                                                          p7502460
!      write(7,7004) ia,ib,ic,id,ha,prma,pprma,tseta,ipr                 p7502470
 7004 format(4i5,4f15.8,8x,'in',i5)                                     p7502480
      do 50 i=1,ia                                                      p7502490
!      write(7,7000) i,aa(i),a(i),c(i),ipr                               p7502500
 7000 format(i4,f20.8,f20.7,f20.6,10x,'a',i5)                           p7502510
   50 continue                                                          p7502520
      do 51 i=1,ib                                                      p7502530
!      write(7,7001) i,bb(i),b(i),d(i),ipr                               p7502540
 7001 format(i4,3f20.10,10x,'b',i5)                                     p7502550
   51 continue                                                          p7502560
      do 52 i=1,ic                                                      p7502570
      baf=bf(i)/pirr                                                    p7502580
      paf=pf(i)/pirr                                                    p7502590
      dpaf=dpf(i)/pirr                                                  p7502600
      s=sa(i)/pirr                                                      p7502610
      d11=ddr(i)/pir                                                    p7502620
c     if(dabs(baf).ge.0.1d0)goto 54
c     if(dabs(paf).ge.1.0d0)goto 54
c     goto 52
!54    write(7,7002) i,baf,paf,dpaf,s,d11,ipr                            p7502630
 7002 format(i4,f15.7,f15.7,f13.7,f13.6,f13.5,' bp',i5)                  p7502640
   52 continue                                                          p7502650
      icc=ic+1                                                          p7502660
      do 53 i=icc,id                                                    p7502670
      paf=pf(i)/pirr                                                    p7502680
      dpaf=dpf(i)/pirr                                                  p7502690
      s=sa(i)/pirr                                                      p7502700
      d11=ddr(i)/pir                                                    p7502710
c     if(dabs(baf).ge.0.1d0)goto 56
c     if(dabs(paf).ge.1.0d0)goto 56
c     goto 53
!   56 write(7,7003) i,paf,dpaf,s,d11,ipr                                p7502720
 7003 format(i4,15x,f15.7,f13.7,f13.6,f13.5,' bp',i5)                    p7502730
   53 continue                                                          p7502740
c   
c	stop                                                            p7502750
c       goto 10
c   test obl prec tprec                                                 p7502760
c                                                                       p7502770
*      write(6,6502)                                                    p7502780
 6502 format(1h1,///,1x,'test valeurs obl prec tprec for dates from 1950p7502790
     *',//)                                                             p7502800
      zer=0.0d0                                                         p7502810
      bvo=zer                                                           p7502820
      pvo=zer                                                           p7502830
      dpvo=zer                                                          p7502840
      do 60 i=1,ic                                                      p7502850
      bvo=bvo+bf(i)*dcos(ddr(i))                                        p7502860
      pvo=pvo+pf(i)*dsin(ddr(i))                                        p7502870
      dpvo=dpvo+dpf(i)*dcos(ddr(i))                                     p7502880
   60 continue                                                          p7502890
      ic1=ic+1                                                          p7502900
      do 61 i=ic1,id                                                    p7502910
      pvo=pvo+pf(i)*dsin(ddr(i))                                        p7502920
      dpvo=dpvo+dpf(i)*dcos(ddr(i))                                     p7502930
   61 continue                                                          p7502940
      bvo=bvo/pir     +ha                                               p7502950
      pvo=pvo/pir     +tseta                                            p7502960
      dpvo=dpvo/pirr+pprma                                              p7502970
      iwr=1950                                                          p7502980
*      write(6,6503) iwr,bvo,pvo,dpvo                                   p7502990
*      write(6,6910)                                                    p7503000
      apa=1000.0d0                                                      p7503010
      tc=0.0d0                                                          p7503020
      do 14 kk=1,6                                                      p7503030
      tcc=(kk-1)*100000.0d0                                             p7503040
      do 11 i=1,10                                                      p7503050
      t=tc-(i-1)*apa-tcc                                                p7503060
      bvo=zer                                                           p7503070
      pvo=pprm *t                                                       p7503080
      dpvo=pprm                                                         p7503090
      do 12 j=1,ic                                                      p7503100
      wa=sa(j)*t+ddr(j)                                                 p7503110
      bvo=bvo+bf(j)*dcos(wa)                                            p7503120
      pvo=pvo+pf(j)*dsin(wa)                                            p7503130
      dpvo=dpvo+dpf(j)*dcos(wa)                                         p7503140
   12 continue                                                          p7503150
      ic1=ic+1                                                          p7503160
      do 13 j=ic1,id                                                    p7503170
      wa=sa(j)*t+ddr(j)                                                 p7503180
      pvo=pvo+pf(j)*dsin(wa)                                            p7503190
      dpvo=dpvo+dpf(j)*dcos(wa)                                         p7503200
   13 continue                                                          p7503210
      iwr=(t-tc)/apa                                                    p7503220
      bvo=bvo/pir                                                       p7503230
      bvo=bvo+ha                                                        p7503240
      pvo=pvo/pir+tseta                                                 p7503250
      dpvo=dpvo/pirr                                                    p7503260
   19 if(dabs(pvo).le.360.0d0) go to 18                                 p7503270
      pvo=dsign(1.0d0,pvo)*(dabs(pvo)-360.0d0)                          p7503280
      go to 19                                                          p7503290
   18 continue                                                          p7503300
*      write(6,6503) iwr,bvo,pvo,dpvo                                   p7503310
 6503 format(i10,f10.4,f12.4,f12.6)                                     p7503320
   11 continue                                                          p7503330

   14 continue                                                          p7503360
c                                                                       p7503370
c   calcul  e*sin(long per mobil)                                       p7503380
c                                                                       p7503390
!      write(6,6951)                                                     p7503400
! 6951 format(1h1,//,1x,'parameters of series expansion',//)             p7503410
!      write(6,6952)                                                     p7503420
! 6952 format(1x,'eccentricity and longitude of the moving perihelion',/)p7503430
!      write(6,6953)                                                     p7503440
! 6953 format(3x,'kk',2x,'j',2x,'i',5x,' amplitude',4x,'   mean rate ',3xp7503450
!     *,5x,'phase  ',3x,6x,'period   ',/)                                p7503460
!      write(6,6954)                                                     p7503470
 6954 format(32x,'(''''/years)',8x,'(degree)',9x,'(years)',/)           p7503480
      do 90 j=1,ia                                                      p7503490
      qaa(j)=aa(j)                                                      p7503500
      qa(j)=a(j)+prma                                                   p7503510
      qc(j)=c(j)+tseta                                                  p7503520
      wq=qc(j)                                                          p7503530
c   wqa est toujours positif                                            p7503540
      call phase(wq)                                                    p7503550
      qc(j)=wq                                                          p7503560
      wqaa=qaa(j)                                                       p7503570
      wqa=qa(j)                                                         p7503580
      wqc=qc(j)                                                         p7503590
      i=0                                                               p7503600
c   premier ind rappelle ind e*sin                                      p7503610
c   deuxieme  ind rapprlle  precession                                  p7503620
c   troisieme ind est le numero d ordre final                           p7503630
      call awri(wqaa,wqa,wqc,j,i,j)                                     p7503640
   90 continue                                                          p7503650
!      write(6,6950)                                                     p7503660
 6950 format(1x)                                                        p7503670
      do 91 j=1,ia                                                      p7503680
      do 92 i=1,ib                                                      p7503690
      j1=ia+(j-1)*ib+i                                                  p7503700
      qaa(j1)=aa(j)*pf(i)*0.5d0                                         p7503710
      qa(j1)=a(j)+b(i)+prma+prma                                        p7503720
      qc(j1)=c(j)+d(i)+2.0d0*tseta                                      p7503730
      wq=qc(j1)                                                         p7503740
      call phase(wq)                                                    p7503750
      qc(j1)=wq                                                         p7503760
      wqaa=qaa(j1)                                                      p7503770
      wqa=qa(j1)                                                        p7503780
      wqc=qc(j1)                                                        p7503790
      call awri(wqaa,wqa,wqc,j,i,j1)                                    p7503800
      j2=ia*ib+j1                                                       p7503810
      qaa(j2)=qaa(j1)                                                   p7503820
      qa(j2)=a(j)-b(i)                                                  p7503830
      qc(j2)=180.0d0-d(i)+c(j)                                          p7503840
      wq=qc(j2)                                                         p7503850
      call phase(wq)                                                    p7503860
      qc(j2)=wq                                                         p7503870
   92 continue                                                          p7503880
   91 continue                                                          p7503890
!      write(6,6950)                                                     p7503900
c                                                                       p7503910
c   impression du troisieme groupe pour ne pas melanger avec le deuxiemep7503920
c                                                                       p7503930
      j3=ia*(ib+1)                                                      p7503940
      j3=j3+1                                                           p7503950
      j4=ia*(2*ib+1)                                                    p7503960
      i=0                                                               p7503970
      j=1                                                               p7503980
      do 93 j1=j3,j4                                                    p7503990
      i=i+1                                                             p7504000
      if(i.le.ib) go to 94                                              p7504010
      i=1                                                               p7504020
      j=j+1                                                             p7504030
   94 waa=qaa(j1)                                                       p7504040
      wa=qa(j1)                                                         p7504050
      wc=qc(j1)                                                         p7504060
      call awri(waa,wa,wc,j,i,j1)                                       p7504070
   93 continue                                                          p7504080
c                                                                       p7504090
c   elimination des valeurs inf a 10(-05)                               p7504100
c                                                                       p7504110
      acons=1296.0d03                                                   p7504120
      tes=1.0d-05                                                       p7504130
      ji=0                                                              p7504140
      zz=0.0d0                                                          p7504150
!      write(6,6957)                                                     p7504160
 6957 format(1h1,///)                                                   p7504170
!      write(6,6952)                                                     p7504180
!      write(6,6955)                                                     p7504190
 6955 format(3x,'i',6x,' amplitude',4x,'   mean rate ',3x,6x,'phase',3x,p7504200
     *7x,'period',/)                                                    p7504210
!      write(6,6956)                                                     p7504220
 6956 format(27x,'(''''/year)',8x,'(degree)',9x,'(years)',/)            p7504230
      jk=0                                                              p7504240
      do 95 j=1,j4                                                      p7504250
      waa=dabs(qaa(j))                                                  p7504260
      if (waa.le.tes) go to 95                                          p7504270
      ji=ji+1                                                           p7504280
      if(qa(j).eq.zz) go to 88                                          p7504290
      go to 89                                                          p7504300
   88 aper=zz                                                           p7504310
   89 aper=acons/dabs(qa(j))                                            p7504320
      jk=jk+1                                                           p7504330
      ind(jk)=j                                                         p7504340
!      write(6,6958) ji,qaa(j),qa(j),qc(j),aper                          p7504350
!      write(7,6958) ji,qaa(j),qa(j),qc(j),aper                          p7504360
 6958 format(1x,i4,2x,d15.8,1x,f13.6,3x,f13.6,3x,f13.0)                 p7504370
   95 continue                                                          p7504380
c                                                                       p7504390
c   test e*sin(xl)                                                      p7504400
c                                                                       p7504410
!      write(6,6504)                                                     p7504420
 6504 format(1h1,///,1x,'test valeurs exc  long per mob for dates from  p7504430
     *1950',//)                                                         p7504440
!      write(7,9502) jk                                                  p7504450
 9502 format(i10)                                                       p7504460
!      write(6,9500) jk                                                  p7504470
! 9500 format(1x,'for limited number of terms    jk=',i8,//)             p7504480
      do 107 kk=1,6                                                     p7504490
      tcc=(kk-1)*100000.0d0                                             p7504500
      do 106 i=1,10                                                     p7504510
      t=tc-(i-1)*apa-tcc                                                p7504520
      es=zer                                                            p7504530
      ec=zer                                                            p7504540
      do 100 jf=1,jk                                                    p7504550
      j=ind(jf)                                                         p7504560
      ta=qa(j)*t*pirr+qc(j)*pir                                         p7504570
      es=es+qaa(j)*dsin(ta)                                             p7504580
      ec=ec+qaa(j)*dcos(ta)                                             p7504590
  100 continue                                                          p7504600
      call excen(es,ec,pi,zer,ee,peva,pir)                              p7504610
      iwr=(t-tc)/apa                                                    p7504620
!      write(6,6505) iwr,ee,peva,es                                      p7504630
  106 continue                                                          p7504640
!      write(6,6910)                                                     p7504650
  107 continue                                                          p7504660
!      write(6,6910)                                                     p7504670
!      write(6,6910)                                                     p7504680
      ie=ia+ia*ib*2                                                     p7504690
!      write(6,9501) ie                                                  p7504700
 9501 format(1x,'for all terms   ie=',i8,//)                            p7504710
      do 17 kk=1,6                                                      p7504720
      tcc=(kk-1)*100000.0d0                                             p7504730
      do 16 i=1,10                                                      p7504740
      t=tc-(i-1)*apa-tcc                                                p7504750
      es=zer                                                            p7504760
      ec=zer                                                            p7504770
      do 15 j=1,ie                                                      p7504780
      ta=qa(j)*t*pirr+qc(j)*pir                                         p7504790
      ec=ec+dcos(ta)*qaa(j)                                             p7504800
      es=es+dsin(ta)*qaa(j)                                             p7504810
   15 continue                                                          p7504820
      call excen(es,ec,pi,zer,ee,peva,pir)                              p7504830
      iwr=(t-tc)/apa                                                    p7504840
!      write(6,6505) iwr,ee,peva,es                                      p7504850
! 6505 format(i10,f12.8,f12.5,f10.6)                                     p7504860
   16 continue                                                          p7504870
!      write(6,6910)                                                     p7504880
   17 continue                                                          p7504890
10    continue
      close(10)
      end subroutine                                                    p7504910 
      subroutine impd(ia,ib,aa,a,c,bb,b,d,kz,ar,cr,br,dr,pip,pirr,pir,anp7504920
     *umi)                                                              p7504930
      implicit double precision (a-h,p-z)                               p7504940
      dimension aa(ia),a(ia),c(ia),ar(ia),cr(ia),bb(ib),b(ib),d(ib),br(ip7504950
     *b),dr(ib),cyc1(80),cyc2(500)                                      p7504960
       common kktilde
c                                                                       p7504970
c impression data + transformation en radians                           p7504980
c anumi=1.0 signifie que on a sin(i)                                    p7504990
c anumi=2.0 signifie que on a sin(i/2)                                  p7505000
c                                                                       p7505010
      z=0.0d0                                                           p7505020
      zz=1.0d-08                                                        p7505030
      emax=z                                                            p7505040
      do 3 k=1,ia                                                       p7505050
      emax=emax+dabs(aa(k))                                             p7505060
      ar(k)=a(k)*pirr                                                   p7505070
      cr(k)=c(k)*pir                                                    p7505080
      aza=dabs(ar(k))                                                   p7505090
      if(aza.le.zz) go to 41                                            p7505100
      cyc1(k)=pip/aza                                                   p7505110
      go to 3                                                           p7505120
c si periode infinie j'ecris 0.0                                        p7505130
   41 cyc1(k)=z                                                         p7505140
    3 continue                                                          p7505150
      aim=z                                                             p7505160
      do 40 k=1,ib                                                      p7505170
      aim=aim+dabs(bb(k))                                               p7505180
      br(k)=b(k)*pirr                                                   p7505190
      dr(k)=d(k)*pir                                                    p7505200
      aza=dabs(br(k))                                                   p7505210
      if (aza.le.zz) go to 42                                           p7505220
      cyc2(k)=pip/aza                                                   p7505230
      go to 40                                                          p7505240
   42 cyc2(k)=z                                                         p7505250
   40 continue                                                          p7505260
      ain=aim-dabs(bb(kz))                                              p7505270
      aim=dasin(aim)                                                    p7505280
      ain=dasin(ain)                                                    p7505290
      aim=aim/pir                                                       p7505300
      ain=ain/pir                                                       p7505310
!      write(6,6009)                                                     p7505320
 6009 format(13x,'aa',13x,'a',13x,'c',12x,'bb',13x,'b',13x,'d',9x,'periop7505330
     *de a',5x,'periode b',/)                                           p7505340
      if (ia.ge.ib) go to 20                                            p7505350
      ix=ia                                                             p7505360
      go to 21                                                          p7505370
   20 ix=ib                                                             p7505380
!   21 write(6,6002)(i,aa(i),a(i),c(i),bb(i),b(i),d(i),cyc1(i),cyc2(i),i=p7505390
!     *1,ix)                                                             p7505400
   21 continue
 6002 format(i6,f14.8,f14.7,f14.6,f14.8,f14.7,f14.6,2f14.1)             p7505410
      if (ia.ge.ib) go to 22                                            p7505420
      l=ia+1                                                            p7505430
!      write(6,6021) (i,bb(i),b(i),d(i),cyc2(i),i=l,ib)                  p7505440
 6021 format(i6,42x,f14.8,f14.7,f14.6,14x,f14.1)                        p7505450
      go to 11                                                          p7505460
   22 if (ia.eq.ib) go to 11                                            p7505470
      l=ib+1                                                            p7505480
!      write(6,6020)(i,aa(i),a(i),c(i),cyc1(i),i=l,ia)                   p7505490
 6020 format(i6,f14.8,f14.7,f14.6,42x,f14.1)                            p7505500
   11 continue
 6007 format(/)                                                         p7505520
      aim=aim*anumi                                                     p7505530
      ain=ain*anumi                                                     p7505540
!      write(6,6010) emax,aim,ain                                        p7505550
! 6010 format(1h ,'max excen =',f9.7,5x,'max incl eclip origine =',f10.6,p7505560
!     *15x,'max incl plan inv =',f10.6)                                  p7505570
      return                                                            p7505580
      end                                                               p7505590
      subroutine trani(ibt,bbt,bt,dt,ibs,bb,b,d,ib)                     p7505600
      implicit double precision (a-h,p-z)                               p7505610
      dimension bbt(ibt),bt(ibt),dt(ibt),bb(ibs),b(ibs),d(ibs)          p7505620
       common kktilde
c                                                                       p7505630
c calcul du developpement de sin(i) from sin(i/2)                       p7505640
c nombre max. de termes presuppose : ibs                                p7505650
c   nombre reel de termes sup  zz   = ib                                p7505660
c bt,b sont en secondes                                                 p7505670
c b,d sont en degres                                                    p7505680
c                                                                       p7505690
!      write(6,6000)                                                     p7505700
! 6000 format(1h1,//,10x,'coeff developpement sin berger-laskar',//)     p7505710
!      write(6,6001)                                                     p7505720
! 6001 format(1h ,'   kk   j   i   k',10x,'bb',12x,'b',13x,'d',/)        p7505730
      zz=1.0d-05                                                        p7505740
      do 1 i=1,ibt                                                      p7505750
      x=0.0d0                                                           p7505760
      y=bbt(i)*bbt(i)                                                   p7505770
      do 2 j=1,ibt                                                      p7505780
      if(j.eq.i) go to 2                                                p7505790
      x=x+bbt(j)*bbt(j)                                                 p7505800
    2 continue                                                          p7505810
      kk=i                                                              p7505820
      bb(i)=bbt(i)*(2.0d0-y-2.0d0*x)                                    p7505830
      b(i)=bt(i)                                                        p7505840
      d(i)=dt(i)                                                        p7505850
!      write(6,6002) i,i,bb(i),b(i),d(i)                                 p7505860
! 6002 format(1h ,i5,i4,8x,f14.8,f14.7,f14.6)                            p7505870
    1 continue                                                          p7505880
!*      write(6,6003)                                                    p7505890
 6003 format(//)                                                        p7505900
      do 3 j=1,ibt                                                      p7505910
      do 4 i=1,ibt                                                      p7505920
      if(i.eq.j) go to 4                                                p7505930
      x=-bbt(i)*bbt(i)                                                  p7505940
      y=2.0d0*bt(i)                                                     p7505950
      z=2.0d0*dt(i)                                                     p7505960
      x=bbt(j)*x                                                        p7505970
      xtra=y-bt(j)                                                      p7505980
      ijcon=0                                                           p7505990
      do 60 ixy=1,kk                                                    p7506000
c   *** il ne peut y avoir qu une seule valeur = xtra car je teste chaqup7506010
      if (xtra.ne.b(ixy)) go to 60                                      p7506020
      bb(ixy)=bb(ixy)+x                                                 p7506030
      ijcon=1                                                           p7506040
      ijl=ixy                                                           p7506050
   60 continue                                                          p7506060
      if (ijcon.eq.0) go to 61                                          p7506070
      ixyz=0                                                            p7506080
*      write(6,6004) ixyz,j,i,x,xtra,d(ijl)                             p7506090
      go to 4                                                           p7506100
   61 continue                                                          p7506110
c   je ne garde pas de termes dont b n est pas dans les 15 premiers et qp7506120
c          inf a zz *** si il est dans les 15 premiers il est ajoute immp7506130
      xabs=dabs(x)                                                      p7506140
      if (xabs.lt.zz) go to 4                                           p7506150
      kk=kk+1                                                           p7506160
      bb(kk)=x                                                          p7506170
      b(kk)=y-bt(j)                                                     p7506180
      d(kk)=z-dt(j)                                                     p7506190
    6 if (d(kk).lt.360.0d0) go to 5                                     p7506200
      d(kk)=d(kk)-360.0d0                                               p7506210
      go to 6                                                           p7506220
    5 if (d(kk)) 8,7,7                                                  p7506230
    8 d(kk)=d(kk)+360.0d0                                               p7506240
      go to 5                                                           p7506250
    7 continue
!       write(6,6004) kk,j,i,bb(kk),b(kk),d(kk)                          p7506260
! 6004 format(1h ,i5,2i4,4x,f14.8,f14.7,f14.6)                           p7506270
    4 continue                                                          p7506280
 6005 format(/)                                                         p7506290
    3 continue                                                          p7506300
!*      write(6,6005)                                                    p7506310
      do 18 j=1,ibt                                                     p7506320
      iii=ibt-1                                                         p7506330
      do 9 i=1,iii                                                      p7506340
      if (i.eq.j) go to 9                                               p7506350
      x=bbt(j)*bbt(i)                                                   p7506360
      y=-bt(j)+bt(i)                                                    p7506370
      z=-dt(j)+dt(i)                                                    p7506380
      ii=i+1                                                            p7506390
      do 10 k=ii,ibt                                                    p7506400
      if (k.eq.j) go to 10                                              p7506410
      x=-2.0d0*bbt(k)*x                                                 p7506420
      xtra=y+bt(k)                                                      p7506430
      ijcon=0                                                           p7506440
      do 62 ixy=1,kk                                                    p7506450
c   *** il ne peut y avoir qu une seule valeur = xtra car je teste chaqup7506460
      if (xtra.ne.b(ixy)) go to 62                                      p7506470
      bb(ixy)=bb(ixy)+x                                                 p7506480
      ijcon=1                                                           p7506490
      ijl=ixy                                                           p7506500
   62 continue                                                          p7506510
      if (ijcon.eq.0) go to 63                                          p7506520
      ixyz=0                                                            p7506530
*      write(6,6006) ixyz,j,i,k,x,xtra,d(ijl)                           p7506540
      go to 10                                                          p7506550
   63 continue                                                          p7506560
c   je ne garde pas de termes dont b n est pas dans les 15 premiers et qp7506570
c          inf a zz *** si il est dans les 15 premiers il est ajoute immp7506580
      xabs=dabs(x)                                                      p7506590
      if (xabs.lt.zz) go to 10                                          p7506600
      kk=kk+1                                                           p7506610
      bb(kk)=x                                                          p7506620
      b(kk)=y+bt(k)                                                     p7506630
      d(kk)=z+dt(k)                                                     p7506640
   12 if (d(kk).lt.360.0d0) go to 14                                    p7506650
      d(kk)=d(kk)-360.0d0                                               p7506660
      go to 12                                                          p7506670
   14 if (d(kk)) 13,11,11                                               p7506680
   13 d(kk)=d(kk)+360.0d0                                               p7506690
      go to 14                                                          p7506700
   11 continue
!      write(6,6006) kk,j,i,k,bb(kk),b(kk),d(kk)                         p7506710
 6006 format(1h ,i5,3i4,f14.8,f14.7,f14.6)                              p7506720
   10 continue                                                          p7506730
    9 continue                                                          p7506740
   18 continue                                                          p7506750
      ib=kk                                                             p7506760
!      write(6,6005)                                                     p7506770
!      write(6,6005)                                                     p7506780
!      write(6,6011)                                                     p7506790
 6011 format(1h ,'test valeurs incl noeud pour qques valeurs t',/)      p7506800
c                                                                       p7506810
c test sur valeurs de incl. pour quelques dates                         p7506820
c                                                                       p7506830
      z=0.0d0                                                           p7506840
      kcon=0                                                            p7506850
      pi=dacos(-1.0d0)                                                  p7506860
      pir=pi/180.0d0                                                    p7506870
      pirr=pir/3600.0d0                                                 p7506880
      t=0.0d0                                                           p7506890
   25 xt1=0.0d0                                                         p7506900
      yt1=0.0d0                                                         p7506910
      do 20 i=1,ibt                                                     p7506920
      xp=bt(i)*pirr                                                     p7506930
      brt=xp*t+dt(i)*pir                                                p7506940
      xt1=xt1+bbt(i)*dsin(brt)                                          p7506950
      yt1=yt1+bbt(i)*dcos(brt)                                          p7506960
   20 continue                                                          p7506970
      xt2=0.0d0                                                         p7506980
      yt2=0.0d0                                                         p7506990
      do 21 i=1,ib                                                      p7507000
      xp=b(i)*pirr                                                      p7507010
      brt=xp*t+d(i)*pir                                                 p7507020
      xt2=xt2+bb(i)*dsin(brt)                                           p7507030
      yt2=yt2+bb(i)*dcos(brt)                                           p7507040
   21 continue                                                          p7507050
      xx=xt1                                                            p7507060
      yy=yt1                                                            p7507070
      icont=0                                                           p7507080
   22 continue                                                          p7507090
      if(yy) 39,35,40                                                   p7507100
   39 ano=datan(xx/yy)+pi                                               p7507110
      go to 34                                                          p7507120
   40 ano=datan(xx/yy)                                                  p7507130
      if (xx.ge.z) go to 34                                             p7507140
      ano=2.0d0*pi+ano                                                  p7507150
      go to 34                                                          p7507160
   35 if (xx) 36,37,370                                                 p7507170
   36 ano=pi+pi/2.0d0                                                   p7507180
      go to 34                                                          p7507190
  370 ano=pi/2.0d0                                                      p7507200
      go to 34                                                          p7507210
   37 ano=z                                                             p7507220
   34 ain=dasin(dsqrt(xx*xx+yy*yy))                                     p7507230
      if (icont.eq.0) go to 23                                          p7507240
      go to 24                                                          p7507250
   23 anoa=ano/pir                                                      p7507260
      aina=ain/pir                                                      p7507270
      aina=aina*2.0d0                                                   p7507280
      xx=xt2                                                            p7507290
      yy=yt2                                                            p7507300
      icont=1                                                           p7507310
      go to 22                                                          p7507320
   24 anoa1=ano/pir                                                     p7507330
      aina1=ain/pir                                                     p7507340
!      write(6,6010) t,anoa,anoa1,aina,aina1                             p7507350
 6010 format(1h ,f10.0,4f16.8,/)                                        p7507360
      kcon=kcon+1                                                       p7507370
      go to (26,27,28,29,31,32,33,54,55,30),kcon                        p7507380
   26 t=1000.0d0                                                        p7507390
      go to 25                                                          p7507400
   27 t=10000.0d0                                                       p7507410
      go to 25                                                          p7507420
   28 t=100000.0d0                                                      p7507430
      go to 25                                                          p7507440
   29 t=500000.0d0                                                      p7507450
      go to 25                                                          p7507460
   31 t=-100.0d0                                                        p7507470
      go to 25                                                          p7507480
   32 t=-1000.0d0                                                       p7507490
      go to 25                                                          p7507500
   33 t=-10000.0d0                                                      p7507510
      go to 25                                                          p7507520
   54 t=-100000.0d0                                                     p7507530
      go to 25                                                          p7507540
   55 t=-500000.0d0                                                     p7507550
      go to 25                                                          p7507560
   30 continue                                                          p7507570
      return                                                            p7507580
      end                                                               p7507590
      subroutine coef(h,prm,tset,apo,al,icim,ia,aa,a,c,ar,cr,ib,bb,b,d, p7507600
     *br,dr,bf,pf,dpf,sa,ddr,ic,id,pi,pir,pirr,hh)                      p7507610
      implicit double precision (a-h,p-z)                               p7507620
      dimension aa(ia),a(ia),c(ia),ar(ia),cr(ia),bb(ib),b(ib),d(ib),br(ip7507630
     *b),dr(ib),bf(ic),pf(id),dpf(id),sa(id),ddr(id),jm(5)              p7507640
       common kktilde
c                                                                       p7507650
c calcul des coefficients des developpements de obliquite mobile        p7507660
c   precession mobile   taux prec mobile                                p7507670
c valeurs sont en radians. - pas impression : icim=0                    p7507680
c impression en secondes, pulsations en secondes, phases en degres icim=p7507690
c bf vecteurs de ic termes
c pf dpf vecteurs de id termes                                          p7507710
c                                                                       p7507720
      if (icim-1) 320,308,320                                           p7507730
  320 continue                                                          p7507740
      ud=0.5d0                                                          p7507750
      un=1.0d0                                                          p7507760
      th=dtan(h)                                                        p7507770
      twt=1.0d0/dtan(h)                                                 p7507780
      do 302 i=1,ib                                                     p7507790
      sa(i)=br(i)+prm                                                   p7507800
      ddr(i)=dr(i)+tset                                                 p7507810
      ii=ib+i                                                           p7507820
      sa(ii)=2.0d0*sa(i)                                                p7507830
      ddr(ii)=2.0d0*ddr(i)                                              p7507840
  302 continue                                                          p7507850
c                                                                       p7507860
c bf,pf ... ne sont dans un premier calcul pas multiples par n ...      p7507870
c                                                                       p7507880
      do 300 i=1,ib                                                     p7507890
      w=prm/sa(i)                                                       p7507900
      ww=w*w                                                            p7507910
      ii=ib+i                                                           p7507920
      bf(i)=w                                                           p7507930
      bf(ii)=0.25d0*(w*(ww-un)*th+ww*twt)                               p7507940
      pf(i)=w*(w-un)*th+w*twt                                           p7507950
      pf(ii)=0.25d0*w*(ww+w-un)+0.125d0*ww*(w-un)*(w-un)*th*th+0.5d0*ww*p7507960
     *twt*twt                                                           p7507970
  300 continue                                                          p7507980
      kt=0                                                              p7507990
      iii=ib-1                                                          p7508000
      do 303 i=1,iii                                                    p7508010
      ii=i+1                                                            p7508020
      w=prm/sa(i)                                                       p7508030
      ww=w*w                                                            p7508040
      aba=bf(i)*bf(i)                                                   p7508050
      do 301 k=ii,ib                                                    p7508060
      kt=kt+1                                                           p7508070
      kk=2*ib+kt                                                        p7508080
      kkk=2*ib+(ib*(ib-1))/2+kt                                         p7508090
      wk=prm/sa(k)                                                      p7508100
      wwk=wk*wk                                                         p7508110
      sa(kk)=sa(i)+sa(k)                                                p7508120
      sa(kkk)=br(i)-br(k)                                               p7508130
      ddr(kk)=ddr(i)+ddr(k)                                             p7508140
      ddr(kkk)=dr(i)-dr(k)                                              p7508150
      bf(kk)=0.5d0*prm/sa(kk)*(th*(ww+wwk-2.0d0)+twt*(w+wk))            p7508160
      bf(kkk)=0.5d0*prm/sa(kkk)*(th*(ww-wwk+2.0d0*(wk- w))+twt*(w-wk))  p7508170
      pf(kk)=prm/sa(kk)*(0.5d0*(ww+wwk+w+wk-w*wk-2.0d0)+bf(kk)*th-0.5d0*p7508180
     *th*th*(w*(w-un)+wk*(wk-un))+twt*twt*(w+wk))                       p7508190
      pf(kkk)=prm/sa(kkk)*(0.5d0*(5.0d0*(w+wk)-(ww+wwk)-w*wk-6.0d0)+bf(kp7508200
     *kk)*th+0.5d0*th*th*(w*(w-un)+wk*(wk-un)))                         p7508210
      pfik=0.5d0*(pf(i)+pf(k))                                          p7508220
      twt2=0.5d0*twt                                                    p7508230
      abb=bf(k)*bf(k)                                                   p7508240
      bf(kk)=bf(kk)-pfik+twt2                                           p7508250
      bf(kkk)=bf(kkk)+un*pfik-ud*twt                                    p7508260
      pf(kk)=pf(kk)-(bf(i)+bf(k)-un)*twt*twt-ud*(aba+abb-un)            p7508270
      pf(kkk)=pf(kkk)-ud*(aba-abb-2.0d0*bf(i)+2.0d0*bf(k))              p7508280
      dpf(kk)=sa(kk)*pf(kk)                                             p7508290
      dpf(kkk)=sa(kkk)*pf(kkk)                                          p7508300
c                                                                       p7508310
c multiplication par nn                                                 p7508320
c                                                                       p7508330
      t=bb(i)*bb(k)                                                     p7508340
      bf(kk)=bf(kk)*t*(-1.0d0)                                          p7508350
      bf(kkk)=bf(kkk)*t*(-1.0d0)                                        p7508360
      pf(kk)=pf(kk)*t                                                   p7508370
      pf(kkk)=pf(kkk)*t                                                 p7508380
      dpf(kk)=dpf(kk)*t                                                 p7508390
      dpf(kkk)=dpf(kkk)*t                                               p7508400
  301 continue                                                          p7508410
  303 continue                                                          p7508420
      hht=0.0d0                                                         p7508430
      do 304 i=1,ib                                                     p7508440
      tw=bb(i)*bb(i)                                                    p7508450
      xba=bf(i)                                                         p7508460
      xbb=pf(i)                                                         p7508470
      hht=hht+tw*(ud*bf(i)*(bf(i)-un)*th+0.25d0*twt*(2.0d0*bf(i)-un))   p7508480
      bf(i)=xba-un                                                      p7508490
      pf(i)=xbb-twt                                                     p7508500
      bf(i)=bf(i)*bb(i)*(-1.0d0)                                        p7508510
      pf(i)=pf(i)*bb(i)                                                 p7508520
      dpf(i)=pf(i)*sa(i)                                                p7508530
      ii=ib+i                                                           p7508540
      tw=bb(i)*bb(i)                                                    p7508550
      bf(ii)=bf(ii)-0.5d0*xbb+0.25d0*twt                                p7508560
      pf(ii)=pf(ii)-(xba-ud)*twt*twt-ud*(xba*xba-ud)                    p7508570
      bf(ii)=bf(ii)*tw*(-1.0d0)                                         p7508580
      pf(ii)=pf(ii)*tw                                                  p7508590
      dpf(ii)=pf(ii)*sa(ii)                                             p7508600
  304 continue                                                          p7508610
      hh=h-hht                                                          p7508620
c                                                                       p7508630
c calcul des contributions en m. ici ampl immediatement multipliee par mp7508640
c                                                                       p7508650
      kt=0                                                              p7508660
      ii=ia-1                                                           p7508670
      tw=3.0d0*apo/al*prm                                               p7508680
      do 305 i=1,ii                                                     p7508690
      ki=i+1                                                            p7508700
      do 306 k=ki,ia                                                    p7508710
      kt=kt+1                                                           p7508720
      ik=ic+kt                                                          p7508730
      sa(ik)=ar(i)-ar(k)                                                p7508740
      pf(ik)=tw/(ar(i)-ar(k))*aa(i)*aa(k)                               p7508750
c
c   prn .ne. pprn
c
      if (kktilde.eq.0)then
       dpf(ik)=-2*tw*aa(i)*aa(k)                                        p7508760
      endif
c
c    prn .eq.pprn
c
c     if (kktilde.eq.1)then
       dpf(ik)=pf(ik)*sa(ik)                                            p7508760
c     endif
c
      ddr(ik)=cr(i)-cr(k)                                               p7508770
  306 continue                                                          p7508780
  305 continue                                                          p7508790
c                                                                       p7508800
c impression                                                            p7508810
c                                                                       p7508820
      if(icim-1) 307,308,307                                            p7508830
!  308 write(6,6000)                                                     p7508840
  308  continue
! 6000 format(1h1,///,10x,'parametres des developpements berger-sharaf 2'p7508850
!     *,/)                                                               p7508860
!      write (6,6001)                                                    p7508870
! 6001 format(1h ,9x,'obliquite mobi',5x,'precession mobi',4x,'taux de prp7508880
!     1ecession mobi',/)                                                 p7508890
!      write(6,6002)                                                     p7508900
! 6002 format(1h ,12x,'amplitude',11x,'amplitude',13x,'amplitude',6x,'pulp7508910
!     1sation',6x,'phase',/)                                             p7508920
      jm(1)=ib                                                          p7508930
      jm(2)=ib                                                          p7508940
      jm(3)=(ib*(ib-1))/2                                               p7508950
      jm(4)=jm(3)                                                       p7508960
      jm(5)=(ia*(ia-1))/2                                               p7508970
      i=0                                                               p7508980
      do 309 k=1,5                                                      p7508990
      jk=jm(k)                                                          p7509000
      do 316 kk=1,jk                                                    p7509010
      i=i+1                                                             p7509020
      s=sa(i)/pirr                                                      p7509030
      d11=ddr(i)/pir                                                    p7509040
c                                                                       p7509050
c phase comprise entre 0 et 360                                         p7509060
c                                                                       p7509070
  313 if (d11-360.0d0) 310,314,312                                      p7509080
  312 d11=d11-360.0d0                                                   p7509090
      go to 313                                                         p7509100
  310 if (d11+360.0d0) 315,314,314                                      p7509110
  315 d11=d11+360.0d0                                                   p7509120
      go to 310                                                         p7509130
  314 paf=pf(i)/pirr                                                    p7509140
      dpaf=dpf(i)/pirr                                                  p7509150
      if (k.eq.5) go to 317                                             p7509160
      baf=bf(i)/pirr                                                    p7509170
      if((dabs(baf).lt.1.0).and.(dabs(paf).lt.1.1))goto 316
!      write(6,6003) i,baf,paf,dpaf,s,d11                                p7509180
! 6003 format(1h ,i4,3x,f16.7,4x,f16.7,4x,f16.7,3x,f12.6,2x,f12.5)       p7509190
      go to 316                                                         p7509200
  317 if(dabs(paf).lt.1.1)goto 316
!      write(6,6004) i,paf,dpaf,s,d11                                    p7509210
 6004 format(1h ,i4,3x,20x,f16.7,4x,f16.7,3x,f12.6,2x,f12.5)            p7509220
  316 continue                                                          p7509230
!      write (6,6005)                                                    p7509240
! 6005 format(/)                                                         p7509250
  309 continue                                                          p7509260
  307 continue                                                          p7509270
      return                                                            p7509280
      end                                                               p7509290
      subroutine fbpf(t,ic,id,bf,pf,sa,ddr,fbf,fpf)                     p7509300
      implicit double precision (a-h,p-z)                               p7509310
      dimension bf(ic),pf(id),sa(id),ddr(id)                            p7509320
       common kktilde
c                                                                       p7509330
c calcul somme des termes obliquite fixe                                p7509340
c calcul somme des termes precession fixe                               p7509350
c ic= nombre termes obl. fixe                                           p7509360
c id= nombre termes prec. fixe                                          p7509370
c                                                                       p7509380
      z=0.0d0                                                           p7509390
      fbf=z                                                             p7509400
      fpf=z                                                             p7509410
      n=ic+1                                                            p7509420
      do 1 i=1,id                                                       p7509430
      w=sa(i)*t+ddr(i)                                                  p7509440
      if (i.ge.n) go to 2                                               p7509450
      x=bf(i)*dcos(w)                                                   p7509460
      fbf=fbf+x                                                         p7509470
    2 y=pf(i)*dsin(w)                                                   p7509480
      fpf=fpf+y                                                         p7509490
    1 continue                                                          p7509500
      return                                                            p7509510
      end                                                               p7509520
      subroutine fdpf1(t,id,dpf,sa,ddr,fdpf)                            p7509530
      implicit double precision (a-h,p-z)                               p7509540
      dimension dpf(id),sa(id),ddr(id)                                  p7509550
       common kktilde
c                                                                       p7509560
c calcul somme termes precession fixe                                   p7509570
c                                                                       p7509580
      fdpf=0.0d0                                                        p7509590
      do 1 i=1,id                                                       p7509600
      yy=dpf(i)*dcos(sa(i)*t+ddr(i))                                    p7509610
      fdpf=fdpf+yy                                                      p7509620
    1 continue                                                          p7509630
      return                                                            p7509640
      end                                                               p7509650
      subroutine dff(h,prm,tset,ib,bb,br,dr,sbfh,sbfk,sbfa,spfh,spfk,   p7509660
     *spfa,sdpfh,sdpfk,sdpfa,ia,aa,ar,cr,apo,al)                        p7509670
      implicit double precision (a-h,p-z)                               p7509680
      dimension bb(ib),br(ib),dr(ib),aa(ia),ar(ia),cr(ia)               p7509690
       common kktilde
c                                                                       p7509700
c constantes                                                            p7509710
c                                                                       p7509720
      z=0.0d0                                                           p7509730
      un=1.0d0                                                          p7509740
      de=2.0d0                                                          p7509750
      ud=0.5d0                                                          p7509760
      t=dtan(h)                                                         p7509770
      tt=t*t                                                            p7509780
      tw=1.0d0/dtan(h)                                                  p7509790
      tww=tw*tw                                                         p7509800
      tut=tt+un                                                         p7509810
      twt=-tww-un                                                       p7509820
      sbfh=un                                                           p7509830
      sbfk=z                                                            p7509840
      sbfa=z                                                            p7509850
      spfh=z                                                            p7509860
      spfk=z                                                            p7509870
      spfa=un                                                           p7509880
      sdpfh=z                                                           p7509890
      sdpfk=un                                                          p7509900
      sdpfa=z                                                           p7509910
      do 1 i=1,ib                                                       p7509920
      dw=dr(i)+tset                                                     p7509930
      dww=de*dw                                                         p7509940
      bc=bb(i)*bb(i)                                                    p7509950
      st1=bb(i)*dsin(dw)                                                p7509960
      st2=bc*dsin(dww)                                                  p7509970
      ct1=bb(i)*dcos(dw)                                                p7509980
      ct2=bc*dcos(dww)                                                  p7509990
c
      xw=br(i)+prm                                                      p7510010
      ai=prm/xw                                                         p7510020
      ! MC dans sommes au bas de la p. 113 dans la these
      ! MC les ai sont les Ci dans la these
      x=ai*ai                                                           p7510030
      aii=0.25d0*t*ai*(x-un)+0.25d0*tw*x                                p7510040
      bi=t*ai*(ai-un)+tw*ai                                             p7510050
      bii=0.25d0*ai*(x+ai-un)+0.125d0*tt*x*(ai-un)*(ai-un)+0.5d0*tww*x  p7510060
      ci=ai-un                                                          p7510070
      cii=aii-ud*bi+0.25d0*tw                                           p7510080
      gi=bi-tw                                                          p7510090
      gii=bii-(ai-ud)*tww-ud*(ai*ai-ud)                                 p7510100
c                                                                       p7510110
      aiih=0.25d0*ai*(x-un)*tut+0.25d0*x*twt                            p7510120
      bih=ai*(ai-un)*tut+ai*twt                                         p7510130
      ciih=aiih-ud*bih+0.25d0*twt                                       p7510140
      sbfh=sbfh-ciih*ct2                                                p7510150
      hhh=ud*ai*(ai-un)*tut+0.25d0*twt*(de*ai-un)                       p7510160
      sbfh=sbfh-hhh*bc                                                  p7510170
      aik=br(i)/xw/xw                                                   p7510180
      bik=aik*(t*(de*ai-un)+tw)                                         p7510190
      cik=aik                                                           p7510200
      wa=de*ai*aik                                                      p7510210
      aiik=aik*(0.25d0*t*(3.0d0*x-un)+tw*de*ai*0.25d0)                  p7510220
      ciik=aiik-ud*bik                                                  p7510230
      sbfk=sbfk-cik*ct1-ciik*ct2                                        p7510240
      hhk=ud*t*(de*ai*aik-aik)+0.25d0*tw*de*aik                         p7510250
      sbfk=sbfk-hhk*bc                                                  p7510260
      sbfa=sbfa+ci*st1+cii*de*st2                                       p7510270
c                                                                       p7510280
      gih=bih-twt                                                       p7510290
      biih=0.125d0*x*(ai-un)*(ai-un)*de*t*tut+x*twt*tw                  p7510300
      giih=biih-(ai-ud)*de*tw*twt                                       p7510310
      spfh=spfh+gih*st1+st2*giih                                        p7510320
      gik=bik                                                           p7510330
      biik=aik*(0.25d0*(3.0d0*ai*ai+de*ai-un)+0.125d0*tt*de*(x-ai)*(de*ap7510340
     *i-un)+tww*ai)                                                     p7510350
      giik=biik-tww*aik-ai*aik                                          p7510360
      spfk=spfk+gik*st1+giik*st2                                        p7510370
      spfa=spfa+gi*ct1+de*gii*ct2                                       p7510380
c                                                                       p7510390
c
c prn .eq. pprn
c
      if (kktilde.eq.1)then
       sdpfh=sdpfh+gih*xw*ct1+giih*de*xw*ct2
       sdpfk=sdpfk+gik*xw*ct1+gi*ct1+giik*de*xw*ct2+de*gii*ct2          p7510410
          endif
c
c prn .ne. pprn
c     MC this is equation 66.5 p. 114 of the thesis
c     tt is tan2 epsilon; t is tan epslion. epsilon is "h" is this code
c     and given the value of 23.4 (hardwired, se above)
c
c     x is ci2

      if (kktilde.eq.0)then
       sdpfh=sdpfh+gih*xw*ct1+giih*de*xw*ct2
     $    +prm*bc*ai*(ai-un)*tut*t
c    $    +al*dsin(h)*(bc*(1.5d0+0.75d0*x-2.5d0*ai-0.5d0*ai*(ai-un)*tt))
       sdpfk=sdpfk+gik*xw*ct1+gi*ct1+giik*de*xw*ct2+de*gii*ct2          p7510410
     $  -bc*(1.5d0+0.75d0*x-2.5d0*ai-0.5d0*ai*(ai-un)*tt)
     $  -prm*bc*(1.5d0*ai-2.5d0-0.5d0*(2.0d0*ai-un)*tt)*aik 
          endif
c
      sdpfa=sdpfa-gi*xw*st1-gii*un*de*de*xw*st2                         p7510420
      if (i.eq.ib) go to 1                                              p7510430
      ii=i+1                                                            p7510440
      do 2 j=ii,ib                                                      p7510450
c                                                                       p7510460
      dv=dr(j)+tset                                                     p7510470
      dwv=dw+dv                                                         p7510480
      bd=bb(i)*bb(j)                                                    p7510490
      st3=bd*dsin(dwv)                                                  p7510500
      st4=bd*dsin(dr(i)-dr(j))                                          p7510510
      ct3=bd*dcos(dwv)                                                  p7510520
      ct4=bd*dcos(dr(i)-dr(j))                                          p7510530
      yw=br(j)+prm                                                      p7510540
      xyw=xw+yw                                                         p7510550
      bijt=br(i)-br(j)                                                  p7510560
      aj=prm/yw                                                         p7510570
      bj=t*aj*(aj-un)+tw*aj                                             p7510580
      y=aj*aj                                                           p7510590
      txz=ai*(ai-un)+aj*(aj-un)                                         p7510600
      txy=x-y+de*aj-de*ai                                               p7510610
      tyz=x+y+ai+aj-ai*aj-de                                            p7510620
      xt=0.5d0*prm/xyw                                                  p7510630
      yt=0.5d0*prm/bijt                                                 p7510640
      aij=xt*(t*(x+y-de)+tw*(ai+aj))                                    p7510650
      cij=aij-ud*(bi+bj)+ud*tw                                          p7510660
      aaij=yt*(t*txy+tw*(ai-aj))                                        p7510670
      ccij=aaij+ud*(bi+bj)-ud*tw                                        p7510680
      bij=xt*de*(0.5d0*tyz+aij*t-0.5d0*tt*txz+tww*(ai+aj))              p7510690
      bbij=yt*de*(0.5d0*(5.0d0*(ai+aj)-(x+y)-ai*aj-6.0d0)+aaij*t        p7510700
     1+0.5d0*tt*txz)                                                    p7510710
      gij=bij-(ai+aj-un)*tww-ud*(ai*ai+aj*aj-un)                        p7510720
      ggij=bbij-ud*(ai*ai-aj*aj-de*ai+de*aj)                            p7510730
c                                                                       p7510740
      aijh=xt*((x+y-de)*tut+(ai+aj)*twt)                                p7510750
      bjh=aj*(aj-un)*tut+aj*twt                                         p7510760
      cijh=aijh-ud*(bih+bjh)-ud*twt                                     p7510770
      aaijh=yt*(tut*txy+(ai-aj)*twt)                                    p7510780
      ccijh=aaijh+ud*(bih+bjh)-ud*twt                                   p7510790
      sbfh=sbfh-cijh*ct3-ccijh*ct4                                      p7510800
      ajk=br(j)/yw/yw                                                   p7510810
      wb=de*aj*ajk                                                      p7510820
      aijk=0.5d0*((t*(x+y-de)+tw*(ai+aj))*(br(i)+br(j))/xyw/xyw+prm/xyw*p7510830
     *(t*(wa+wb)+tw*(aik+ajk)))                                         p7510840
      bjk=ajk*(t*(de*aj-un)+tw)                                         p7510850
      cijk=aijk-ud*(bik+bjk)                                            p7510860
      aaijk=0.5d0*((t*txy+tw*(ai-aj))/bijt+prm/bijt*(t*(wa-wb-de*(aik-ajp7510870
     *k))+tw*(aik-ajk)))                                                p7510880
      ccijk=aaijk+ud*(bik+bjk)                                          p7510890
      sbfk=sbfk-cijk*ct3-ccijk*ct4                                      p7510900
      sbfa=sbfa+de*cij*st3                                              p7510910
c                                                                       p7510920
      bijh=xt*de*(aij*tut-txz*t*tut+t*aijh+(ai+aj)*de*tw*twt)           p7510930
      gijh=bijh-(ai+aj-un)*de*tw*twt                                    p7510940
      bbijh=yt*de*(aaij*tut+t*tut*txz+t*aaijh)                          p7510950
      ggijh=bbijh                                                       p7510960
      b2ijh=0.0d0                                                       p7510970
      g2ijh=b2ijh                                                       p7510980
      spfh=spfh+gijh*st3+ggijh*st4                                      p7510990
      uvw=(br(i)+br(j))/xyw/xyw                                         p7511000
      bijk=uvw*(bij/xt/de)+xt*de*(0.5d0*(wa+wb+aik+ajk                  p7511010
     *-ai*ajk-aik*aj)+t*aijk-0.5d0*tt*(wa+wb-aik-ajk)+tww*(aik+ajk))    p7511020
      gijk=bijk-tww*(aik+ajk)-ai*aik-aj*ajk                             p7511030
      bbijk=bbij/prm+de*yt*(0.5d0*(5.0d0*(aik+ajk)-(wa+wb)-aik*aj-ajk*aip7511040
     *)+t*aaijk+0.5d0*tt*(wa+wb-aik-ajk))                               p7511050
      ggijk=bbijk-(ai*aik-aj*ajk-aik+ajk)                               p7511060
      spfk=spfk+gijk*st3+ggijk*st4                                      p7511070
      spfa=spfa+gij*de*ct3                                              p7511080
      sdpfh=sdpfh+gijh*xyw*ct3+ggijh*bijt*ct4                           p7511090
      sdpfk=sdpfk+gijk*xyw*ct3+de*gij*ct3+ggijk*bijt*ct4                p7511100
      sdpfa=sdpfa-gij*xyw*de*st3                                        p7511110
    2 continue                                                          p7511120
    1 continue                                                          p7511130
c                                                                       p7511140
      tpl=3.0d0*apo/al                                                  p7511150
      iai=ia-1                                                          p7511160
      do 3 i=1,iai                                                      p7511170
      ii=i+1                                                            p7511180
      do 4 j=ii,ia                                                      p7511190
      b2ijk=tpl/(ar(i)-ar(j))                                           p7511200
      g2ijk=b2ijk                                                       p7511210
      tab=aa(i)*aa(j)*g2ijk                                             p7511220
      spfk=spfk+tab*dsin(cr(i)-cr(j))                                   p7511230
c
c  prm .eq. pprm
c
      if (kktilde.eq.1)then
       sdpfk=sdpfk+tab*dcos(cr(i)-cr(j))*(ar(i)-ar(j))
         endif
c
c  prm .ne. pprm
c
      if (kktilde.eq.0)then
      sdpfk=sdpfk+tab*dcos(cr(i)-cr(j))*(ar(i)-ar(j))
     $     -tpl*aa(i)*aa(j)*dcos(cr(i)-cr(j))                           p7511240
c     sdpfh=sdpfh+al*dsin(h)*tpl*aa(i)*aa(j)*dcos(cr(i)-cr(j))          p7511240
          endif
c
    4 continue                                                          p7511250
    3 continue                                                          p7511260
      return                                                            p7511270
      end                                                               p7511280
      subroutine conpr(pprm,h,al,apo,prm,ib,bb,br,ia,aa,cr)             p7511290
      implicit double precision (a-h,p-z)                               p7511300
      dimension bb(ib),br(ib),aa(ia),cr(ia)                             p7511310
       common kktilde
      x=prm                                                             p7511320
      tlp=3.0d0*apo/al                                                  p7511330
      th=dtan(h)                                                        p7511340
      thh=th*th                                                         p7511350
      ta=al*dcos(h)                                                     p7511360
      tc=0.5d0*thh-2.5d0                                                p7511370
      td=0.75d0-0.50d0*thh                                              p7511380
c
c  prm .eq. pprm
c
      if (kktilde.eq.1)then
    6 continue                                                          p7511390
      fk=0.0d0                                                          p7511400
      fki=0.0d0                                                         p7511410
      do 1 i=1,ib                                                       p7511420
      tb=br(i)+x                                                        p7511430
      t=x/tb                                                            p7511440
      ti=br(i)/tb/tb                                                    p7511450
      tt=t*t                                                            p7511460
      fkt=1.5d0+t*tc+tt*td                                              p7511470
      fkti=tc+2.0d0*t*td                                                p7511480
      fk=fk+fkt*bb(i)*bb(i)                                             p7511490
      fki=fki+fkti*bb(i)*bb(i)*ti                                       p7511500
    1 continue                                                          p7511510
      fkk=0.0d0                                                         p7511520
      ia1=ia-1                                                          p7511530
      do 2 i=1,ia1                                                      p7511540
      ii=i+1                                                            p7511550
      do 3 j=ii,ia                                                      p7511560
      t=aa(i)*aa(j)*dcos(cr(i)-cr(j))                                   p7511570
      fkk=fkk+t                                                         p7511580
    3 continue                                                          p7511590
    2 continue                                                          p7511600
      fk=fk+tlp*fkk                                                     p7511610
      fk=(1.0d0-fk)*ta-x                                                p7511620
      fkp=-1.0d0-ta*fki                                                 p7511630
      x1=x-fk/fkp                                                       p7511640
      ax=dabs(x1-x)                                                     p7511650
      if (ax-1.0d-16) 4,4,5                                             p7511660
    5 x=x1                                                              p7511670
      go to 6                                                           p7511680
    4 pprm=x                                                            p7511690
       endif
c
c  prm .ne. pprm
c
      if (kktilde.eq.0)then
      fk=0.0d0                                                          p7511400
      fki=0.0d0                                                         p7511410
      do 11 i=1,ib                                                      p7511420
      tb=br(i)+x                                                        p7511430
      t=x/tb                                                            p7511440
      ti=br(i)/tb/tb                                                    p7511450
      tt=t*t                                                            p7511460
      fkt=1.5d0+t*tc+tt*td                                              p7511470
      fkti=tc+2.0d0*t*td                                                p7511480
      fk=fk+fkt*bb(i)*bb(i)                                             p7511490
      fki=fki+fkti*bb(i)*bb(i)*ti                                       p7511500
   11 continue                                                          p7511510
      fkk=0.0d0                                                         p7511520
      ia1=ia-1                                                          p7511530
      do 12 i=1,ia1                                                     p7511540
      ii=i+1                                                            p7511550
      do 13 j=ii,ia                                                     p7511560
      t=aa(i)*aa(j)*dcos(cr(i)-cr(j))                                   p7511570
      fkk=fkk+t                                                         p7511580
   13 continue                                                          p7511590
   12 continue                                                          p7511600
      fk=fk+tlp*fkk                                                     p7511610
      fk=(1.0d0-fk)*ta-x                                                p7511620
      pprm=fk+x
         endif
c
      return                                                            p7511700
      end                                                               p7511710
      subroutine simq(a,b,n,ks)                                         p7511720
c                                                                       p7511730
c simq from the scientific subroutine package **************************p7511740
c                                                                       p7511750
c     ..................................................................p7511760
c                                                                       p7511770
c        subroutine simq                                                p7511780
c                                                                       p7511790
c        purpose                                                        p7511800
c           obtain solution of a set of simultaneous linear equations,  p7511810
c           ax=b                                                        p7511820
c                                                                       p7511830
c        usage                                                          p7511840
c           call simq(a,b,n,ks)                                         p7511850
c                                                                       p7511860
c        description of parameters                                      p7511870
c           a - matrix of coefficients stored columnwise.  these are    p7511880
c               destroyed in the computation.  the size of matrix a is  p7511890
c               n by n.                                                 p7511900
c           b - vector of original constants (length n). these are      p7511910
c               replaced by final solution values, vector x.            p7511920
c           n - number of equations and variables. n must be .gt. one.  p7511930
c           ks - output digit                                           p7511940
c                0 for a normal solution                                p7511950
c                1 for a singular set of equations                      p7511960
c                                                                       p7511970
c        remarks                                                        p7511980
c           matrix a must be general.                                   p7511990
c           if matrix is singular , solution values are meaningless.    p7512000
c           an alternative solution may be obtained by using matrix     p7512010
c           inversion (minv) and matrix product (gmprd).                p7512020
c                                                                       p7512030
c        subroutines and function subprograms required                  p7512040
c           none                                                        p7512050
c                                                                       p7512060
c        method                                                         p7512070
c           method of solution is by elimination using largest pivotal  p7512080
c           divisor. each stage of elimination consists of interchangingp7512090
c           rows when necessary to avoid division by zero or small      p7512100
c           elements.                                                   p7512110
c           the forward solution to obtain variable n is done in        p7512120
c           n stages. the back solution for the other variables is      p7512130
c           calculated by successive substitutions. final solution      p7512140
c           values are developed in vector b, with variable 1 in b(1),  p7512150
c           variable 2 in b(2),........, variable n in b(n).            p7512160
c           if no pivot can be found exceeding a tolerance of 0.0,      p7512170
c           the matrix is considered singular and ks is set to 1. this  p7512180
c           tolerance can be modified by replacing the first statement. p7512190
c                                                                       p7512200
c     ..................................................................p7512210
c                                                                       p7512220
      implicit double precision(a-h,p-z)                                p7512230
      dimension a(9),b(3)                                               p7512240
c                                                                       p7512250
c        forward solution                                               p7512260
c                                                                       p7512270
      tol=0.0d0                                                         p7512280
      ks=0                                                              p7512290
      jj=-n                                                             p7512300
      do 65 j=1,n                                                       p7512310
      jy=j+1                                                            p7512320
      jj=jj+n+1                                                         p7512330
      biga=0.0d0                                                        p7512340
      it=jj-j                                                           p7512350
      do 30 i=j,n                                                       p7512360
c                                                                       p7512370
c        search for maximum coefficient in column                       p7512380
c                                                                       p7512390
      ij=it+i                                                           p7512400
      if(dabs(biga)-dabs(a(ij))) 20,30,30                               p7512410
   20 biga=a(ij)                                                        p7512420
      imax=i                                                            p7512430
   30 continue                                                          p7512440
c                                                                       p7512450
c        test for pivot less than tolerance (singular matrix)           p7512460
c                                                                       p7512470
      if(dabs(biga)-tol) 35,35,40                                       p7512480
   35 ks=1                                                              p7512490
      return                                                            p7512500
c                                                                       p7512510
c        interchange rows if necessary                                  p7512520
c                                                                       p7512530
   40 i1=j+n*(j-2)                                                      p7512540
      it=imax-j                                                         p7512550
      do 50 k=j,n                                                       p7512560
      i1=i1+n                                                           p7512570
      i2=i1+it                                                          p7512580
      save=a(i1)                                                        p7512590
      a(i1)=a(i2)                                                       p7512600
      a(i2)=save                                                        p7512610
c                                                                       p7512620
c        divide equation by leading coefficient                         p7512630
c                                                                       p7512640
   50 a(i1)=a(i1)/biga                                                  p7512650
      save=b(imax)                                                      p7512660
      b(imax)=b(j)                                                      p7512670
      b(j)=save/biga                                                    p7512680
c                                                                       p7512690
c        eliminate next variable                                        p7512700
c                                                                       p7512710
      if(j-n) 55,70,55                                                  p7512720
   55 iqs=n*(j-1)                                                       p7512730
      do 65 ix=jy,n                                                     p7512740
      ixj=iqs+ix                                                        p7512750
      it=j-ix                                                           p7512760
      do 60 jx=jy,n                                                     p7512770
      ixjx=n*(jx-1)+ix                                                  p7512780
      jjx=ixjx+it                                                       p7512790
   60 a(ixjx)=a(ixjx)-(a(ixj)*a(jjx))                                   p7512800
   65 b(ix)=b(ix)-(b(j)*a(ixj))                                         p7512810
c                                                                       p7512820
c        back solution                                                  p7512830
c                                                                       p7512840
   70 ny=n-1                                                            p7512850
      it=n*n                                                            p7512860
      do 80 j=1,ny                                                      p7512870
      ia=it-j                                                           p7512880
      ib=n-j                                                            p7512890
      ic=n                                                              p7512900
      do 80 k=1,j                                                       p7512910
      b(ib)=b(ib)-a(ia)*b(ic)                                           p7512920
      ia=ia-n                                                           p7512930
   80 ic=ic-1                                                           p7512940
      return                                                            p7512950
        end                                                             p7512960
      subroutine phase(wq)                                              p7512970
      double precision z,wq                                             p7512980
      z=360.0d0                                                         p7512990
    4 if(wq.lt.z) go to 3                                               p7513000
      wq=wq-z                                                           p7513010
      go to 4                                                           p7513020
    3 if(wq) 5,6,6                                                      p7513030
    5 wq=wq+z                                                           p7513040
      go to 3                                                           p7513050
    6 return                                                            p7513060
      end                                                               p7513070
      subroutine awri(x,y,z,j,i,kk)                                     p7513080
      double precision x,y,z,aper,zz,a                                  p7513090
      a=1296.0d03                                                       p7513100
      zz=0.0d0                                                          p7513110
      if(y.eq.zz) go to 8                                               p7513120
      aper=a/dabs(y)                                                    p7513130
      go to 9                                                           p7513140
    8 aper=zz                                                           p7513150
    9 continue
*      write(6,6000) kk,j,i,x,y,z,aper                                  p7513160
!*      write(7,6000) kk,j,i,x,y,z,aper                                  p7513170
 6000 format(1x,i4,i3,i3,5x,f10.7,3x,f13.6,3x,f13.6,3x,f13.0)           p7513180
      return                                                            p7513190
      end                                                               p7513200
      subroutine excen(es,ec,pi,zer,ee,peva,pir)                        p7513210
      double precision es,ec,pi,zer,ee,peva,pir                         p7513220
      if(ec) 29,25,30                                                   p7513230
   29 pev=datan(es/ec)+pi                                               p7513240
      go to 24                                                          p7513250
   30 pev=datan(es/ec)                                                  p7513260
      if(es.ge.zer) go to 24                                            p7513270
      pev=2.0d0*pi+pev                                                  p7513280
      go to 24                                                          p7513290
   25 if(es) 26,27,28                                                   p7513300
   26 pev=pi+pi/2.0d0                                                   p7513310
      go to 24                                                          p7513320
   28 pev=pi/2.0d0                                                      p7513330
      go to 24                                                          p7513340
   27 pev=zer                                                           p7513350
   24 ee=dsqrt(ec*ec+es*es)                                             p7513360
      peva=pev/pir                                                      p7513370
      return                                                            p7513380
c     end                                                               p7513390
      end  subroutine 


