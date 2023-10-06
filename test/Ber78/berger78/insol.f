
      IMPLICIT REAL*8(A-H,P-Z)                                          PAL00010
      DIMENSION AE(19),BE(19),CE(19),AOB(47),BOB(47),COB(47),AOP(78),BOPPAL00020
     *(78),COP(78)                                                      PAL00030
C                                                                       PAL00040
C                                                                       PAL00050
C     THIS SOLUTION OF BERGER 1978 IS VALID ONLY FOR 1.000.000 YEARS    PAL00060
C     CENTERED ON PRESENT-DAY.                                          PAL00070
C     FOR LONGER PERIOD THE SOLUTION 1990 MUST BE USED.                 PAL00080
C     (CONTACT BERGER FOR THIS 1990 SOLUTION)                           PAL00090
C                                                                       PAL00100
C                                                                       PAL00110
C   PLEASE REFER TO :                                                   PAL00120
C      BERGER A. 1978. A SIMPLE ALGORITHM TO COMPUTE LONG TERM          PAL00130
C                      VARIATIONS OF DAILY OR MONTHLY INSOLATION        PAL00140
C                      CONTR. 18  INST OF ASTRONOMY AND GEOPHYSICS      PAL00150
C                      UNIVERSITE CATHOLIQUE DE LOUVAIN.                PAL00160
C                      LOUVAIN-LA-NEUVE    BELGIUM.                     PAL00170
C                                                                       PAL00180
C      BERGER A. 1978. LONG TERM VARIATIONS OF DAILY INSOLATION AND     PAL00190
C                      QUATERNARY CLIMATIC CHANGES                      PAL00200
C                      J. OF ATMOSPHERIC SCIENCES  35  2362-2367        PAL00210
C                                                                       PAL00220
C   THE READ AND WRITE STATEMENTS MIGHT HAVE TO BE CHANGED.             PAL00230
C   THE FUNCTION VALUE RETURNED BY DATAN IS ASSUMED TO BE A REAL*8      PAL00240
C    RANGING FROM -PI/2 TO PI/2                                         PAL00250
C                                                                       PAL00260
C                                                                       PAL00270
C                                                                       PAL00280
C   THE INPUT DATA ARE GIVEN IN ANNEX                                   PAL00290
C                                                                       PAL00300
C*******************************************                            PAL00310
C   DAILY INSOLATION - LONG TERM VARIATION *                            PAL00320
C*******************************************                            PAL00330
C                                                                       PAL00340
C   CONSTANT                                                            PAL00350
C                                                                       PAL00360
      PI=3.14159265358979D0                                             PAL00370
      PIR=PI/180.0D0                                                    PAL00380
      PIRR=PIR/3600.0D0                                                 PAL00390
      STEP=360.0D0/365.25D0                                             PAL00400
      TEST=0.0001D0                                                     PAL00410
c                                                                       PAL00411
      WRITE(6,5996)                                                     PAL00412
 5996 FORMAT(1X,'THIS PROGRAMME COMPUTES THE TOTAL DAILY IRRADIATION RECPAL00413
     .EIVED AT THE TOP '/'OF THE ATMOSPHERE FOR A GIVEN LATITUDE AND TIMPAL00414
     .E IN THE YEAR (IN KJ M-2)')                                       PAL00415
      OPEN(UNIT=8,STATUS='OLD',FILE='INSOL.IN')                         PAL00416
C                                                                       PAL00420
C   1.EARTH ORBITAL ELEMENTS : ECCENTRICITY           ECC   TABLE 1     PAL00430
C***************************   PRECESSIONAL PARAMETER PRE               PAL00440
C                              OBLIQUITY              XOB   TABLE 2     PAL00450
C                              GENERAL PRECESSION     PRG               PAL00460
C                              LONGITUDE PERIHELION   PERH  TABLE 3     PAL00470
C                                                                       PAL00480
C   READ NAME OF DATA                                                   PAL00490
C                                                                       PAL00500
      READ(8,5001)                                                      PAL00510
5001  FORMAT(/////)                                                     PAL00520
C         READ AMPLITUDE A  MEAN RATE B  PHASE C                        PAL00530
C              THEY ARE IMMEDIATELY CONVERTED IN RADIANS                PAL00540
C                                                                       PAL00550
C         NEF  NOB  NOP  MAY BE REDUCED TO  19  18  9                   PAL00560
C         BUT THE INPUT DATA MUST BE CHANGED ACCORDINGLY                PAL00570
C                                                                       PAL00580
C   ECCENTRICITY                                                        PAL00590
C                                                                       PAL00600
      NEF=19                                                            PAL00610
      DO 1 I=1,NEF                                                      PAL00620
      READ(8,5000) AE(I),Y,Z                                            PAL00630
 5000 FORMAT (13X,F11.8,F20.7,F20.6)                                    PAL00640
      BE(I)=Y*PIRR                                                      PAL00650
      CE(I)=Z*PIR                                                       PAL00660
    1 CONTINUE                                                          PAL00670
C                                                                       PAL00680
C   OBLIQUITY                                                           PAL00690
C                                                                       PAL00700
      XOD=23.320556D0                                                   PAL00710
      NOB=47                                                            PAL00720
      DO 3 I=1,NOB                                                      PAL00730
      READ(8,5002) AOB(I),Y,Z                                           PAL00740
 5002 FORMAT(7X,F13.7,2X,F10.6,2X,F10.4)                                PAL00750
      BOB(I)=Y*PIRR                                                     PAL00760
      COB(I)=Z*PIR                                                      PAL00770
    3 CONTINUE                                                          PAL00780
C                                                                       PAL00790
C   GENERAL PRECESSION IN LONGITUDE                                     PAL00800
C                                                                       PAL00810
      XOP=3.392506D0                                                    PAL00820
      PRM=50.439273D0                                                   PAL00830
      NOP=78                                                            PAL00840
      DO 31 I=1,NOP                                                     PAL00850
      READ(8,5002) AOP(I),Y,Z                                           PAL00860
      BOP(I)=Y*PIRR                                                     PAL00870
      COP(I)=Z*PIR                                                      PAL00880
   31 CONTINUE                                                          PAL00890
C                                                                       PAL00900
  100 CONTINUE                                                          PAL00910
C                                                                       PAL00920
C    2.INPUT PARAMETERS : LATITUDE PHI - TIME T - OPTION IOPT           PAL00930
C**********************                                                 PAL00940
C         IF IOPT=1  TRUE LONG.SUN TLS                                  PAL00950
C         IF IOPT=2  MONTH MA - DAY JA                                  PAL00960
C                                                                       PAL00970
C      NEFF  NOBB  NOPP  NUMBER OF TERMS KEPT FOR COMPUTATION OF        PAL00980
C                         EARTH ORBITAL ELEMENTS                        PAL00990
C      THEY CAN BE REDUCED TO 19,18,9 RESPECTIVELY                      PAL01000
C                 FOR A MINIMUM ACCURACY                                PAL01010
C                                                                       PAL01020
      WRITE(6,5999)                                                     PAL01721
 5999 FORMAT(1X,'ENTER SUCCESSIVELY : '/ 'PHI : LATITUDE'/'T   : TIME (IPAL01722
     .N 1000 YEAR'/'      0 FOR THE PRESENT'/'       NEGATIVE IN THE PASPAL01723
     .T'/'       POSITIVE IN THE FUTURE'/                               PAL01724
     .'IOPT: = 1 TRUE LONG. SUN TLS'/'      = 2  MONTH MA - DAY J'/     PAL01725
     .'NEFF  NOBB  NOPP :  NUMBER  OF TERMS KEPT FOR COMPUTATION        PAL01726
     .OF'/'THE EARTH ORBITAL ELEMENTS.'/'THEY ARE 19 47 78 AT THE MAXIMUPAL01727
     .M'/'THEY CAN BE REDUCED TO 19 18 9 RESPECTIVELY FOR A MINIMUM ACCUPAL01728
     .RACY')                                                            PAL01729
C                                                                       PAL01030
      READ(5,*) PHI,T,IOPT,NEFF,NOBB,NOPP                               PAL01035
C                                                                       PAL01040
C   3.NUMERICAL VALUE FOR ECC PRE XOB                                   PAL01050
C************************************                                   PAL01060
C       T IS NEGATIVE FOR THE PAST   T IS IN 1000 YEARS                 PAL01070
C                                                                       PAL01080
      T=T*1000.0D0                                                      PAL01090
      XES=0.0D0                                                         PAL01100
      XEC=0.0D0                                                         PAL01110
      DO 4 I=1,NEFF                                                     PAL01120
      ARG=BE(I)*T+CE(I)                                                 PAL01130
      XES=XES+AE(I)*DSIN(ARG)                                           PAL01140
      XEC=XEC+AE(I)*DCOS(ARG)                                           PAL01150
    4 CONTINUE                                                          PAL01160
      ECC=DSQRT(XES*XES+XEC*XEC)                                        PAL01170
      TRA=DABS(XEC)                                                     PAL01180
      IF(TRA.LE.1.0D-08) GO TO 10                                       PAL01190
      RP=DATAN(XES/XEC)                                                 PAL01200
      IF(XEC) 11,10,12                                                  PAL01210
   11 RP=RP+PI                                                          PAL01220
      GO TO 13                                                          PAL01230
   12 IF(XES) 14,13,13                                                  PAL01240
   14 RP=RP+2.0D0*PI                                                    PAL01250
      GO TO 13                                                          PAL01260
   10 IF(XES) 15,16,17                                                  PAL01270
   15 RP=1.5D0*PI                                                       PAL01280
      GO TO 13                                                          PAL01290
   16 RP=0.0D0                                                          PAL01300
      GO TO 13                                                          PAL01310
   17 RP=PI/2.0D0                                                       PAL01320
   13 PERH=RP/PIR                                                       PAL01330
C                                                                       PAL01340
      PRG=PRM*T                                                         PAL01350
      DO 5 I=1,NOP                                                      PAL01360
      ARG=BOP(I)*T+COP(I)                                               PAL01370
      PRG=PRG+AOP(I)*DSIN(ARG)                                          PAL01380
    5 CONTINUE                                                          PAL01390
      PRG=PRG/3600.0D0+XOP                                              PAL01400
      PERH=PERH+PRG                                                     PAL01410
   54 IF(PERH) 51,55,53                                                 PAL01420
   51 PERH=PERH+360.0D0                                                 PAL01430
      GO TO 54                                                          PAL01440
   53 IF(PERH.LT.360.0D0) GO TO 55                                      PAL01450
      PERH=PERH-360.0D0                                                 PAL01460
      GO TO 53                                                          PAL01470
   55 CONTINUE                                                          PAL01480
C                                                                       PAL01490
      PRE=ECC*DSIN(PERH*PIR)                                            PAL01500
C                                                                       PAL01510
      XOB=XOD                                                           PAL01520
      DO 6 I=1,NOBB                                                     PAL01530
      ARG=BOB(I)*T+COB(I)                                               PAL01540
      XOB=XOB+AOB(I)/3600.0D0*DCOS(ARG)                                 PAL01550
    6 CONTINUE                                                          PAL01560
C                                                                       PAL01570
C   4.DAILY INSOLATION                                                  PAL01580
C*********************                                                  PAL01590
C                                                                       PAL01600
C   OPTION SOLAR DATE - CALENDAR DATE                                   PAL01610
C       DAILY INSOLATION IN LY DAY(-1)    OR    KJ M(-2) DAY(-1)        PAL01620
C                  IF S0 IN LY MIN(-1)    OR    W M(-2)                 PAL01630
C                     TAU = 24*60 MIN     OR    24*60*60 SEC / 1000     PAL01640
C                                                                       PAL01650
C                        IN W M(-2)                                     PAL01660
C                  IF S0 IN W M(-2)      AND    TAU=1.0                 PAL01670
C                                                                       PAL01680
      SS=1353.0D0                                                       PAL01690
      TAU=86.4D0                                                        PAL01700
      SF=TAU*SS/PI                                                      PAL01710
      SO=DSIN(XOB*PIR)                                                  PAL01720
      XL=PERH+180.0D0                                                   PAL01730
      WRITE(6,6003) PHI                                                 PAL01740
 6003 FORMAT(1X,'LATITUDE =',F6.1)                                      PAL01750
      WRITE(6,6000) T,ECC,PRE,PERH,XOB                                  PAL01760
 6000 FORMAT(1X,'TIME =',F10.1,3X,'ECCENTRICITY =',F8.5,/,20X,'PREC. PARPAL01770
     *AM. =',F8.5,/,20X,'LONG. PERH. =',F7.1,/,20X,'OBLIQUITY =',F7.3)  PAL01780
C                                                                       PAL01790
      IF(IOPT.EQ.2) GO TO 20                                            PAL01800
C                                                                       PAL01810
C   4.1 SOLAR DATE                                                      PAL01820
C-----------------                                                      PAL01830
C       CONSTANT INCREMENT OF TRUE LONGITUDE OF SUN TLS                 PAL01840
C       ORIGIN IS VERNAL EQUINOX                                        PAL01850
C       IF TLS=I*30 MID-MONTH IS NOW AROUND 21                          PAL01860
C       TLS=0,30,...300,330 RESPECTIVELY FOR MARCH ... JANUARY, FEBRUARYPAL01870
C                                                                       PAL01880
      WRITE(6,5998)                                                     PAL01882
 5998 FORMAT(/1X,'ENTER TRUE LONGITUDE OF SUN (TLS)'/'TLS=0,30,...300,33PAL01884
     .0 RESPECTIVELY FOR MARCH ... JANUARY, FEBRUARY')                  PAL01886
C                                                                       PAL01888
      READ(5,*) TLS                                                     PAL01890
C                                                                       PAL01900
      CALL DAYINS(ECC,XL,SO,TLS,PHI,PIR,PI,TEST,SF,WW,DAYL)             PAL01910
      WRITE(6,6001) TLS                                                 PAL01920
 6001 FORMAT(1X,'TRUE LONG. SUN =',F7.1)                                PAL01930
      WRITE(6,6002) WW,DAYL                                             PAL01940
 6002 FORMAT(1X,'DAY INSOL. =',F7.0,1X,'KJ M(-2) DAY(-1) ',4X,'LENGTH DAPAL01950
     * =',F6.2,1X,'HOURS')                                              PAL01960
      GO TO 101                                                         PAL01970
C                                                                       PAL01980
   20 CONTINUE                                                          PAL01990
C                                                                       PAL02000
C   4.2 CALENDAR DATE  MA-JA                                            PAL02010
C--------------------                                                   PAL02020
C      ND  NUMBER OF THIS DAY IN A YEAR OF 365 DAYS                     PAL02030
C      XLAM = MEAN LONG. SUN FOR TRUE LONG. = 0                         PAL02040
C      DLAMM = MEAN LONG. SUN FOR MA-JA                                 PAL02050
C                                                                       PAL02052
      WRITE(6,5997)                                                     PAL02054
 5997 FORMAT(/1X,'ENTER CALENDAR DATE  (MONTH MA - DAY JA) ')           PAL02056
C                                                                       PAL02060
      READ(5,*) MA,JA                                                   PAL02070
C                                                                       PAL02080
      CALL NDAY(MA,JA,ND)                                               PAL02090
      XLLP=XL*PIR                                                       PAL02100
      XEE=ECC*ECC                                                       PAL02110
      XSE=DSQRT(1.0D0-XEE)                                              PAL02120
      XLAM=(ECC/2.0D0+ECC*XEE/8.0D0)*(1.0D0+XSE)*DSIN(XLLP)-XEE/4.0D0*(0PAL02130
     1.5D0+XSE)*DSIN(2.0D0*XLLP)+ECC*XEE/8.0D0*(1.0D0/3.0D0+XSE)*DSIN(3.PAL02140
     20D0*XLLP)                                                         PAL02150
      XLAM=2.0D0*XLAM/PIR                                               PAL02160
      DLAMM=XLAM+(ND-80)*STEP                                           PAL02170
      ANM=DLAMM-XL                                                      PAL02180
      RANM=ANM*PIR                                                      PAL02190
      XEC=XEE*ECC                                                       PAL02200
      RANV=RANM+(2.0D0*ECC-XEC/4.0D0)*DSIN(RANM)+5.0D0/4.0D0*ECC*ECC*   PAL02210
     1DSIN(2.0D0*RANM)+13.0D0/12.0D0*XEC*DSIN(3.0D0*RANM)               PAL02220
      ANV=RANV/PIR                                                      PAL02230
      TLS=ANV+XL                                                        PAL02240
      CALL DAYINS(ECC,XL,SO,TLS,PHI,PIR,PI,TEST,SF,WW,DAYL)             PAL02250
      WRITE(6,6010) MA,JA,TLS                                           PAL02260
 6010 FORMAT(1X,'MONTH =',I3,3X,'DAY =',I3,3X,'TLS =',F7.1)             PAL02270
      WRITE(6,6002) WW,DAYL                                             PAL02280
C                                                                       PAL02290
  101 WRITE(6,6100)                                                     PAL02300
 6100 FORMAT(/70('*')/)                                                 PAL02310
C                                                                       PAL02320
      GO TO 100                                                         PAL02330
      STOP                                                              PAL02340
      END                                                               PAL02350
      SUBROUTINE NDAY(MA,JA,ND)                                         PAL02360
      IMPLICIT REAL*8(A-H,P-Z)                                          PAL02370
      DIMENSION NJM(12)                                                 PAL02380
      DATA NJM/31,28,31,30,31,30,31,31,30,31,30,31/                     PAL02390
      ND=0                                                              PAL02400
      M=MA-1                                                            PAL02410
      IF (M.EQ.0) GO TO 2                                               PAL02420
      DO 1 I=1,M                                                        PAL02430
    1 ND=ND+NJM(I)                                                      PAL02440
    2 ND=ND+JA                                                          PAL02450
      RETURN                                                            PAL02460
      END                                                               PAL02470
      SUBROUTINE DAYINS(ECC,XL,SO,DLAM,PHI,PIR,PI,TEST,SF,WW,DAYL)      PAL02480
      IMPLICIT REAL*8(A-H,P-Z)                                          PAL02490
C                                                                       PAL02500
C   OUTPUT : WW=LY/DAY  OR  KJ M(-2) DAY(-1)  DAYL=LENGTH OF DAY (HOURS)PAL02510
C                                                                       PAL02520
      RPHI=PHI*PIR                                                      PAL02530
      RANV=(DLAM-XL)*PIR                                                PAL02540
      RAU=(1.0D0-ECC*ECC)/(1.0D0+ECC*DCOS(RANV))                        PAL02550
      S=SF/RAU/RAU                                                      PAL02560
      RLAM=DLAM*PIR                                                     PAL02570
      SD=SO*DSIN(RLAM)                                                  PAL02580
      CD=DSQRT(1.0D0-SD*SD)                                             PAL02590
      RDELTA=DATAN(SD/CD)                                               PAL02600
      DELTA=RDELTA/PIR                                                  PAL02610
      SP=SD*DSIN(RPHI)                                                  PAL02620
      CP=CD*DCOS(RPHI)                                                  PAL02630
      APHI=DABS(PHI)                                                    PAL02640
      ADELTA=DABS(DELTA)                                                PAL02650
C                                                                       PAL02660
C   SINGULARITY FOR APHI=90 AND DELTA=0                                 PAL02670
C   PARTICULAR CASES FOR PHI=0  OR  DELTA=0                             PAL02680
C                                                                       PAL02690
      TT=DABS(APHI-90.0D0)                                              PAL02700
      IF ((TT.LE.TEST).AND.(ADELTA.LE.TEST)) GO TO 2                    PAL02710
      IF(ADELTA.LE.TEST) GO TO 6                                        PAL02720
      IF(APHI.LE.TEST) GO TO 7                                          PAL02730
C                                                                       PAL02740
C   LABEL 2 : POLAR CONTINUAL NIGHT OR W=0  DAYL=0                      PAL02750
C   LABEL 4 : POLAR CONTINUAL DAY                                       PAL02760
C   LABEL 3 : DAILY SUNRISE AND SUNSET                                  PAL02770
C   LABEL 6 : EQUINOXES                                                 PAL02780
C   LABEL 7 : EQUATOR                                                   PAL02790
C                                                                       PAL02800
      AT=90.0D0-ADELTA                                                  PAL02810
      SPD=PHI*DELTA                                                     PAL02820
      IF (APHI.LE.AT) GO TO 3                                           PAL02830
      IF (SPD) 2,3,4                                                    PAL02840
    2 DAYL=0.00D0                                                       PAL02850
      WW=0.00D0                                                         PAL02860
      GO TO 5                                                           PAL02870
    4 DAYL=24.00D0                                                      PAL02880
      WW=S*SP*PI                                                        PAL02890
      GO TO 5                                                           PAL02900
    3 TP=-SP/CP                                                         PAL02910
      STP=DSQRT(1.0D0-TP*TP)                                            PAL02920
      RDAYL=DACOS(TP)                                                   PAL02930
      DAYL=24.0D0*RDAYL/PI                                              PAL02940
      WW=S*(RDAYL*SP+CP*STP)                                            PAL02950
      GO TO 5                                                           PAL02960
    6 DAYL=12.0D0                                                       PAL02970
      WW=S*DCOS(RPHI)                                                   PAL02980
      GO TO 5                                                           PAL02990
    7 DAYL=12.0D0                                                       PAL03000
      WW=S*DCOS(RDELTA)                                                 PAL03010
    5 RETURN                                                            PAL03020
      END                                                               PAL03030
