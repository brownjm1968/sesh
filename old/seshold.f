CC ---------------------------------------------------------------------    
CC Modified for AIX xlf Fortran compiler  ( Klaus Berthold, IRMM 4/94)
CC ---------------------------------------------------------------------    
CC
CC(Revised and extended by Konstantin Volev, IRMM 5/2002 )
CC
CC MAIN :  107 FORMAT Statements  (no more using Hollerith count)
CC         109 
CC         113  
CC         114
CC         115
CC         116
CC
CC         COMMON BLOCK : CHARACTER*26 CDATE
CC 
CC SUBROUTINE MUSC :  3 CONTINUE (one removed)
CC                    3 CONTINUE 
CC 
CC SUBROUTINE MOCT :    PLOT(....) commented out
CC  
CC FUNCTION RHO    :    CBRT not intrinsic so added FUNCTION CBRT
CC
CC SUBROUTINE SPACE:    RANDU not intrinsic so added SUBROUTINE RANDU 
CC 
CC ---------------------------------------------------------------------    
C     SESH MAIN PROGRAM (KFK VERSION 1975, EXPERIMENTAL COPY 1994)
C
C     MANUAL: F.H. FROEHNER,
C             "SESH - A FORTRAN IV CODE FOR CALCULATING THE SELF-
C             SHIELDING AND MULTIPLE SCATTERING EFFECTS FOR
C             NEUTRON CROSS SECTION DATA INTERPRETATION
C             IN THE UNRESOLVED RESONANCE REGION",
C             REPORT GA-8380 (1968)
C
C     GILBERT-CAMERON COMPOSITE LEVEL DENSITY AND GIANT DIPOLE RESONANCE
C     MODEL ARE USED FOR THE ENERGY DEPENDENCE OF LEVEL SPACINGS AND
C     RADIATION WIDTHS IN VERSION 1975. THE INPUT REMAINS AS DESCRIBED
C     IN GA-8380 WITH TWO EXCEPTIONS: (1) THE NUCLEAR TEMPERATURE IS
C     REPLACED BY THE GILBERT-CAMERON PAIRING ENERGY OF THE COMPOUND
C     NUCLEUS, (2) ONLY THE S-WAVE LEVEL SPACING MUST BE GIVEN (SPACING
C     INPUT FOR L>0 IS IGNORED). FURTHERMORE, CHANNEL RADIUS AND
C     EFFECTIVE NUCLEAR RADIUS ARE DISTINCT FOR ALL PARTIAL WAVES:
C     THE CHANNEL RADIUS IS TAKEN AS RC(I)=(1.23*A**(1/3)+0.80) FM,
C     THE EFFECTIVE NUCLEAR RADII R(L,I) ARE INPUT NUMBERS WITH THE
C     DEFAULT VALUES R(L,I)=RC(I).
C
C     INPUT LIMITATIONS:
C             UP TO      10 ISOTOPES
C             "  "        4 PARTIAL WAVES
C             "  "        5 SAMPLE GEOMETRIES
C             "  "       24 NEUTRON ENERGIES
C             "  "   100000 MONTE CARLO HISTORIES PER ENERGY
C
      COMMON COMM(18),ZL(4),SGI(10,100),SGL(10,100,4),AP(10),UM(10),
     1A(11),AB(11),BE(11),PE(11),TEFF(11),SPIN(11),NL(10),AA(10),RC(10),
     2GG(5,11),D(5,11),S(5,11),SI(5,11),R(5,11),EFF(5,11),J2X(4,10),
     3J2N(4,10),SUMGJ(4,10),XN(6),E(100),SG(100),SP(100),SC(100),ST(100)
     4,G(8,4,10),GNR(8,4,10),GNRIN(8,4,10),DJL(8,4,10),SSCF,DSSCF,N,K,
     5I,L,J,NN,NE,NI,LX,RN(6),RX(6),ZH(100),PO,PS,DPO,DPS,ZP,SGM,STM,
     6DSGM,DSTM,TMC,DTMC,TAU,DTAU,SNC,XNSS,ITYPE,IHIST
     7,LJX(2,8,10),LJN(2,8,10),JX2(10),JN2(10),JMX(10)
      COMMON/RANDM/ IY
      DIMENSION STI(10,100),SPI(10,100),STL(10,100,4),SPL(10,100,4)
c      CHARACTER*26 CDATE
c FG051198
      CHARACTER*30 CDATE
c modification for HP on bonsai3
      open (unit=5,file='input.data',status='old')
      open (unit=8,file='output.data',status='unknown')
c end modification for HP on bonsai3
C
C             READ INPUT
C             IHIST IS A TAG FOR HISTOGRAMS OF ST, SG, T, F0;
C             IHIST=1  MEANS HISTOGRAM VALUES ARE CALCULATED AND PRINTED
C             IHIST=0  MEANS NO SUCH VALUES ARE CALCULATED OR PRINTED
C     CALL ERRSET(217,0  ,-1,1,0,100)
C     CALL ERRSET(208,300,-1,0,1,208)
      IY=4751
    1 READ(5,100,END=999)COMM,IHIST
  100 FORMAT(18A4,1I8)
      SUM=0.
      I=0
    2 I=I+1
C             I-TH ISOTOPE (MAXIMUM OF 10 ISOTOPES):
C             NUCLEON NUMBER A(I), ABUNDANCE AB(I), BINDING ENERGY BE(I) [MeV], PAIRING
C             ENERGY PE(I) [MeV], EFFECTIVE SAMPLE TEMPERATURE  TEFF(I) [DEG. K],
C             TARGET SPIN QUANTUM NUMBER SPIN(I);
      READ(5,101)A(I),AB(I),BE(I),PE(I),TEFF(I),SPIN(I)
  101 FORMAT(6E10.5)
      L=0
    3 L=L+1
C             L-TH PARTIAL WAVE (UP TO F-WAVE):
C             RADIATION WIDTH GG(L,I) [eV], AVERAGE LEVEL SPACING D(L,I) [eV], STRENGTH
C             FUNCTION FOR ELASTIC SCATTERING S(L,I), STRENGTH FUNCTION FOR
C             INELASTIC SCATTERING SI(L,I), NUCLEAR RADIUS R(L,I) [fm], DETECTION
C             EFFICIENCY
      READ(5,102)GG(L,I),D(L,I),S(L,I),SI (L,I),R(L,I),EFF(L,I)
  102 FORMAT(6E10.5)
C             CHECK FOR SAMPLE THICKNESS CARD  (FIRST WORD ZERO)
      IF(GG(L,I).EQ.0.0)GO TO 4
C             CHECK FOR LAST PARTIAL WAVE
      IF(S(L,I).LT.1.)GO TO 3
C             LAST CARD WAS ISOTOPE CARD. STORE CORRESPONDINGLY
      NL(I)=L-1
      I=I+1
      A(I)=GG(L,I-1)
      AB(I)=D(L,I-1)
      BE(I)=S(L,I-1)
      PE(I)=SI(L,I-1)
      TEFF(I)=R(L,I-1)
      SPIN(I)=EFF(L,I-1)
      GG(L,I-1)=0.
      D(L,I-1) =0.
      S(L,I-1) =0.
      SI(L,I-1)=0.
      R(L,I-1) =0.
      EFF(L,I-1)=0.
      L=0
      GO TO 3

C             LAST CARD WAS SAMPLE THICKNESS CARD (NUCLEI/B),
C             STORE CORRESPONDINGLY.

    4 NI=I
      NL(I)=L-1
      XN(1)=D(L,I)
      XN(2)=S(L,I)
      XN(3)=SI(L,I)
      XN(4)=R(L,I)
      XN(5)=EFF(L,I)
      XN(6)=0.

C             FIND NUMBER OF SAMPLE THICKNESSES
C             FOR TRANSMISSION AND CAPTURE DATA

      DO 5 N=2,6
      IF(XN(N).NE.0.)GO TO 5
      NN=N-1
      GO TO 6
    5 CONTINUE

C             READ OUTER RADII (NUCLEI/BARN)

    6 READ(5,103)(RX(N),N=1,5)
      RX(6)=0.
      IF(NN.GT.1)GO TO 8

C             FIND NUMBER OF SHELL THICKNESSES FOR SHELL TRANSMISSION
C             DATA

      DO 7 N=2,6
      IF(RX(N).NE.0.)GO TO 7
      NN=N-1
      GO TO 8
    7 CONTINUE


C             READ INNER RADII (NUCLEI/BARN) FOR SHELL TRANSMISSION
C             CALCULATIONS. RN(N)=0. MEANS CYLINDRICAL SAMPLE


    8 READ(5,103)(RN(N),N=1,5)

  103 FORMAT(E20.5,4E10.5)


C             READ NUMBER OF RESONANCE PAIRS


      READ(5,104)ZP
  104 FORMAT(E10.5)
      DO 9 M=1,25
      M4=(M-1)*4

C             READ ENERGIES (KEV) AND NUMBERS OF MONTE CARLO HISTORIES

      READ(5,105)E(M4+1),ZH(M4+1),E(M4+2),ZH(M4+2)
     *          ,E(M4+3),ZH(M4+3),E(M4+4),ZH(M4+4)
  105 FORMAT(8E10.5)


C             BLANK CARD SIGNALS END OF INPUT

      IF(E(M4+1).EQ.0.)GO TO 10
    9 CONTINUE


C             FIND NUMBER OF ENERGIES, NE

   10 DO 11 J=1,4
      MJ=M4-4+J
      IF(E(MJ).GT.0.)NE=MJ
   11 CONTINUE

C             WRITE INPUT HEADING
CC    VMS SPECIFIC DATE AND TIME PROCEDURES
c      CALL FDATE_(CDATE)
c FG051198
c modification for HP on bonsai3
c     CALL FDATE(CDATE)
c     WRITE(8,106)CDATE,COMM
c 106 FORMAT(1H1,//,'[UNIX-SESH ',A24,']     ',18A4)
c end modification for HP on bonsai3
      WRITE(8,107)
  107 FORMAT(1H0,//,
     *' I N P U T '/
     *' ========= '//
     *' NUCLEON   ABUN-    BINDING  PAIRING EFF.      NUCL.  ',
     *'ORB.ANG.  AV. RAD.    AV. LEVEL   STRENGTH    STRENGTH',
     *'    NUCLEAR     EFFI- '/
     *' NUMBER    DANCE    ENERGY   ENERGY  TEMP.     SPIN   ',
     *' MOM.     WIDTH       SPACING     FCT. FOR    FCT. FOR',
     *'    RADIUS      CIENCY'/
     *'                    (MEV)    (MEV)   (DEG.K)   Q.NO.  ',
     *' Q.NO.    (EV)        (EV)        EL. SCATT.  INEL. SC',
     *'.   (FM)              '/)


C             PREPARE ENERGY-INDEPENDENT PARAMETERS AND PRINT INPUT

      DO 16 I=1,NI

C             CHANNEL RADIUS

      RC(I)=1.23*A(I)**(1./3.)+0.8

C             GILBERT-CAMERON MATCHING ENERGY

      AI=FLOAT(INT(A(I)+1.5))
      write (6,*)'a= ',a(i),ai
      UM(I)=2.5+150./AI

C             FIND FERMI GAS MODEL A-PARAMETER FROM S-WAVE SPACING

      QI=SPIN(I)+.5
      U=BE(I)-PE(I)
      CC=ALOG(0.2367E6*AI*U/(D(1,I)*QI))
      XO=CC
      DO 12 M=1,12
      VARJ2=0.1460*XO*AI**(2./3.)
      FJ=.5*(EXP(-(QI-1.)/VARJ2)-EXP(-(QI+1.)/VARJ2))*VARJ2/QI
      XX=CC-ALOG(FJ)+2.*ALOG(XO)
      IF(ABS(XX-XO).LT.1.E-6)GO TO 13
      XO=XX
   12 CONTINUE
   13 AP(I)=XX**2/(4.*U)

C             DETERMINE MINIMUM AND MAXIMUM COMPOUND SPIN POSSIBLE

      X2I=2.*SPIN(I)
      I2=X2I+.01
      LX=NL(I)
      DO 15 L=1,LX
      J2X(L,I)=I2+2*L-1
      J2N(L,I)=1
      IF(I2.GT.2*L-2)J2N(L,I)=I2-2*L+1
      IF(I2.LT.2*L-2)J2N(L,I)=2*L-I2-3

C             CALCULATE SUM OF SPIN FACTORS

      SUMGJ(L,I)=FLOAT((J2X(L,I)+J2N(L,I)+2)*(J2X(L,I)-J2N(L,I)+2))/8./
     1(X2I+1.)

C             LEVEL DENSITY FOR GIVEN L

      IF(L.GT.1)D(L,I)=D(1,I)/SUMGJ(L,I)
      J=0
      J2MN=J2N(L,I)
      J2MX=J2X(L,I)
      DO 14 J2=J2MN,J2MX,2
      J=J+1

C             CALCULATE SPIN FACTOR

      X2J=FLOAT(J2)
      G(J,L,I)=(X2J+1.)/2./(X2I+1.)

C             LEVEL DENSITY FOR GIVEN L AND J

      DJL(J,L,I)=D(1,I)/G(J,L,I)
      ZL(L)=FLOAT(L-1)
   14 CONTINUE
      IF(R(L,I).EQ.0.)R(L,I)=RC(I)
   15 CONTINUE
      WRITE(8,108)A(I),AB(I),BE(I),PE(I),TEFF(I),SPIN(I),(ZL(L),GG(L,I)
     1,D(L,I),S(L,I),SI (L,I),R(L,I),EFF(L,I),L=1,LX)
  108 FORMAT(F6.1,F11.4,F8.3,F9.3,2F8.1,F7.0,1PE16.4,1P,4E12.4,0PF9.4/
     1                           (50X,0PF7.0,1PE16.4,1P,4E12.4,0PF9.4))
   16 CONTINUE


C             FIRST PART
C             ANALYTICAL CROSS SECTION CALCULATION FOR ALL ENERGIES
C             BEGIN K-LOOP (ENERGIES)

      DO 21 K=1,NE
      SG(K)=0.
      SC(K)=0.
      SP(K)=0.
      ST(K)=0.

C             BEGIN ISOTOPE LOOP

      DO 20 I=1,NI
      SGI(I,K)=0.
      STI(I,K)=0.
      SPI(I,K)=0.
      DO 17 L=1,4
   17 SGL(I,K,L)=0.
      AA(I)=(1.+1./A(I))**2

C             GET ENERGY DEPENDENCE FACTORS FOR LEVEL SPACINGS AND
C             RADIATION WIDTHS
C

         CALL ENDEP(E(K),A(I),BE(I),PE(I),AP(I),UM(I),EDGG,EDD)

C
C             BEGIN LOOP OF PARTIAL WAVES
      LX=NL(I)
      DO 19 L=1,LX

C             GET PENETRABILITIES, SHIFTS, HARD SPHERE PHASE SHIFTS
      AK=RC(I)*SQRT(E(K)/AA(I))/143.92

C

         CALL PEPS(AK,L,DUMMY,VL)

C

      RK=AK*R(L,I)/RC(I)

C

         CALL PEPS(RK,L,PL,DUMMY)

C
C             CALCULATE SUM OVER ALL POSSIBLE COMPOUND SPINS

      S3=0.
      J=0
      J2MN=J2N(L,I)
      J2MX=J2X(L,I)
      DO 18 J2=J2MN,J2MX,2
      J=J+1

C             FIND NUMBER OF POSSIBLE CHANNEL SPINS FOR GIVEN J AND L
      EJL=1.
      IF(J2.LT.J2X(L,I).AND.J2.GT.J2N(L,I).OR.2*L.EQ.I2+2.AND.L.NE.1.AND
     1.J2.EQ.J2N(L,I))EJL=2.

C             REDUCED WIDTHS FOR ELASTIC AND INELASTIC SCATTERING
      GNR(J,L,I)  =S(L,I) *EJL*DJL(J,L,I)
      GNRIN(J,L,I)=SI(L,I)*EJL*DJL(J,L,I)

C             NEUTRON WIDTHS FOR ELASTIC AND INELASTIC SCATTERING

      SQ=SQRT(E(K)*1000.)*VL
      GN  =GNR  (J,L,I)*SQ
      GNIN=GNRIN(J,L,I)*SQ

C             FIND PSI-FUNCTION GIVING PORTER-THOMAS AVERAGE

      ETA=SQRT((GG(L,I)*EDGG+GNIN)/(2.*GN*EDD))

C

         CALL PFCN(0.,ETA,U,V,KZ)

C

      PSI0=1.77245  *ETA*U

C             FORM J-SUM (COMPOUND SPINS)
      S3=S3+G(J,L,I)**2*(1.-PSI0)
   18 CONTINUE

C             CONTRIBUTION TO EFFECTIVE CAPTURE CROSS SECTION
      SGL(I,K,L)=4.09E3*AA(I)*AB(I)*GG(L,I)*EDGG*S3/E(K)/D(L,I)/EDD/SUMG
     1J(L,I)*EFF(L,I)

C             PURE CAPTURE CROSS SECTION

      SC(K)=SC(K)+SGL(I,K,L)/EFF(L,I)

C             CONTRIBUTION TO POTENTIAL SCATTERING CROSS SECTION

      SPL(I,K,L)=2.605E3/E(K)*AA(I)*AB(I)*(2.*FLOAT(L)-1.)*SIN(PL)**2

C             CONTRIBUTION TO TOTAL CROSS SECTION

      STL(I,K,L)=SPL(I,K,L)+4.09E6/SQRT(1000.*E(K))*AA(I)*AB(I)*S(L,I)*
     1(2.*FLOAT(L)-1.)*VL*COS(2.*PL)

C             FORM L-SUMS (PARTIAL WAVES)

      SGI(I,K)=SGI(I,K)+SGL(I,K,L)
      SPI(I,K)=SPI(I,K)+SPL(I,K,L)
      STI(I,K)=STI(I,K)+STL(I,K,L)
   19 CONTINUE

C             FORM I-SUMS (ISOTOPES)
C             AVERAGE EFFECTIVE CAPTURE CROSS SECTION, K-TH ENERGY

      SG(K)=SG(K)+SGI(I,K)

C             AVERAGE TOTAL CROSS SECTION, K-TH ENERGY

      ST(K)=ST(K)+STI(I,K)


C             POTENTIAL SCATTERING CROSS SECTION, K-TH ENERGY

      SP(K)=SP(K)+SPI(I,K)
   20 CONTINUE
   21 CONTINUE
      WRITE(8,109)
  109 FORMAT(1H0,//,
     *' A N A L Y T I C A L   R E S U L T S '/
     *' =================================== '//
     *' AVERAGE CROSS SECTIONS                               ',
     *'ISOTOPIC CONTRIBUTIONS                    PARTIAL-WAVE',
     *' CONTRIBUTIONS      '//
     *' NEUTRON   TOTAL    POTENTIAL     PURE    EFFECTIVE   ',
     *'NUCLEON  TOTAL    POTENTIAL  EFFECTIVE    L    TOTAL  ',
     *'  POTENTIAL  EFFECTIVE'/
     *'  ENERGY           SCATTERING    CAPTURE   CAPTURE    ',
     *'NUMBER           SCATTERING   CAPTURE                 ',
     *' SCATTERING   CAPTURE'/
     *'  (KEV)    (BARN)     (BARN)     (BARN)    (BARN)     ',
     *'         (BARN)     (BARN)    (BARN)           (BARN) ',
     *'    (BARN)    (BARN) ')
      DO 25 K=1,NE
      WRITE(8,110)E(K),ST(K),SP(K),SC(K),SG(K),A(1),STI(1,K),SPI(1,K),SG
     1I(1,K),STL(1,K,1),SPL(1,K,1),SGL(1,K,1)
  110 FORMAT(/,0PF7.1,1P,4E11.3,0PF7.0,1P,3E11.3,4H   0,1P,3E11.3)
      DO 24 I=1,NI
      IF(I.EQ.1)GO TO 22
      WRITE(8,111)                             A(I),STI(I,K),SPI(I,K),SG
     1I(I,K),STL(I,K,1),SPL(I,K,1),SGL(I,K,1)
  111 FORMAT(/,051X           ,0PF7.0,1P,3E11.3,4H   0,1P,3E11.3)
   22 IF(NL(I).EQ.1)GO TO 24
      LX=NL(I)
      DO 23 L=2,LX
      LL=L-1
      WRITE(8,112)LL,STL(I,K,L),SPL(I,K,L),SGL(I,K,L)
  112 FORMAT(94X,I1,1P,3E11.3)
   23 CONTINUE
   24 CONTINUE
   25 CONTINUE
C             END OF ANALYTICAL CALCULATION
      IF(ZH(1).LT.1.)GO TO 1
C             SECOND PART
C             MONTE CARLO CALCULATION OF CROSS SECTIONS AND OF
C             PROBABILITY FOR DETECTED CAPTURE
      IF(XN(1).GT.0.0.AND.RX(1).GT.0.0.AND.RN(1).EQ.0.0)ITYPO=1
      IF(XN(1).EQ.0.0.AND.RX(1).GT.0.0.AND.RN(1).GT.0.0)ITYPO=2
      IF(XN(1).GT.0.0.AND.RX(1).EQ.0.0.AND.RN(1).EQ.0.0)ITYPO=3
      IF(XN(1).GT.0.0.AND.RX(1).GT.0.0.AND.RN(1).GT.0.0)ITYPO=4
      IF(ITYPO.EQ.1)WRITE(8,113)
      IF(ITYPO.EQ.2)WRITE(8,114)
      IF(ITYPO.EQ.3)WRITE(8,115)
      IF(ITYPO.EQ.4)WRITE(8,116)
  113 FORMAT(1H1,//,
     *' M O N T E   C A R L O   R E S U L T S '/
     *' ===================================== '//
     *' CAPTURE DATA, CIRCULAR DISC SAMPLE    '//
     *' SAMPLE    NEUTRON  AV. EFF.  AVERAGE   FIRST-COLL.  P',
     *'ROB. FOR   SELF-SHIELD.  AVERAGE     NUMBER OF  NUMBER',
     *' OF  SAMPLE  AV.   '/
     *' THICK-    ENERGY   CAPTURE   TOTAL     CAPTURE      C',
     *'APTURE     CORRECTION    NUMBER OF   HISTORIES  RESONA',
     *'NCE  RADIUS  TRANS-'/
     *' NESS               CROSS     CROSS     PROBAB./     A',
     *'FTER 1     FACTOR/       COLLISIONS             PAIRS ',
     *'             MISS./'/
     *'                    SECTION/  SECTION/  UNCERT.      C',
     *'OLLISION/  UNCERT.                                    ',
     *'             UNCERT.'/
     *'                    UNCERT.   UNCERT.                U',
     *'NCERT. '/
     *' (NUC./B)   (KEV)   (B)/(B)   (B)/(B)                 ',
     *'             /(PERCENT)                               ',
     *'   (NUC./B)' )
  114 FORMAT(1H1,//,
     *' M O N T E   C A R L O   R E S U L T S '/
     *' ===================================== '//
     *' CAPTURE DATA, SPHERICAL SHELL'//
     *' SHELL     NEUTRON  AV. EFF.  AVERAGE   FIRST-COLL.  P',
     *'ROB. FOR   SHELL         AVERAGE     NUMBER OF  NUMBER',
     *' OF  SAMPLE  AV.    '/
     *' THICK-    ENERGY   CAPTURE   TOTAL     CAPTURE      C',
     *'APTURE     TRANS-        NUMBER OF   HISTORIES  RESONA',
     *'NCE  RADIUS  TRANS- '/
     *' NESS               CROSS     CROSS     PROBAB./     A',
     *'FTER 1     MISSION/      COLLISIONS             PAIRS ',
     *'             MISS./ '/
     *'                    SECTION/  SECTION/  UNCERT.      C',
     *'OLLISION/  UNCERT.                                    ',
     *'             UNCERT.'/
     *'                    UNCERT.   UNCERT.                U',
     *'NCERT. '/
     *' (NUC./B)   (KEV)   (B)/(B)   (B)/(B)                 ',
     *'                           (NUC./B)' )
  115 FORMAT(1H1,//,
     *' M O N T E   C A R L O   R E S U L T S '/
     *' ===================================== '//
     *' TRANSMISSION DATA '//
     *' SAMPLE    NEUTRON  AV. EFF.  AVERAGE   AVERAGE   SELF',
     *'-SHIELD.  NUMBER OF  NUMBER OF '/
     *' THICK-    ENERGY   CAPTURE   TOTAL     TRANS.    CORR',
     *'ECTION    HISTORIES  RESONANCE '/
     *' NESS               CROSS     CROSS     MISSION/  FACT',
     *'OR/                  PAIRS     '/
     *'                    SECTION/  SECTION/  UNCERT.   UNCE',
     *'RT. '/
     *'                    UNCERT.   UNCERT.                 ',
     *'    '/
     *' (NUC./B)   (KEV)   (B)/(B)   (B)/(B)               /(',
     *'PERCENT) ')
  116 FORMAT(1H1,//,
     *' M O N T E   C A R L O   R E S U L T S '/
     *' ===================================== '//
     *' SELF-INDICATION DATA, DISC SAMPLES '//
     *' CAPTURE   NEUTRON  AV. EFF.  AVERAGE   FIRST-COLL.  P',
     *'ROB. FOR   SELF-SHIELD.  AVERAGE     NUMBER OF  NUMBER',
     *' OF  SAMPLE  FIRST '/
     *' SAMPLE    ENERGY   CAPTURE   TOTAL     CAPTURE      C',
     *'APTURE     CORRECTION    NUMBER OF   HISTORIES  RESONA',
     *'NCE  RADIUS  SAMPLE'/
     *' THICK.             CROSS     CROSS     PROBAB./     A',
     *'FTER 1     FACTOR/       COLLISIONS             PAIRS ',
     *'             THICK.'/
     *'                    SECTION/  SECTION/  UNCERT.      C',
     *'OLLISION/  UNCERT. '/
     *'                    UNCERT.   UNCERT.                U',
     *'NCERT. '/
     *' (NUC./B)   (KEV)   (B)/(B)   (B)/(B)                 ',
     *'           /(PERCENT)                                 ',
     *'    (NUC./B) (NUC./B)')
C             DETERMINATION OF EXTREME VALUES OF L FOR EACH PARITY AND J
      DO 36 I=1,NI
      LX=NL(I)
C             TWICE EXTREME COMPOUND SPIN QUANTUM NUMBERS
      I2=2.*SPIN(I)
      JX2(I)=I2+2*LX-1
      JN2(I)=I2-2*LX+1
      IF(I2.GT.6)GO TO 26
      IF(MOD(I2,2).EQ.0)JN2(I)=1
      IF(MOD(I2,2).EQ.1)JN2(I)=0
C             JMX(I): NUMBER OF POSSIBLE J VALUES
   26 JMX(I)=(JX2(I)-JN2(I))/2+1
      JX=JMX(I)
C             J: COUNTER FOR COMPOUND SPIN IN ASCENDING ORDER
      DO 35 J=1,JX
      J2=JN2(I)+2*J-2
C             M: PARITY LABEL, M=1 FOR PARITY OF TARGET GROUND STATE
      M=1
C             EXTREMAL L VALUES FOR GIVEN PARITY, J, AND ISOTOPE
   27 IF(LX.GT.M+1)GO TO 30
      IF(J2N(M,I).LE.J2.AND.J2X(M,I).GE.J2)GO TO 29
   28 LJX(M,J,I)=-1
      LJN(M,J,I)=-1
      GO TO 34
   29 LJX(M,J,I)=M
      LJN(M,J,I)=M
      GO TO 34
   30 IF(J2N(M  ,I).LE.J2.AND.J2X(M  ,I).GE.J2)GO TO 31
      IF(J2N(M+2,I).LE.J2.AND.J2X(M+2,I).GE.J2)GO TO 33
      GO TO 28
   31 IF(J2N(M+2,I).LE.J2.AND.J2X(M+2,I).GE.J2)GO TO 32
      GO TO 29
   32 LJX(M,J,I)=M+2
      LJN(M,J,I)=M
      GO TO 34
   33 LJX(M,J,I)=M+2
      LJN(M,J,I)=M+2
   34 IF(M.GE.2)GO TO 35
      M=M+1
      IF(M.LE.LX)GO TO 27
      LJX(M,J,I)=-1
      LJN(M,J,I)=-1
   35 CONTINUE
   36 CONTINUE
C             BEGIN N-LOOP (SAMPLE THICKNESSES)
   37 DO 46 NQ=1,NN
      N=NQ
C             BEGIN K-LOOP (ENERGIES)
      DO 45 KQ=1,NE
      K=KQ
      IF(ZH(KQ).EQ.0.)GO TO 45
      IF(ZH(KQ).GT.0..AND.ZH(KQ).LT.100.)ZH(KQ)=100.
      IF(XN(N).GT.0.0 .AND. RX(N).NE.0.0 .AND. RN(N).EQ.0.0)ITYPE=1
      IF(XN(N).EQ.0.0 .AND. RX(N).NE.0.0 .AND. RN(N).NE.0.0)ITYPE=2
      IF(XN(N).GT.0.0 .AND. RX(N).EQ.0.0 .AND. RN(N).EQ.0.0)ITYPE=3
      IF(XN(N).GT.0.0 .AND. RX(N).NE.0.0 .AND. RN(N).NE.0.0)ITYPE=4
      IF(ITYPE.EQ.ITYPO)GO TO 38
      IF(ITYPE.EQ.1)WRITE(8,113)
      IF(ITYPE.EQ.2)WRITE(8,114)
      IF(ITYPE.EQ.3)WRITE(8,115)
      IF(ITYPE.EQ.4)WRITE(8,116)
      ITYPO=ITYPE
   38 CONTINUE
      GO TO (39,40,42,43),ITYPE
C
   39    CALL MUSC
C
      GO TO 41
C
   40    CALL MUSS
C
      XN(N)=RX(N)-RN(N)
   41 IF(K.EQ.1)
     *WRITE(8,117)XN(N),E(K),SGM,STM,PO,PS,SSCF,SNC, ZH(K),ZP,RX(N),
     *TMC  ,DSGM,DSTM,DPO,DPS,DSSCF,DTMC
      IF(K.GT.1)
     *WRITE(8,118)      E(K),SGM,STM,PO,PS,SSCF,SNC, ZH(K),ZP,RX(N),
     *TMC  ,DSGM,DSTM,DPO,DPS,DSSCF,DTMC
      GO TO 44
C
   42    CALL MOCT
C
      IF(K.EQ.1)
     *WRITE(8,119)XN(N),E(K),SGM,STM,TMC,TAU,ZH(K),ZP,DSGM,DSTM,DTMC,
     *DTAU
      IF(K.GT.1)
     *WRITE(8,120)      E(K),SGM,STM,TMC,TAU,ZH(K),ZP,DSGM,DSTM,DTMC,
     *DTAU
      GO TO 44
C
   43    CALL MUSC
C
      IF(K.EQ.1)
     *WRITE(8,121)XN(N),E(K),SGM,STM,PO,PS,SSCF,SNC, ZH(K),ZP,RX(N),
     *RN(N),DSGM,DSTM,DPO,DPS,DSSCF
      IF(K.GT.1)
     *WRITE(8,122)      E(K),SGM,STM,PO,PS,SSCF,SNC, ZH(K),ZP,RX(N),
     *RN(N),DSGM,DSTM,DPO,DPS,DSSCF
  117 FORMAT(' '/
     *1PE10.3,0PF7.3,1PE11.3,2E10.3,3E13.3,0P,2F10.0,1X,1P,2E10.3/
     *           17X,1PE11.3,2E10.3,2E13.3,44X        ,1P E10.3/)
  118 FORMAT(0PF17.3,1PE11.3,2E10.3,3E13.3,0P,2F10.0,1X,1P,2E10.3/
     *           17X,1PE11.3,2E10.3,2E13.3,44X        ,1P E10.3/)
  119 FORMAT(' '/
     *1PE10.3,0PF7.3,1PE11.3,2E10.3,E13.4,0P,2F11.0/
     *           17X,1PE11.3,2E10.3,E13.4/)
  120 FORMAT(0PF17.3,1PE11.3,2E10.3,E13.4,0P,2F11.0/
     *           17X,1PE11.3,2E10.3,E13.4/)
  121 FORMAT(' '/
     *1PE10.3,0PF7.3,1PE11.3,2E10.3,3E13.3,0P,2F10.0,1X,1P,2E10.3/
     *           17X,1PE11.3,2E10.3,2E13.3/)
  122 FORMAT(0PF17.3,1PE11.3,2E10.3,3E13.3,0P,2F10.0,1X,1P,2E10.3/
     *           17X,1PE11.3,2E10.3,2E13.3/)
   44 CONTINUE
   45 CONTINUE
   46 CONTINUE
      GO TO 1
  999 STOP  
      END
C
      SUBROUTINE MUSC
C             MUSC  YIELDS THE MULTIPLE SCATTERING CORRECTION FOR A
C             CYLINDRICAL CAPTURE SAMPLE (MONTE CARLO)
C                XN(N)      CAPTURE SAMPLE THICKNESS
C                RX(N)      CAPTURE SAMPLE RADIUS
C                RN(N)      FILTER  SAMPLE THICKNESS
      COMMON COMM(18),ZL(4),SGI(10,100),SGL(10,100,4),AP(10),UM(10),
     1A(11),AB(11),BE(11),PE(11),TEFF(11),SPIN(11),NL(10),AA(10),RC(10),
     2GG(5,11),D(5,11),S(5,11),SI(5,11),R(5,11),EFF(5,11),J2X(4,10),
     3J2N(4,10),SUMGJ(4,10),XN(6),E(100),SG(100),SP(100),SS(100),STT(100
     4),G(8,4,10),GNR(8,4,10),GNRIN(8,4,10),DJL(8,4,10),SSCF,DSSCF,N,K,
     5I,L,J,NN,NE,NI,LX,RN(6),RX(6),ZH(100),PO,PS,DPO,DPS,ZP,SGM,STM,
     6DSGM,DSTM,TMC,DTMC,TAU,DTAU,SNC,XNSS,ITYPE,IHIST
     7,LJX(2,8,10),LJN(2,8,10),JX2(10),JN2(10),JMX(10)
      COMMON/A/ PSM(1000),POM(1000),DOP(10,10),PLE(10,4,10),GGE(10,4,10)
     1,DE(10,8,4,10),GNE(10,8,4,10),GNINE(10,8,4,10),EE(10),
     2STMC(1000),SGMC(1000),TM(1000)
      DIMENSION ZST(101),ZSG(101),ZT(101),ZF0(101),GNL(8),
     1COSE(10,4,10),SINE(10,4,10)
      IF(IHIST.NE.1)GO TO 51
      DO 50 IZ=1,101
      ZST(IZ)=0.
      ZSG(IZ)=0.
      ZT (IZ)=0.
   50 ZF0(IZ)=0.
   51 CONTINUE
      NS=1
      MM=1
      SNC=0.
      DTMC=0.
      SUM=0.
      SOM=0.
      SAM=0.
      SEM=0.
      SIM=0.
      SYM=0.
      ZSIR=0.
C     A1=(1.-EXP(-STT(K)*XN(N)))/STT(K)
C     A2=A1*SG(K)
C     IF(ITYPE.EQ.4)A1=A1*EXP(-STT(K)*RN(N))
      FE=1.-2./(A(1)*AA(1))
      DO 20 NC=1,10
      IF(NC.EQ.1)EE(NC)=E(K)
      IF(NC.GT.1)EE(NC)=EE(NC-1)*FE
      DO 20 I=1,NI
      DOP(NC,I)=SQRT(0.72531*A(I)/TEFF(I)/EE(NC))
C
         CALL ENDEP(EE(NC),A(I),BE(I),PE(I),AP(I),UM(I),EDGG,EDD)
C
      LX=NL(I)
      DO 20 L=1,LX
      GGE(NC,L,I)=GG(L,I)*EDGG
      AK=RC(I)*SQRT(EE(NC)/AA(I))/143.92
C
         CALL PEPS(AK,L,DUMMY,VL)
C
      RK=AK*R(L,I)/RC(I)
C
         CALL PEPS(RK,L,PLE(NC,L,I),DUMMY)
C
      SINE(NC,L,I)=SIN(2.*PLE(NC,L,I))
      COSE(NC,L,I)=COS(2.*PLE(NC,L,I))
      JX=(J2X(L,I)-J2N(L,I))/2+1
      DO 20 JL=1,JX
      J=(J2N(L,I)-JN2(I))/2+JL
      DE(NC,J,L,I)=DJL(JL,L,I)*EDD
      SQ=SQRT(EE(NC)*1000.)*VL*EDD
      GNE  (NC,J,L,I)=GNR  (JL,L,I)*SQ
      GNINE(NC,J,L,I)=GNRIN(JL,L,I)*SQ
   20 CONTINUE
C             BEGIN LOOP OF NEUTRON HISTORIES
      NH=ZH(K)
      NP=ZP
      DO 1 M=1,NH
C             INITIALIZE
      U=0.
      V=0.
      W=1.
      X=RX(N)*SQRT(RANDOM(D))
      Y=0.
      Z=0.
      TX=XN(N)
      NC=0
C             BEGIN OF COLLISION LOOP
    2 IF(NC.LT.10)NC=NC+1
C             FIND CROSS SECTIONS FOR COLLISION
      SN=0.
      SC=0.
      ST=SP(K)
      DO 3 I=1,NI
      JX=JMX(I)
      DO 83 J=1,JX
C             LP=1: PARITY OF TARGET GROUND STATE
C             LP=2: OPPOSITE PARITY
      DO 93 LP=1,2
      L1=LJN(LP,J,I)
      L2=LJX(LP,J,I)
      IF(L1.LT.0)GO TO 93
      JL=J-(J2N(L1,I)-JN2(I))/2
C        4.617E3=2.605E3*1.772454
      C0=4.617E3/EE(NC)*AA(I)*AB(I)*G(JL,L1,I)
C             MP: COUNTER FOR RESONANCE PAIRS
      MP=1
C             FIND CENTRAL LEVEL INTERVAL AND ENERGY IN IT
C
         CALL SPACE(DE(NC,J,L1,I),DD)
C
      H=DD*RANDOM(R)
C             FIND REDUCED NEUTRON WIDTH
   17 GNS=0.
      GNINS=0.
      DO 40 L=L1,L2,2
C
         CALL PORTER(GNE(NC,J,L,I),GNL(L))
C
      GNS=GNS+GNL(L)
   40 GNINS=GNINS+GNINE(NC,J,L,I)
      GT=GNS+GNINS+GGE(NC,L1,I)
      ETA=GT*DOP(NC,I)
      XI=2.*ETA*H/GT
      CT=C0*ETA/GT
      CG=CT*GNS*GGE(NC,L1,I)/GT
C
         CALL PFCN(XI,ETA,UU,VV,KZ)
C
      SCC=CG*UU
      SN=SN+SCC
C             EFFECTIVE CAPTURE CROSS SECTION
      SC=SC+SCC*EFF(L1,I)
C             TOTAL CROSS SECTION
      DO 41 L=L1,L2,2
   41 ST=ST+CT*GNL(L)*(UU*COSE(NC,L,I)+VV*SINE(NC,L,I))
      IF(MP.EQ.1)GO TO 18
C             CHECK IF SECOND MEMBER OF PAIR IS ALREADY INCLUDED
      IF(H.LT.0.)GO TO 16
C
         CALL WIGNER(DE(NC,J,L1,I),DD)
C
      HNEG=HNEG-DD
      H=HNEG
      GO TO 17
C             CHECK IF SECOND MEMBER OF PAIR IS ALREADY INCLUDED
   18 IF(H.LT.0.)GO TO 16
      HPOS=H
      HNEG=H-DD
      H=HNEG
      GO TO 17
C             INCLUDE ANOTHER PAIR OF RESONANCES IF REQUIRED
   16 IF(MP.EQ.NP)GO TO 93
      MP=MP+1
C
         CALL WIGNER(DE(NC,J,L1,I),DD)
C
      HPOS=HPOS+DD
      H=HPOS
      GO TO 17
   93 CONTINUE
   83 CONTINUE
    3 CONTINUE
      SN=ST-SN
      T=EXP(-TX*ST)
C             INTERACTING FRACTION
      WI=1.-T
      WG=WI*SC
C             SURVIVING FRACTION
      WN=WI*SN/ST
      IF(NC.GT.1)GO TO 6
      TF=EXP(-RN(N)*ST)
      WG=WG*TF
      WN=WN*TF
C             CALCULATE HISTOGRAM VALUES IF REQUIRED
      IF(IHIST.NE.1)GO TO 44
      IZ=20.*ST/STT(K)+1.
      IF(IZ.LT.1)  IZ=1
      IF(IZ.GT.101)IZ=101
      ZST(IZ)=ZST(IZ)+1.
      IZ=20.*SC/SG(K)+1.
      IF(IZ.LT.1)  IZ=1
      IF(IZ.GT.101)IZ=101
      ZSG(IZ)=ZSG(IZ)+1.
      IZ=100.*T+1.
      IF(IZ.LT.1)  IZ=1
      IF(IZ.GT.101)IZ=101
      ZT(IZ)=ZT(IZ)+1.
      IZ=20.*(WG/ST)/(XN(N)*SG(K)*EXP(-RN(N)*STT(K)))
      IF(IZ.LT.1)  IZ=1
      IF(IZ.GT.101)IZ=101
      ZF0(IZ)=ZF0(IZ)+1.
      ZSIR=ZSIR+TF*SC
   44 CONTINUE
      POMG=WG/ST
      SOM=SOM+POMG
      PSMG=0.
C
         CALL AVERT(TM(NS),T,SGMC(NS),SC,STMC(NS),ST,NS,1)
C
      GO TO 7
    6 WG=WG/ST
      SUM=SUM+WG
      PSMG=PSMG+WG
C             PLAY RUSSIAN ROULETTE WITH SURVIVAL PROBABILITY WN
    7 IF(RANDOM(X).GT.WN)GO TO 19
C             FREE PATH SAMPLING
      FP=-ALOG(1.-WI*RANDOM(X))/ST
C             COORDINATES OF COLLISION
      X=X+U*FP
      Y=Y+V*FP
      Z=Z+W*FP
C             SCATTERING ANGLES
      PHI=6.28318 *RANDOM(D)
      CTHC=2.*RANDOM(D)-1.
      CTH=(1.+A(1)*CTHC)/SQRT(1.+A(1)**2+2.*A(1)*CTHC)
      STH=SQRT(1.-CTH**2)
C             NEW DIRECTION COSINES
      RHO=SQRT(1.-W**2)
      IF(RHO.GT.0.001)GO TO 8
      U=W*STH*COS(PHI)
      V=W*STH*SIN(PHI)
      W=W*CTH
      GO TO 9
    8 O=CTH*U+STH*(U*W*COS(PHI)-V*SIN(PHI))/RHO
      V=CTH*V+STH*(V*W*COS(PHI)+U*SIN(PHI))/RHO
      W=CTH*W-STH*COS(PHI)*RHO
      U=O
C             PATH LENGTH TO SAMPLE SURFACE
    9 B1=U**2+V**2
      B2=U*X+V*Y
      B3=RX(N)**2-X**2-Y**2
      IF(B2.LT.0.)GO TO 10
      TX=B3/(SQRT(B2**2+B1*B3)+B2)
      GO TO 11
   10 TX=(SQRT(B2**2+B1*B3)-B2)/B1
   11 IF(W.LT.-0.000001)GO TO 12
      IF(W.GT.+0.000001)GO TO 13
      GO TO 2
   12 TP=-Z/W
      GO TO 14
   13 TP=(XN(N)-Z)/W
   14 IF(TP.LT.TX)TX=TP
      GO TO 2
   19 SNC=SNC+FLOAT(NC)
C
    1    CALL AVERP(POM(MM),POMG,PSM(MM),PSMG,MM,1)
C
C             END OF HISTORY LOOP
C             MONTE CARLO RESULTS
C
         CALL AVERP(POM(MM),POMG,PSM(MM),PSMG,MM,2)
         CALL AVERT(TM(NS),T,SGMC(NS),SC,STMC(NS),ST,NS,2)
C
      TMC=0.
      SGM=0.
      STM=0.
      DO 22 M=1,NS
      TMC=TMC+TM(M)
      SGM=SGM+SGMC(M)
   22 STM=STM+STMC(M)
      TMC=TMC/FLOAT(NS)
      SGM=SGM/FLOAT(NS)
      STM=STM/FLOAT(NS)
      PO =SOM/ZH(K)
      PS =SUM/ZH(K)
      SNC=SNC/ZH(K)
      ZSIR=ZSIR/(ZH(K)*SGM)
C             STATISTICAL UNCERTAINTIES
      DO 23 M=1,NS
      DTMC=DTMC+(TM(M)-TMC)**2
      SIM=SIM+(SGMC(M)-SGM)**2
   23 SYM=SYM+(STMC(M)-STM)**2
      B11=0.
      B12=0.
      B23=0.
      BBOTH=0.
      DO 15 M=1,MM
      BBOTH=BBOTH+(POM(M)+PSM(M)-PO-PS)*(SGMC(M)-SGM)
      IF(ITYPE.EQ.4)B12=B12+(POM(M)+PSM(M)-PO-PS)*(STMC(M)-STM)
      IF(ITYPE.EQ.4)B23=B23+(STMC(M)-STM)*(SGMC(M)-SGM)
                    B11=B11+(POM(M)-PO)*(PSM(M)-PS)
      SAM=SAM+(POM(M)-PO)**2
   15 SEM=SEM+(PSM(M)-PS)**2
      DYN=FLOAT(NS*(NS-1))
      DTMC=SQRT(DTMC/DYN)
      DSGM=SQRT(SIM/DYN)
      DSTM=SQRT(SYM/DYN)
      DEN=FLOAT(MM*(MM-1))
      DPO=SQRT(SAM/DEN)
      DPS=SQRT(SEM/DEN)
      BBOTH=BBOTH/DEN
      B12=B12/DEN
      B23=B23/DEN
      B11=B11/DEN
      B11=2.*B11+DPO**2+DPS**2
C             SELF-SHIELDING CORRECTION FACTOR AND ITS PERCENT ERROR
      SSCF=(PO+PS)/(XN(N)*SGM)
      IF(ITYPE.EQ.4)EXN1=EXP(RN(N)*STM)
      IF(ITYPE.EQ.4)SSCF=SSCF*EXN1
      IF(ITYPE.EQ.1)DSSCF=SQRT(B11/(XN(N)*SGM)**2-2.*BBOTH*SSCF/(XN(N)*
     1SGM**2)+(SSCF/SGM)**2*DSGM**2)/SSCF*100.0
      IF(ITYPE.EQ.4)DSSCF=SQRT((EXN1/(XN(N)*SGM))**2*B11+2.*(EXN1/(XN(N)
     1*SGM))*RN(N)*SSCF*B12-2.*(EXN1/(XN(N)*SGM))*SSCF/SGM*BBOTH+(RN(N)*
     2SSCF*DSTM)**2-2.*RN(N)*SSCF**2/SGM*B23+(SSCF/SGM*DSGM)**2)/
     3SSCF*100.0
C             WRITE HISTOGRAM VALUES IF REQUIRED
      IF(IHIST.NE.1)GO TO 53
      SUMST=0.
      SUMSG=0.
      SUMT =0.
      SUMF0=0.
      DO 54 IZ=1,101
      SUMST=SUMST+ZST(IZ)
      SUMSG=SUMSG+ZSG(IZ)
      SUMT =SUMT + ZT(IZ)
      SUMF0=SUMF0+ZST(IZ)
   54 CONTINUE
      DO 55 IZ=1,101
      ZST(IZ)=ZST(IZ)/SUMST
      ZSG(IZ)=ZSG(IZ)/SUMSG
      ZT(IZ)=ZT(IZ)/SUMT
      ZF0(IZ)=ZF0(IZ)/SUMF0
   55 CONTINUE
      WRITE(8,198)ZSIR
  198 FORMAT(/1PE15.3/)
      WRITE(8,199)(IZ,ZST(IZ),ZSG(IZ),ZT(IZ),ZF0(IZ),IZ=1,101)
  199 FORMAT(1H1,///,' NO.      RELATIVE FREQUENCIES OF       '        /
     1               '          ST     SG     T      F0       '       //
     2(I4,0PF8.4,3F7.4))
   53 CONTINUE
      RETURN
      END
C
      SUBROUTINE MUSS
C             MUSS YIELDS THE MULTIPLE SCATTERING CORRECTION FOR A
C             SPHERICAL SHELL (MONTE CARLO)
      COMMON COMM(18),ZL(4),SGI(10,100),SGL(10,100,4),AP(10),UM(10),
     1A(11),AB(11),BE(11),PE(11),TEFF(11),SPIN(11),NL(10),AA(10),RC(10),
     2GG(5,11),D(5,11),S(5,11),SI(5,11),R(5,11),EFF(5,11),J2X(4,10),
     3J2N(4,10),SUMGJ(4,10),XN(6),E(100),SG(100),SP(100),SS(100),STT(100
     4),G(8,4,10),GNR(8,4,10),GNRIN(8,4,10),DJL(8,4,10),SSCF,DSSCF,N,K,
     5I,L,J,NN,NE,NI,LX,RN(6),RX(6),ZH(100),PO,PS,DPO,DPS,ZP,SGM,STM,
     6DSGM,DSTM,TMC,DTMC,TAU,DTAU,SNC,XNSS,ITYPE,IHIST
     7,LJX(2,8,10),LJN(2,8,10),JX2(10),JN2(10),JMX(10)
      COMMON/A/ PSM(1000),POM(1000),DOP(10,10),PLE(10,4,10),GGE(10,4,10)
     1,DE(10,8,4,10),GNE(10,8,4,10),GNINE(10,8,4,10),EE(10),
     3STMC(1000),SGMC(1000),TM(1000)
      DIMENSION ZST(101),ZSG(101),ZT(101),ZF0(101),GNL(8),
     1COSE(10,4,10),SINE(10,4,10)
      IF(IHIST.NE.1)GO TO 51
      DO 50 IZ=1,101
      ZST(IZ)=0.
      ZSG(IZ)=0.
      ZT (IZ)=0.
   50 ZF0(IZ)=0.
   51 CONTINUE
      NS=1
      MM=1
      SNC=0.
      SUM=0.
      SOM=0.
      SAM=0.
      SEM=0.
      SIM=0.
      SYM=0.
      DTMC=0.
      XNSS=RX(N)-RN(N)
      A1=(1.-EXP(-STT(K)*XNSS))/STT(K)
      A2=A1*SG(K)
      IF(ITYPE.EQ.4)A1=A1*EXP(-STT(K)*RN(N))
      FE=1.-2./(A(1)*AA(1))
      DO 20 NC=1,10
      IF(NC.EQ.1)EE(NC)=E(K)
      IF(NC.GT.1)EE(NC)=EE(NC-1)*FE
      DO 20 I=1,NI
      DOP(NC,I)=SQRT(0.72531*A(I)/TEFF(I)/EE(NC))
      LX=NL(I)
      DO 20 L=1,LX
C
         CALL ENDEP(EE(NC),A(I),BE(I),PE(I),AP(I),UM(I),EDGG,EDD)
C
      GGE(NC,L,I)=GG(L,I)*EDGG
      RK=RC(I)*SQRT(EE(NC)/AA(I))/143.92
C
         CALL PEPS(RK,L,PLE(NC,L,I),VL)
C
      PLE(NC,L,I)=PLE(NC,L,I)*R(L,I)/RC(I)
      SINE(NC,L,I)=SIN(2.*PLE(NC,L,I))
      COSE(NC,L,I)=COS(2.*PLE(NC,L,I))
      JX=(J2X(L,I)-J2N(L,I))/2+1
      DO 20 JL=1,JX
      J=(J2N(L,I)-JN2(I))/2+JL
      DE(NC,J,L,I)=DJL(JL,L,I)*EDD
      SQ=SQRT(EE(NC)*1000.)*VL*EDD
      GNE  (NC,J,L,I)=GNR  (JL,L,I)*SQ
      GNINE(NC,J,L,I)=GNRIN(JL,L,I)*SQ
   20 CONTINUE
C             BEGIN LOOP OF NEUTRON HISTORIES
      NH=ZH(K)
      NP=ZP
      DO 1 M=1,NH
C             INITIALIZE
      Q=0.
      CTH=-1.
      STH=0.
      Z=RN(N)
      TX=XNSS
      W=1.
      NC=0
C             BEGIN OF COLLISION LOOP
    2 IF(NC.LT.10)NC=NC+1
C             FIND CROSS SECTIONS FOR COLLISION
      SN=0.
      SC=0.
      ST=SP(K)
      DO 3 I=1,NI
      JX=JMX(I)
      DO 83 J=1,JX
      DO 93 LP=1,2
      L1=LJN(LP,J,I)
      L2=LJX(LP,J,I)
      IF(L1.LT.0)GO TO 93
C             4.617E3=2.605E3*1.772454
      JL=J-(J2N(L1,I)-JN2(I))/2
      C0=4.617E3/EE(NC)*AA(I)*AB(I)*G(JL,L1,I)
      MP=1
C             FIND CENTRAL LEVEL INTERVAL AND ENERGY IN IT
C
         CALL SPACE(DE(NC,J,L1,I),DD)
C
      H=DD*RANDOM(R)
   17 GNS=0.
      GNINS=0.
      DO 40 L=L1,L2,2
C             FIND REDUCED NEUTRON WIDTH
C
         CALL PORTER(GNE(NC,J,L,I),GNL(L))
C
      GNS=GNS+GNL(L)
   40 GNINS=GNINS+GNINE(NC,J,L,I)
      GT=GNS+GNINS+GGE(NC,L1,I)
      ETA=GT*DOP(NC,I)
      XI=2.*ETA*H/GT
      CT=C0*ETA/GT
      CG=CT*GNS*GGE(NC,L1,I)/GT
C
         CALL PFCN(XI,ETA,UU,VV,KZ)
C
C             CAPTURE CROSS SECTION
      SCC=CG*UU
      SN=SN+SCC
C             EFFECTIVE CAPTURE CROSS SECTION
      SC=SC+SCC*EFF(L1,I)
C             TOTAL CROSS SECTION
      DO 41 L=L1,L2,2
   41 ST=ST+CT*GNL(L)*(UU*COSE(NC,L,I)+VV*SINE(NC,L,I))
      IF(MP.EQ.1)GO TO 18
C             CHECK IF SECOND MEMBER OF PAIR IS ALREADY INCLUDED
      IF(H.LT.0.)GO TO 16
C
         CALL WIGNER(DE(NC,J,L1,I),DD)
C
      HNEG=HNEG-DD
      H=HNEG
      GO TO 17
C             CHECK IF SECOND MEMBER OF PAIR IS ALREADY INCLUDED
   18 IF(H.LT.0.)GO TO 16
      HPOS=H
      HNEG=H-DD
      H=HNEG
      GO TO 17
C             INCLUDE ANOTHER PAIR OF RESONANCES IF REQUIRED
   16 IF(MP.EQ.NP)GO TO 93
      MP=MP+1
C
         CALL WIGNER(DE(NC,J,L1,I),DD)
C
      HPOS=HPOS+DD
      H=HPOS
      GO TO 17
   93 CONTINUE
   83 CONTINUE
    3 CONTINUE
      SN=ST-SN
      T=EXP(-TX*ST)
C             INTERACTING FRACTION
      WI=1.-T
C             CAPTURED AND DETECTED FRACTION
      WG=WI*SC
      IF(ITYPE.EQ.4)WG=WG*EXP(-RN(N)*ST)
C             SURVIVING FRACTION
      WN=WI*SN/ST
      IF(NC.GT.1)GO TO 6
C             CALCULATE HISTOGRAM VALUES IF REQUIRED
      IF(IHIST.NE.1)GO TO 52
      IZ=20.*ST/STT(K)+1.
      IF(IZ.LT.1)  IZ=1
      IF(IZ.GT.101)IZ=101
      ZST(IZ)=ZST(IZ)+1.
      IZ=20.*SC/SG(K)+1.
      IF(IZ.LT.1)  IZ=1
      IF(IZ.GT.101)IZ=101
      ZSG(IZ)=ZSG(IZ)+1.
      IZ=100.*T+1.
      IF(IZ.LT.1)  IZ=1
      IF(IZ.GT.101)IZ=101
      ZT(IZ)=ZT(IZ)+1.
      IZ=20.*(WG/ST)/(EXP(-XNSS*STT(K))*XNSS*SG(K))+1.
      IF(IZ.LT.1)  IZ=1
      IF(IZ.GT.101)IZ=101
      ZF0(IZ)=ZF0(IZ)+1.
   52 CONTINUE
      POMG=WG*(1./ST-A1/WI)+A2
      SOM=SOM+POMG
      PSMG=0.
C
         CALL AVERT(TM(NS),T,SGMC(NS),SC,STMC(NS),ST,NS,1)
C
      GO TO 7
    6 WG=WG/ST
      PSMG=PSMG+WG
      SUM=SUM+WG
C             PLAY RUSSIAN ROULETTE WITH SURVIVAL PROBABILITY WN
    7 IF(RANDOM(X).GT.WN)GO TO 19
C             CALCULATE RADIAL COORDINATE (Z) OF NC-TH COLLISION
      T=-ALOG(1.-WI*RANDOM(D))/ST
      IF(T.GT.(Z*CTH-Q))  T=T+2.*Q
      Z=SQRT(Z**2+  T**2-2.*Z*  T*CTH)
C             ANGLE WITH RESPECT TO RADIUS AFTER COLLISION
      CTH=2.*RANDOM(D)-1.
      STH=SQRT(1.-CTH**2)
C             MAXIMUM PATH LENGTH IN MATTER
      Q=0.
      IF(CTH.GT.SQRT(1.-(RN(N)/Z)**2))Q=SQRT(RN(N)**2-(Z*STH)**2)
      TX=Z*CTH+SQRT(RX(N)**2-(Z*STH)**2)-2.*Q
      GO TO 2
   19 SNC=SNC+FLOAT(NC)
C
    1    CALL AVERP(POM(MM),POMG,PSM(MM),PSMG,MM,1)
C
C             END OF HISTORY LOOP
C             MONTE CARLO RESULTS
C
         CALL AVERP(POM(MM),POMG,PSM(MM),PSMG,MM,2)
         CALL AVERT(TM(NS),T,SGMC(NS),SC,STMC(NS),ST,NS,2)
C
      TMC=0.
      SGM=0.
      STM=0.
      DO 22 M=1,NS
      TMC=TMC+TM(M)
      SGM   =SGM   +SGMC(M)
   22 STM   =STM   +STMC(M)
      SGM   =SGM   /FLOAT(NS)
      TMC=TMC/FLOAT(NS)
      STM   =STM   /FLOAT(NS)
      PO=SOM/ZH(K)
      PS=SUM/ZH(K)
      SNC=SNC/ZH(K)
C             STATISTICAL UNCERTAINTIES
      DO 23 M=1,NS
      DTMC=DTMC+(TM(M)-TMC)**2
      SIM=SIM+(SGMC(M)-SGM   )**2
   23 SYM=SYM+(STMC(M)-STM   )**2
      B11=0.
      B12=0.
      B23=0.
      BBOTH=0.
      DO 15 M=1,MM
      BBOTH=BBOTH+(POM(M)+PSM(M)-PO-PS)*(SGMC(M)-SGM)
                    B11=B11+(POM(M)-PO)*(PSM(M)-PS)
      SAM=SAM+(POM(M)-PO     )**2
   15 SEM=SEM+(PSM(M)-PS     )**2
      DYN=FLOAT(NS*(NS-1))
      DTMC=SQRT(DTMC/DYN)
      DSGM   =SQRT(SIM/DYN)
      DSTM   =SQRT(SYM/DYN)
      DEN=FLOAT(MM*(MM-1))
      DPO     =SQRT(SAM/DEN)
      DPS     =SQRT(SEM/DEN)
      BBOTH=BBOTH/DEN
      B12=B12/DEN
      B11=B11/DEN
      B11=2.*B11+DPO**2+DPS**2
C             SELF-SHIELDING CORRECTION FACTOR AND ITS PERCENT ERROR
      SSCF=1.-PO-PS
      DSSCF=SQRT(DPO**2+DPS**2)
C             WRITE HISTOGRAM VALUES IF REQUIRED
      IF(IHIST.NE.1)GO TO 53
      SUMST=0.
      SUMSG=0.
      SUMT =0.
      SUMF0=0.
      DO 54 IZ=1,101
      SUMST=SUMST+ZST(IZ)
      SUMSG=SUMSG+ZSG(IZ)
      SUMT =SUMT + ZT(IZ)
      SUMF0=SUMF0+ZST(IZ)
   54 CONTINUE
      DO 55 IZ=1,101
      ZST(IZ)=ZST(IZ)/SUMST
      ZSG(IZ)=ZSG(IZ)/SUMSG
      ZT(IZ)=ZT(IZ)/SUMT
      ZF0(IZ)=ZF0(IZ)/SUMF0
   55 CONTINUE
      WRITE(8,199)(IZ,ZST(IZ),ZSG(IZ),ZT(IZ),ZF0(IZ),IZ=1,101)
  199 FORMAT(1H1,///,' NO.    TALLY FOR HISTOGRAM OF        '         /
     1               '          ST     SG     T      F0     '         //
     2(I4,0PF8.4,3F7.4))
   53 CONTINUE
      RETURN
      END
C
      SUBROUTINE MOCT
C             MOCT  YIELDS THE SELF-SHIELDING CORRECTION FACTOR AND
C             TRANSMISSION FOR A CYLINDRICAL SAMPLE
      COMMON COMM(18),ZL(4),SGI(10,100),SGL(10,100,4),AP(10),UM(10),
     1A(11),AB(11),BE(11),PE(11),TEFF(11),SPIN(11),NL(10),AA(10),RC(10),
     2GG(5,11),D(5,11),S(5,11),SI(5,11),R(5,11),EFF(5,11),J2X(4,10),
     3J2N(4,10),SUMGJ(4,10),XN(6),E(100),SG(100),SP(100),SS(100),STT(100
     4),G(8,4,10),GNR(8,4,10),GNRIN(8,4,10),DJL(8,4,10),SSCF,DSSCF,N,K,
     5I,L,J,NN,NE,NI,LX,RN(6),RX(6),ZH(100),PO,PS,DPO,DPS,ZP,SGM,STM,
     6DSGM,DSTM,TMC,DTMC,TAU,DTAU,SNC,XNSS,ITYPE,IHIST
     7,LJX(2,8,10),LJN(2,8,10),JX2(10),JN2(10),JMX(10)
      COMMON /B/ TMCG(1000),SGMC(1000),STMC(1000),DOP(10),GGE(4,10),DE(8
     1,4,10),VL(4,10),PL(4,10),EDGG(10,100),EDD(10,100),GN(8,4,10),GNIN
     2(8,4,10)
      DIMENSION TN(101),ZT(101),GNL(8),SINE(4,10),COSE(4,10)
      DIMENSION ZST(101), ZSG(101)
      MM=1
      SIM=0.
      SYM=0.
      TMC=0.
      TN(1)=1.
      ZT(1)=0.
      ZST(1)=0.
      ZSG(1)=0.
      DT=0.02
      TSUM=0.
      TNUM=0.
      TSSM=0.
      DO 2 IT=2,101
      TN(IT)=TN(IT-1)+1.
      ZST(IT)=0.
      ZSG(IT)=0.
    2 ZT(IT)=0.
      DO 20 I=1,NI
      DOP(I)=SQRT(0.72531*A(I)/TEFF(I)/E(K))
C
        CALL ENDEP(E(K),A(I),BE(I),PE(I),AP(I),UM(I),EDGG(I,K),EDD(I,K))
C
      LX=NL(I)
      DO 20 L=1,LX
      GGE(L,I)=GG(L,I)*EDGG(I,K)
      AK=RC(I)*SQRT(E(K)/AA(I))/143.92
C
         CALL PEPS(AK,L,DUMMY,VL(L,I))
C
      RK=AK*R(L,I)/RC(I)
C
         CALL PEPS(RK,L,PL(L,I),DUMMY)
C
      SINE(L,I)=SIN(2.*PL(L,I))
      COSE(L,I)=COS(2.*PL(L,I))
      JX=(J2X(L,I)-J2N(L,I))/2+1
      DO 20 JL=1,JX
      J=(J2N(L,I)-JN2(I))/2+JL
      DE(J,L,I)=DJL(JL,L,I)*EDD(I,K)
      GN  (J,L,I)=GNR  (JL,L,I)*SQRT(E(K)*1000.)*VL(L,I)*EDD(I,K)
      GNIN(J,L,I)=GNRIN(JL,L,I)*SQRT(E(K)*1000.)*VL(L,I)*EDD(I,K)
   20 CONTINUE
C             BEGIN LOOP OF NEUTRON HISTORIES
      NH=ZH(K)
      NP=ZP
      DO 1 M=1,NH
      NHIST=0
      NHIST=NHIST+1
C             FIND CROSS SECTIONS FOR COLLISION
      SN=0.
      SC=0.
      ST=SP(K)
      DO 3 I=1,NI
      JX=JMX(I)
      DO 83 J=1,JX
      DO 93 LP=1,2
      L1=LJN(LP,J,I)
      L2=LJX(LP,J,I)
      IF(L1.LT.0)GO TO 93
C             4.617E3=2.605E3*1.772454
      JL=J-(J2N(L1,I)-JN2(I))/2
      C0=4.617E3/E(K)*AA(I)*AB(I)*G(JL,L1,I)
      MP=1
C             FIND CENTRAL LEVEL INTERVAL AND ENERGY IN IT
C
         CALL SPACE(DE(J,L1,I),DD)
C
      H=DD*RANDOM(R)
   17 GNS=0.
      GNINS=0.
      DO 40 L=L1,L2,2
C             FIND REDUCED NEUTRON WIDTH
C
         CALL PORTER(GN(J,L,I),GNL(L))
C
      GNS=GNS+GNL(L)
   40 GNINS=GNINS+GNIN(J,L,I)
      GT=GNS+GNINS+GGE(L1,I)
      ETA=GT*DOP(I)
      XI=2.*ETA*H/GT
      CT=C0*ETA/GT
      CG=CT*GNS*GGE(L1,I)/GT
C
         CALL PFCN(XI,ETA,UU,VV,KZ)
C
C             CAPTURE CROSS SECTION
      SCC=CG*UU
      SN=SN+SCC
C             EFFECTIVE CAPTURE CROSS SECTION
      SC=SC+SCC*EFF(L1,I)
C             TOTAL CROSS SECTION
      DO 41 L=L1,L2,2
   41 ST=ST+CT*GNL(L)*(UU*COSE(L,I)+VV*SINE(L,I))
      IF(MP.EQ.1)GO TO 18
C             CHECK IF SECOND MEMBER OF PAIR IS ALREADY INCLUDED
      IF(H.LT.0.)GO TO 16
C
         CALL WIGNER(DE(J,L1,I),DD)
C
      HNEG=HNEG-DD
      H=HNEG
      GO TO 17
C             CHECK IF SECOND MEMBER OF PAIR IS ALREADY INCLUDED
   18 IF(H.LT.0.)GO TO 16
      HPOS=H
      HNEG=H-DD
      H=HNEG
      GO TO 17
C             INCLUDE ANOTHER PAIR OF RESONANCES IF REQUIRED
   16 IF(MP.EQ.NP)GO TO 93
      MP=MP+1
C
         CALL WIGNER(DE(J,L1,I),DD)
C
      HPOS=HPOS+DD
      H=HPOS
      GO TO 17
   93 CONTINUE
   83 CONTINUE
    3 CONTINUE
      SN=ST-SN
      T=EXP(-XN(N)*ST)
      TS=T*ST
      IF(IHIST.NE.1)GO TO 52
      IZ=20.*ST/STT(K)+1.
      IF(IZ.LT.1)  IZ=1
      IF(IZ.GT.101)IZ=101
      ZST(IZ)=ZST(IZ)+1.
      IZ=20.*SC/SG(K)+1.
      IF(IZ.LT.1)  IZ=1
      IF(IZ.GT.101)IZ=101
      ZSG(IZ)=ZSG(IZ)+1.
   52 IZ=2.*T/DT+1.
      IF(IZ.LT.1  )IZ=1
      IF(IZ.GT.101)IZ=101
      ZT(IZ)=ZT(IZ)+1.
      TNUM=TNUM+1.
      TSUM=TSUM+T
      TSSM=TSSM+TS
C
         CALL AVERT(TMCG(MM),T,SGMC(MM),SC,STMC(MM),ST,MM,1)
C
    1 CONTINUE
C             END OF HISTORY LOOP
      AVT=TSUM/TNUM
      AVTS=TSSM/TNUM
C             MONTE CARLO RESULTS
C
         CALL AVERT(TMCG(MM),T,SGMC(MM),SC,STMC(MM),ST,MM,2)
C
      TMC=0.
      SGM=0.
      STM=0.
      DO 22 M=1,MM
      SGM=SGM+SGMC(M)
      STM=STM+STMC(M)
   22 TMC=TMC+TMCG(M )
      TMC=TMC/FLOAT(MM)
      SGM=SGM/FLOAT(MM)
      STM=STM/FLOAT(MM)
C             STATISTICAL UNCERTAINTY
      DTMC=0.
      DTS=0.
      DO 23 M=1,MM
      SIM=SIM+(SGMC(M)-SGM)**2
      SYM=SYM+(STMC(M)-STM)**2
      DTS=DTS+(STMC(M)-STM)*(TMCG(M)-TMC)
   23 DTMC=DTMC+(TMCG(M)-TMC)**2
      DEN=FLOAT(MM*(MM-1))
      DSGM=SQRT(SIM/DEN)
      DSTM=SQRT(SYM/DEN)
      DTMC=SQRT(DTMC/DEN)
      DTS2=DTS/DEN
      TAU=EXP(XN(N)*STM)*TMC
      DTAU=(XN(N)*DSTM)**2+2.*XN(N)*DTS2/TMC+(DTMC/TMC)**2
      IF(DTAU.LT.0.0)DTAU=0.0
      DTAU=100.*SQRT(DTAU)
C             WRITE HISTOGRAM VALUES IF REQUIRED
      IF(IHIST.NE.1)GO TO 53
      SUMST=0.
      SUMSG=0.
      SUMT =0.
      DO 54 IZ=1,101
      SUMST=SUMST+ZST(IZ)
      SUMSG=SUMSG+ZSG(IZ)
      SUMT =SUMT + ZT(IZ)
   54 CONTINUE
      DO 55 IZ=1,101
      ZST(IZ)=ZST(IZ)/SUMST
      ZSG(IZ)=ZSG(IZ)/SUMSG
      ZT(IZ)=ZT(IZ)/SUMT
   55 CONTINUE
      WRITE(8,199)(IZ,ZST(IZ),ZSG(IZ),ZT(IZ),IZ=1,101)
  199 FORMAT(1H1,///,' NO.    TALLY FOR HISTOGRAM OF      '           /
     1               '          ST     SG     T           '           //
     2(I4,0PF8.4,2F7.4))
C
CC     CALL PLOT(TN,ZT,100,3,3,1,1,1,1,100.,0.,0.125,100.,0.,0.125,'DISTR
CC   1IBUTION OF TRANSMISSION VALUES',1)
   53 CONTINUE
      RETURN
      END
C
      SUBROUTINE PEPS(RK,L,PL,VL)
C             PEPS CALCULATES CENTRIFUGAL-BARRIER PENETRABILITIES, VL,
C             AND HARD-SPHERE PHASE SHIFTS, PL.
      GO TO (1,2,3,4),L
    1 PL=RK
      VL=1.
      RETURN
    2 PL=RK-ATAN(RK)
      VL=RK**2/(RK**2+1.)
      RETURN
    3 PL=RK-ATAN(3.*RK/(3.-RK**2))
      VL=RK**4/(RK**4+3.*RK**2+9.)
      RETURN
    4 PL=RK-ATAN(RK*(15.-RK**2)/3./(5.-2.*RK**2))
      VL=RK**6/(RK**6+6.*RK**4+45.*RK**2+225.)
      RETURN
      END
C
      SUBROUTINE ENDEP(E,A,B,P,AP,UM,EDGG,EDD)
C             ENDEP CALCULATES THE FACTORS WHICH DESCRIBE THE ENERGY
C             DEPENDENCE OF LEVEL SPACINGS AND RADIATION WIDTHS
      DIMENSION Y1(21),Y2(21)
      E=0.001*E
      A1=A+1.
C             LEVEL DENSITY FACTOR
      U =B+E-P
      UB=B-P
      EDD=RHO(UB,A,P,AP,UM)/RHO(U,A,P,AP,UM)
C             RADIATION WIDTH FACTOR
      ER=35./A1**(1./6.)
      GR=33./A1**(1./3.)
      U1X=B+E
      U2X=B
      DO 1 N=1,21
      U1=U1X*FLOAT(N-1)/20.
      U2=U2X*FLOAT(N-1)/20.
      E1=U1X-U1
      E2=U2X-U2
      F1=(E1**2-ER**2)**2+(GR*ER)**2
      F2=(E2**2-ER**2)**2+(GR*ER)**2
      Y1(N)=E1**4*RHO(U1,A,P,AP,UM)/F1
      Y2(N)=E2**4*RHO(U2,A,P,AP,UM)/F2
    1 CONTINUE
C
         CALL SIMP(Y1,1,21,SUM1)
         CALL SIMP(Y2,1,21,SUM2)
C
      EDGG=EDD*SUM1/SUM2
      E=1000.*E
      RETURN
      END
C
      FUNCTION RHO(U,A,PE,AP,UM)
C
C         GILBERT-CAMERON COMPOSITE LEVEL DENSITY
C         (ALL SPINS AND PARITIES)
C         U:   EXCITATION ENERGY (MEV)
C         A:   ATOMIC WEIGHT OF TARGET NUCLEUS (AMU)
C         PE:  PAIRING CORRECTION (MEV)
C         AP:  A-PARAMETER (1/MEV)
C         UM:  MATCHING ENERGY (MEV)
C
      AI=FLOAT(INT(A+1.5))
      IF(U.GT.UM)GO TO 1
C         LOW ENERGIES: CONSTANT-TEMPERATURE FORMULA
      UEFF=UM-PE
      AU4=AP*UEFF*4.
      VARJ2=0.1460*SQRT(AU4)*CBRT(AI*AI)
      T=1./(SQRT(AP/UEFF)-1.5/UEFF)
      RHO=(AP/3.)*EXP(SQRT(AU4))/AU4**1.25*SQRT(2./VARJ2)
      RHO=RHO*EXP((U-UM)/T)
      RETURN
C         HIGH ENERGIES: FERMI-GAS FORMULA
    1 AU4=AP*(U-PE)*4.
      VARJ2=0.1460*SQRT(AU4)*CBRT(AI*AI)
      RHO=(AP/3.)*EXP(SQRT(AU4))/AU4**1.25*SQRT(2./VARJ2)
      RETURN
      END
C
      SUBROUTINE AVER(SGMC,SC,STMC,ST,NS,I)
      DATA N,J/100,0/
C        N   GROUP SIZE
C        I=1 SUPPLY DATA TO BE GROUPED
C        I=2 FINISH LAST GROUP,WHICH MAY CONTAIN LESS THAN N MEMBERS
      W=1./FLOAT(N)
      GO TO (1,3),I
    1 IF(J.GT.0)GO TO 2
      SGMC=0.
      STMC=0.
      J=1
    2 SGMC=SGMC+SC
      STMC=STMC+ST
      J=J+1
      IF(J.LE.N)RETURN
      J=0
      SGMC=SGMC*W
      STMC=STMC*W
      NS=NS+1
      RETURN
    3 IF(J.EQ.0)GO TO 4
      CORR=1.0/FLOAT(J-1)
      SGMC=SGMC*CORR
      STMC=STMC*CORR
      J=0
      RETURN
    4 NS=NS-1
    5 RETURN
      END
C
      SUBROUTINE AVERP(PZ,PZG,PS,PSG,MM,I)
      DATA N,J/100,0/
C        N   GROUP SIZE
C        I=1 SUPPLY DATA TO BE GROUPED
C        I=2 FINISH LAST GROUP,WHICH MAY CONTAIN LESS THAN N MEMBERS
      W=1./FLOAT(N)
      GO TO (1,3),I
    1 IF(J.GT.0)GO TO 2
      PZ=0.
      PS=0.
      J=1
    2 PZ=PZ+PZG
      PS=PS+PSG
      J=J+1
      IF(J.LE.N)RETURN
      J=0
      PZ=PZ*W
      PS=PS*W
      MM=MM+1
      RETURN
    3 IF(J.EQ.0)GO TO 4
      CORR=1.0/FLOAT(J-1)
      PZ=PZ*CORR
      PS=PS*CORR
      J=0
      RETURN
    4 MM=MM-1
    5 RETURN
      END
C
      SUBROUTINE AVERT(TM,T,SGM,SC,STM,ST,MM,I)
      DATA N,J/100,0/
C        N   GROUP SIZE
C        I=1 SUPPLY DATA TO BE GROUPED
C        I=2 FINISH LAST GROUP,WHICH MAY CONTAIN LESS THAN N MEMBERS
      W=1./FLOAT(N)
      GO TO (1,3),I
    1 IF(J.GT.0)GO TO 2
      TM=0.
      SGM=0.
      STM=0.
      J=1
    2 TM=TM+T
      SGM=SGM+SC
      STM=STM+ST
      J=J+1
      IF(J.LE.N)RETURN
      J=0
      TM=TM*W
      SGM=SGM*W
      STM=STM*W
      MM=MM+1
      RETURN
    3 IF(J.EQ.0)GO TO 4
      CORR=1.0/FLOAT(J-1)
      TM=TM*CORR
      SGM=SGM*CORR
      STM=STM*CORR
      J=0
      RETURN
    4 MM=MM-1
    5 RETURN
      END
C
      SUBROUTINE PFCN(X,Y,U,V,K)
      DATA PI,PI2,N,IS/3.14159 ,6.28318 ,5,0/
      DIMENSION EMN(10)
      IF(Y.GT.0.01)GO TO 40
      IF(IS.NE.0)GO TO 5
      IS=1
      DO 2 I=2,N
    2 EMN(I)=EXP(-FLOAT(I-1)**2)/PI
    5 CONTINUE
      XY=X/Y
      IF(1.+X*X.GT.100.*Y*Y)GO TO 30
      IF(Y.LT.0.001)GO TO 40
      Y2=Y**2
      D=X**2+Y2
      SU=Y/(D*PI)
      SV=X/(D*PI)
      DO 10 I=2,N
      AM=I-1
      XNP=X-AM
      XNM=X+AM
      DP=XNP**2+Y2
      DM=XNM**2+Y2
      SU=SU+EMN(I)*Y*(1./DP+1./DM)
      SV=SV+EMN(I)*(XNP/DP+XNM/DM)
   10 CONTINUE
      IF(Y.GT.PI)GO TO 20
      P=2.0
      IF(Y.EQ.PI)P=1.0
      XY2=X*Y*2.
      SXY2=SIN(XY2)
      CXY2=COS(XY2)
      XP2=X*PI2
      SXP2=SIN(XP2)
      CXP2=COS(XP2)
      EYP2=EXP(Y*PI2)
      EYX=EXP(Y2-X**2)
      D=1.-2.*EYP2*CXP2+EYP2**2
      SU=SU+P*EYX*(CXY2-EYP2*(CXP2*CXY2+SXP2*SXY2))/D
      SV=SV-P*EYX*(SXY2+EYP2*(SXP2*CXY2-CXP2*SXY2))/D
   20 U=SU
      V=SV
      RETURN
   30 F=1.77245  *Y
      U=1./((1.+XY**2)*F)
      V=XY*U
      RETURN
C
   40    CALL PFCNP(X,Y,U,V,K)
C
      RETURN
      END
C
      SUBROUTINE PFCNP(X,Y,U,V,L)
C
C     PFCN YIELDS REAL AND IMAGINARY PART OF THE COMPLEX
C     PROBABILITY INTEGRAL
C
      DIMENSION W287(4),W283(4)
      DATA W283/1.65068,.524648,-.524648,-1.65068/
      DATA W287/.025883,.256212,.256212,.025883/
      II=1
      ASSIGN 244 TO J
      C5=X
      C6=Y
      IF(C5.LT.0.0)GO TO 8
      IF(C6.LT.0.0)GO TO 287
      GO TO 11
    8 IF(C6.GE.0.0)GO TO 14
      ASSIGN 245 TO I
      GO TO 20
   11 ASSIGN 257 TO I
      GO TO 46
   14 ASSIGN 255 TO I
      GO TO 46
   20 Z=C6*C6-C5*C5
      CO=EXP(Z)
      C7=CO+CO
      CO=C5*C6
      C9=CO+CO
      C8=-C7*SIN(C9)
      C7=C7*COS(C9)
   46 C5=ABS(C5)
      C6=ABS(C6)
      IF(C5.GE.6.0)GO TO 219
   50 IF(C6.LE.0.5)GO TO 65
      IF(C6.LE.3.0)GO TO 61
      IF(C6.GT.6.0)GO TO 219
      C9=0.5
      GO TO 73
   61 IF(C6.LE.1.5)GO TO 71
      C9=0.25
      GO TO 73
   65 C10=C6
      C6=0.5
      ASSIGN 128 TO J
   71 C9=0.09375
   73 C11=0.0
      C17=0.0
      C18=0.0
      ASSIGN 123 TO K
   79 C21=C5-C11
      C19=C21*C21
      C20=C6*C6+C19
      T=C11*C11
      C19=EXP(-T)/C20*0.31830*C9
      C17=C19*C6+C17
      C18=C21*C19+C18
  107 GO TO K,(108,123)
  108 II=3-II
      IF(II.EQ.1)GO TO 114
      C11=-C11
      GO TO 79
  114 IF(-C11-4.0.GT.0.0)GO TO J,(128,244)
      C11=-C11+C9
      GO TO 79
  123 II=1
      ASSIGN 108 TO K
      C11=C9
      GO TO 79
  128 C11=C17
      C12=C18
      C9=2.0
      C6=C10-0.5
      C6=C6+C6
      C10=C11/2.0
      C13=(C5*C12+C10-0.56419)*C6
      C10=C12/2.0
      C14=(-C5*C11+C10)*C6
      C17=C11+C13
      C18=C12+C14
  165 C10=C6/C9
      C19=C13/2.0
      C19=C5*C14+C19
      C15=(C6/2.0*C11+C19)*C10
      C17=C15+C17
      T1=C5*C13
      C19=(C6*C12+C14)/2.0
      C16=(-T1+C19)*C10
      C18=C16+C18
      T1=C17+C15
      IF((T1-C17).NE.0.0)GO TO 207
      T1=C18+C16
      IF((T1-C18).EQ.0.0)GO TO 244
  207 C11=C13
      C12=C14
      C13=C15
      C14=C16
      C9=C9+1.0
      GO TO 165
  219 C17=0.0
      C18=0.0
      DO 230 M=1,4
      C12=C5-W283(M)
      C11=C12*C12
      C11=C6*C6+C11
      C11=W287(M)/C11
      C17=C11*C6+C17
      C18=C11*C12+C18
  230 CONTINUE
  244 GO TO I,(245,249,255,257)
  245 C8=-C8
      C18=-C18
  249 C17=C7-C17
      C18=C8-C18
  255 C18=-C18
  257 U=C17
      V=C18
      L=0
      RETURN
  287 C5=-C5
      ASSIGN 249 TO I
      GO TO 20
      END
C
      SUBROUTINE WIGNER(DAV,D)
C             WIGNER YIELDS LEVEL SPACINGS WITH THE WIGNER DISTRIBUTION.
C             DAV IS THE AVERAGE LEVEL SPACING.
    1 WWW=RANDOM(A)
      IF(WWW.EQ.0.) GO TO 1
      D=1.12837  *DAV*SQRT(-ALOG(WWW))
      RETURN
      END
C
      SUBROUTINE PORTER(GNAV,GN)
C             PORTER YIELDS NEUTRON WIDTHS WITH THE PORTER-THOMAS
C             DISTRIBUTION. GNAV IS THE AVERAGE WIDTH.
C             V. NEUMANN SCHEME
    1 R=RANDOM(D)*0.929
      RR=1./(1.-R)**2
      X=R**2*RR
      IF(X.GT.170.) GO TO 1
      IF(RANDOM(D).GT.0.560*EXP(-X)*RR)GO TO 1
      GN=2.*GNAV*X
      RETURN
      END
C
      SUBROUTINE SPACE(DAV,D)
C             SPACE YIELDS LEVEL SPACINGS WITH A D-WEIGHTED WIGNER
C             DISTRIBUTION. DAV IS THE AVERAGE OF THE WIGNER DISTRIBU-
C             TION. V. NEUMANN SCHEME
    1 R=RANDOM(A)
      X=(R/(1.-R))**0.3333
      IF(X.GT.13.) GO TO 1
      IF(RANDOM(A).GT.0.494*EXP(-X**2)/(1.-R)**2)GO TO 1
      D=1.12837  *DAV*X
      RETURN
      END
      FUNCTION RANDOM(DUMMY)
      COMMON/RANDM/ IY
      DATA EPS/1.E-15/
      IX=IY
C
         CALL RANDU(IX,IY,RA)
C
      IF(RA.LT.EPS)RA=EPS
      RANDOM=RA
      RETURN
      END
C
      SUBROUTINE SIMP(Y,M,N,SUM)
C
C         SIMPSON'S RULE
C
      DIMENSION Y(21)
      SUM=0.
      K1=M+1
      K2=N-1
      DO 1 K=K1,K2,2
    1 SUM=SUM+Y(K-1)+4.*Y(K)+Y(K+1)
      SUM=SUM/3.
      RETURN
      END
C
      FUNCTION CBRT(MANT)
C
      EXPO=1./3.
      CBRT=MANT**EXPO
      RETURN
      END
C
C
C
C        SUBROUTINE RANDU
C
C        PURPOSE
C           COMPUTES UNIFORMLY DISTRIBUTED RANDOM REAL NUMBERS BETWEEN
C           0 AND 1.0 AND RANDOM INTEGERS BETWEEN ZERO AND
C           2**31. EACH ENTRY USES AS INPUT AN INTEGER RANDOM NUMBER
C           AND PRODUCES A NEW INTEGER AND REAL RANDOM NUMBER.
C
C        USAGE
C           CALL RANDU(IX,IY,YFL)
C
C        DESCRIPTION OF PARAMETERS
C           IX - FOR THE FIRST ENTRY THIS MUST CONTAIN ANY ODD INTEGER
C                NUMBER WITH NINE OR LESS DIGITS. AFTER THE FIRST ENTRY,
C                IX SHOULD BE THE PREVIOUS VALUE OF IY COMPUTED BY THIS
C                SUBROUTINE.
C           IY - A RESULTANT INTEGER RANDOM NUMBER REQUIRED FOR THE NEXT
C                ENTRY TO THIS SUBROUTINE. THE RANGE OF THIS NUMBER IS
C                BETWEEN ZERO AND 2**31
C           YFL- THE RESULTANT UNIFORMLY DISTRIBUTED, FLOATING POINT,
C                RANDOM NUMBER IN THE RANGE 0 TO 1.0
C
C        REMARKS
C           THIS SUBROUTINE IS SPECIFIC TO SYSTEM/360 AND WILL PRODUCE
C           2**29 TERMS BEFORE REPEATING.  THE REFERENCE BELOW DISCUSSES
C           SEEDS (65539 HERE), RUN PROBLEMS, AND PROBLEMS CONCERNING
C           RANDOM DIGITS USING THIS GENERATION SCHEME.  MACLAREN AND
C           MARSAGLIA, JACM 12, P. 83-89, DISCUSS CONGRUENTIAL
C           GENERATION METHODS AND TESTS.  THE USE OF TWO GENERATORS OF
C           THE RANDU TYPE, ONE FILLING A TABLE AND ONE PICKING FROM THE
C           TABLE, IS OF BENEFIT IN SOME CASES.  65549 HAS BEEN
C           SUGGESTED AS A SEED WHICH HAS BETTER STATISTICAL PROPERTIES
C           FOR HIGH ORDER BITS OF THE GENERATED DEVIATE.
C           SEEDS SHOULD BE CHOSEN IN ACCORDANCE WITH THE DISCUSSION
C           GIVEN IN THE REFERENCE BELOW.  ALSO, IT SHOULD BE NOTED THAT
C           IF FLOATING POINT RANDOM NUMBERS ARE DESIRED,AS ARE
C           AVAILABLE FROM RANDU, THE RANDOM CHARACTERISTICS OF THE
C           FLOATING POINT DEVIATES ARE MODIFIED AND IN FACT THESE
C           DEVIATES HAVE HIGH PROBABILITY OF HAVING A TRAILING LOW
C           ORDER ZERO BIT IN THEIR FRACTIONAL PART.
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           POWER RESIDUE METHOD DISCUSSED IN IBM MANUAL C20-8011,
C           RANDOM NUMBER GENERATION AND TESTING
C
C     ..................................................................
C
      SUBROUTINE RANDU(IX,IY,YFL)
      IY=IX*65539
      IF(IY)5,6,6
    5 IY=IY+2147483647+1
    6 YFL=IY
      YFL=YFL*.4656613E-9
      RETURN
      END


