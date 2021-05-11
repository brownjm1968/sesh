program sesh
!C ---------------------------------------------------------------------    
!C Modified for AIX xlf Fortran compiler  ( Klaus Berthold, IRMM 4/94)
!C ---------------------------------------------------------------------    
!C
!C(Revised and extended by Konstantin Volev, IRMM 5/2002 )
!C
!C MAIN :  107 FORMAT Statements  (no more using Hollerith count)
!C         109 
!C         113  
!C         114
!C         115
!C         116
!C
!C         COMMON BLOCK : CHARACTER*26 CDATE
!C 
!C SUBROUTINE MUSC :  3 CONTINUE (one removed)
!C                    3 CONTINUE 
!C 
!C SUBROUTINE MOCT :    PLOT(....) commented out
!C  
!C FUNCTION RHO    :    CBRT not intrinsic so added FUNCTION CBRT
!C
!C 
!C ---------------------------------------------------------------------    
!     SESH MAIN PROGRAM (KFK VERSION 1975, EXPERIMENTAL COPY 1994)
!
!     MANUAL: F.H. FROEHNER,
!             "SESH - A FORTRAN IV CODE FOR CALCULATING THE SELF-
!             SHIELDING AND MULTIPLE SCATTERING EFFECTS FOR
!             NEUTRON CROSS SECTION DATA INTERPRETATION
!             IN THE UNRESOLVED RESONANCE REGION",
!             REPORT GA-8380 (1968)
!
!     GILBERT-CAMERON COMPOSITE LEVEL DENSITY AND GIANT DIPOLE RESONANCE
!     MODEL ARE USED FOR THE ENERGY DEPENDENCE OF LEVEL SPACINGS AND
!     RADIATION WIDTHS IN VERSION 1975. THE INPUT REMAINS AS DESCRIBED
!     IN GA-8380 WITH TWO EXCEPTIONS: (1) THE NUCLEAR TEMPERATURE IS
!     REPLACED BY THE GILBERT-CAMERON PAIRING ENERGY OF THE COMPOUND
!     NUCLEUS, (2) ONLY THE S-WAVE LEVEL SPACING MUST BE GIVEN (SPACING
!     INPUT FOR L>0 IS IGNORED). FURTHERMORE, CHANNEL RADIUS AND
!     EFFECTIVE NUCLEAR RADIUS ARE DISTINCT FOR ALL PARTIAL WAVES:
!     THE CHANNEL RADIUS IS TAKEN AS RC(I)=(1.23*A**(1/3)+0.80) FM,
!     THE EFFECTIVE NUCLEAR RADII R(L,I) ARE INPUT NUMBERS WITH THE
!     DEFAULT VALUES R(L,I)=RC(I).
!
!     INPUT LIMITATIONS:
!             UP TO      10 ISOTOPES
!             "  "        4 PARTIAL WAVES
!             "  "        5 SAMPLE GEOMETRIES
!             "  "       24 NEUTRON ENERGIES
!             "  "   100000 MONTE CARLO HISTORIES PER ENERGY
    use mc_routines

    implicit none
!
    real(4), allocatable, dimension(:) :: COMM,ZL,AP,UM,A,SC,ST,AB,BE,PE,TEFF, &
                                          SPIN,AA,RC,XN,E,SG,SP,RN,RX,ZH,RB
    integer, allocatable, dimension(:) :: JX2,JN2,JMX,NL
    real(4), allocatable, dimension(:,:) :: SGI,GG,D,S,SI,R,EFF,SUMGJ,STI,SPI
    integer, allocatable, dimension(:,:) :: J2X,J2N
    real(4), allocatable, dimension(:,:,:) :: SGL,G,GNR,GNRIN,DJL,STL,SPL
    integer, allocatable, dimension(:,:,:) :: LJX,LJN

    real(4) :: AI,AK,CC,DPO,DPS,DSGM,DSSCF,DSTM,DTAU,DTMC,DUMMY,EDD,EDGG,EJL, &
               ETA,FJ,GN,GNIN,PL,PO,PS,PSI0,QI,RK,S3,SGM,SNC,SQ,SSCF,STM,SUM, &
               TAU,TMC,U,V,VARJ2,VL,X2I,X2J,XNSS,XO,XX,ZP

    integer :: I,I2,IHIST,ITYPE,ITYPO,IY,J,J2,J2MN,J2MX,JX,K,KQ,KZ,L,LL,LX,M, &
               M4,MJ,N,NE,NI,NN,NQ

!    CHARACTER*26 CDATE
! FG051198
    CHARACTER*30 CDATE
! modification for HP on bonsai3
! ---- Jesse Mod --------------------------------------------------------------
    CHARACTER (LEN=100) :: inp_file_name
    CHARACTER (LEN=100) :: out_file_name
    CHARACTER (LEN=100) :: results_file_name
    CHARACTER (LEN=100) :: cor_file_name

    ! --- dimension variables now -------
    allocate(COMM(18),ZL(4),SGI(10,100),SGL(10,100,4),AP(10),UM(10),A(11),      &
             SC(100),ST(100),AB(11),BE(11),PE(11),TEFF(11),SPIN(11),NL(10),     &
             AA(10),RC(10),GG(5,11),D(5,11),S(5,11),SI(5,11),R(5,11),EFF(5,11), &
             J2X(4,10),J2N(4,10),SUMGJ(4,10),XN(6),E(100),SG(100),SP(100),      &
             G(8,4,10),GNR(8,4,10),GNRIN(8,4,10),DJL(8,4,10),RN(6),RX(6),       &
             ZH(100),LJX(2,8,10),LJN(2,8,10),JX2(10),JN2(10),JMX(10),RB(6),     &
             STI(10,100),SPI(10,100),STL(10,100,4),SPL(10,100,4))

    ! -----------------------------------

    print *, "Input file name?"
    read(*,*) inp_file_name
    print *, "Output file name?"
    read(*,*) out_file_name
    print *, "Analytical results file name?"
    read(*,*) results_file_name
    print *, "Correction file name?"
    read(*,*) cor_file_name
! -----------------------------------------------------------------------------
    open (unit=5,file=inp_file_name,status='old')           ! -- Input file
    open (unit=8,file=out_file_name,status='unknown')       ! -- Main output file
    open (unit=11,file=results_file_name,status='unknown')  ! -- Formatted partial wave cross section output (analytical calc.))
    open (unit=12, file=cor_file_name, status='unknown')    ! -- Formatted self-shielding correction factor output
    call random_seed
! end modification for HP on bonsai3
!
!             READ INPUT
!             IHIST IS A TAG FOR HISTOGRAMS OF ST, SG, T, F0;
!             IHIST=1  MEANS HISTOGRAM VALUES ARE CALCULATED AND PRINTED
!             IHIST=0  MEANS NO SUCH VALUES ARE CALCULATED OR PRINTED
!     CALL ERRSET(217,0  ,-1,1,0,100)
!     CALL ERRSET(208,300,-1,0,1,208)
    IY=4751
1   READ(5,100,END=999)COMM,IHIST
100 FORMAT(18A4,1I8)
    SUM=0.
    I=0
    2 I=I+1
!             I-TH ISOTOPE (MAXIMUM OF 10 ISOTOPES):
!             NUCLEON NUMBER A(I), ABUNDANCE AB(I), BINDING ENERGY BE(I) [MeV], PAIRING
!             ENERGY PE(I) [MeV], EFFECTIVE SAMPLE TEMPERATURE  TEFF(I) [DEG. K],
!             TARGET SPIN QUANTUM NUMBER SPIN(I);
    READ(5,101)A(I),AB(I),BE(I),PE(I),TEFF(I),SPIN(I)
101 FORMAT(6E10.5)
    L=0
3   L=L+1
!             L-TH PARTIAL WAVE (UP TO F-WAVE):
!             RADIATION WIDTH GG(L,I) [eV], AVERAGE LEVEL SPACING D(L,I) [eV], STRENGTH
!             FUNCTION FOR ELASTIC SCATTERING S(L,I), STRENGTH FUNCTION FOR
!             INELASTIC SCATTERING SI(L,I), NUCLEAR RADIUS R(L,I) [fm], DETECTION
!             EFFICIENCY
    READ(5,102)GG(L,I),D(L,I),S(L,I),SI (L,I),R(L,I),EFF(L,I)
102 FORMAT(6E10.5)
!             CHECK FOR SAMPLE THICKNESS CARD  (FIRST WORD ZERO)
    IF(GG(L,I).EQ.0.0)GO TO 4
!             CHECK FOR LAST PARTIAL WAVE
    IF(S(L,I).LT.1.)GO TO 3
!             LAST CARD WAS ISOTOPE CARD. STORE CORRESPONDINGLY
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
!             LAST CARD WAS SAMPLE THICKNESS CARD (NUCLEI/B),
!             STORE CORRESPONDINGLY.
4   NI=I
    NL(I)=L-1
    XN(1)=D(L,I)
    XN(2)=S(L,I)
    XN(3)=SI(L,I)
    XN(4)=R(L,I)
    XN(5)=EFF(L,I)
    XN(6)=0.

!             FIND NUMBER OF SAMPLE THICKNESSES
!             FOR TRANSMISSION AND CAPTURE DATA

    DO 5 N=2,6
    IF(XN(N).NE.0.)GO TO 5
    NN=N-1
    GO TO 6
5   CONTINUE

!             READ OUTER RADII (NUCLEI/BARN)

6   READ(5,103)(RX(N),N=1,5)
    RX(6)=0.
    IF(NN.GT.1)GO TO 8

!             FIND NUMBER OF SHELL THICKNESSES FOR SHELL TRANSMISSION
!             DATA

    DO 7 N=2,6
    IF(RX(N).NE.0.)GO TO 7
    NN=N-1
    GO TO 8
7   CONTINUE


!             READ INNER RADII (NUCLEI/BARN) FOR SHELL TRANSMISSION
!             CALCULATIONS. RN(N)=0. MEANS CYLINDRICAL SAMPLE


8   READ(5,103)(RN(N),N=1,5)

103 FORMAT(E20.5,4E10.5)
      
!         BJM Modification 10/20/2015-- READ BEAM RADIUS (NUCLEI/BARN)
!                         FOR CAPTURE SAMPLES RB(N)=0.0 MEANS 
!                         BEAM RADIUS >= SAMPLE RADIUS      
!      
88  READ(5,103)(RB(N),N=1,5)

!             READ NUMBER OF RESONANCE PAIRS


    READ(5,104)ZP
104 FORMAT(E10.5)
    DO 9 M=1,25
    M4=(M-1)*4

!             READ ENERGIES (KEV) AND NUMBERS OF MONTE CARLO HISTORIES

    READ(5,105)E(M4+1),ZH(M4+1),E(M4+2),ZH(M4+2) &
                ,E(M4+3),ZH(M4+3),E(M4+4),ZH(M4+4)
105 FORMAT(8E10.5)


!             BLANK CARD SIGNALS END OF INPUT

    IF(E(M4+1).EQ.0.)GO TO 10
9   CONTINUE


!             FIND NUMBER OF ENERGIES, NE

10  DO 11 J=1,4
        MJ=M4-4+J
        IF(E(MJ).GT.0.)NE=MJ
11  CONTINUE

!             WRITE INPUT HEADING
!C    VMS SPECIFIC DATE AND TIME PROCEDURES
!      CALL FDATE_(CDATE)
! FG051198
! modification for HP on bonsai3
!     CALL FDATE(CDATE)
!     WRITE(8,106)CDATE,COMM
! 106 FORMAT(1H1,//,'[UNIX-SESH ',A24,']     ',18A4)
! end modification for HP on bonsai3
    WRITE(8,107)
107 FORMAT( &
    ' I N P U T '/ &
    ' ========= '// &
    ' NUCLEON   ABUN-    BINDING  PAIRING EFF.      NUCL.  ', &
    'ORB.ANG.  AV. RAD.    AV. LEVEL   STRENGTH    STRENGTH', &
    '    NUCLEAR     EFFI- '/ &
    ' NUMBER    DANCE    ENERGY   ENERGY  TEMP.     SPIN   ', &
    ' MOM.     WIDTH       SPACING     FCT. FOR    FCT. FOR', &
    '    RADIUS      CIENCY'/ &
    '                    (MEV)    (MEV)   (DEG.K)   Q.NO.  ', &
    ' Q.NO.    (EV)        (EV)        EL. SCATT.  INEL. SC', &
    '.   (FM)              '/)


!             PREPARE ENERGY-INDEPENDENT PARAMETERS AND PRINT INPUT

    DO 16 I=1,NI

!             CHANNEL RADIUS

        RC(I)=1.23*A(I)**(1./3.)+0.8

!             GILBERT-CAMERON MATCHING ENERGY

        AI=FLOAT(INT(A(I)+1.5))
        write (6,*)'a= ',a(i),ai
        UM(I)=2.5+150./AI

!             FIND FERMI GAS MODEL A-PARAMETER FROM S-WAVE SPACING

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
12      CONTINUE
13      AP(I)=XX**2/(4.*U)

!             DETERMINE MINIMUM AND MAXIMUM COMPOUND SPIN POSSIBLE

        X2I=2.*SPIN(I)
        I2=X2I+.01
        LX=NL(I)
        DO 15 L=1,LX
        J2X(L,I)=I2+2*L-1
        J2N(L,I)=1
        IF(I2.GT.2*L-2)J2N(L,I)=I2-2*L+1
        IF(I2.LT.2*L-2)J2N(L,I)=2*L-I2-3

!             CALCULATE SUM OF SPIN FACTORS

        SUMGJ(L,I)=FLOAT((J2X(L,I)+J2N(L,I)+2)*(J2X(L,I)-J2N(L,I)+2))/8./(X2I+1.)

!             LEVEL DENSITY FOR GIVEN L
        IF(L.GT.1)D(L,I)=D(1,I)/SUMGJ(L,I) 
!       IF(L.GT.1)D(L,I)=D(L,I)-- modified by BJM 2/5/2015 to override auto-calculation of level spacings
!                                       and accept user input
!    -- End BJM	  
        J=0
        J2MN=J2N(L,I)
        J2MX=J2X(L,I)
        DO 14 J2=J2MN,J2MX,2
        J=J+1

!             CALCULATE SPIN FACTOR

        X2J=FLOAT(J2)
        G(J,L,I)=(X2J+1.)/2./(X2I+1.)

!             LEVEL DENSITY FOR GIVEN L AND J

        DJL(J,L,I)=D(1,I)/G(J,L,I) 
!       DJL(J,L,I)=D(L,I)/G(J,L,I)-- modified by BJM 2/5/2015 to override auto-calculation of level spacings
!                                    and accept user input
!       PRINT*, DJL(J,L,I)
!       -- End BJM 
        ZL(L)=FLOAT(L-1)
14      CONTINUE
        IF(R(L,I).EQ.0.)R(L,I)=RC(I)
15      CONTINUE
        WRITE(8,108)A(I),AB(I),BE(I),PE(I),TEFF(I),SPIN(I),(ZL(L),GG(L,I),D(L,I),S(L,I),SI (L,I),R(L,I),EFF(L,I),L=1,LX)
108     FORMAT(F6.1,F11.4,F8.3,F9.3,2F8.1,F7.0,1PE16.4,1P,4E12.4,0PF9.4/(50X,0PF7.0,1PE16.4,1P,4E12.4,0PF9.4))
16  CONTINUE


!             FIRST PART
!             ANALYTICAL CROSS SECTION CALCULATION FOR ALL ENERGIES
!             BEGIN K-LOOP (ENERGIES)
!
    PRINT*, "Running analytical calculations..."

    DO 21 K=1,NE
        SG(K)=0.
        SC(K)=0.
        SP(K)=0.
        ST(K)=0.

!             BEGIN ISOTOPE LOOP
        DO 20 I=1,NI
            SGI(I,K)=0.
            STI(I,K)=0.
            SPI(I,K)=0.
            DO 17 L=1,4
17          SGL(I,K,L)=0.
            AA(I)=(1.+1./A(I))**2

!             GET ENERGY DEPENDENCE FACTORS FOR LEVEL SPACINGS AND
!             RADIATION WIDTHS
!
            CALL ENDEP(E(K),A(I),BE(I),PE(I),AP(I),UM(I),EDGG,EDD)
! -- Modified by BJM on 2/6/2015 to remove energy dependence
!      EDD=1
!      EDGG=1
! -- End BJM	 
!             BEGIN LOOP OF PARTIAL WAVES
            LX=NL(I)
            DO 19 L=1,LX

!                 GET PENETRABILITIES, SHIFTS, HARD SPHERE PHASE SHIFTS
                AK=RC(I)*SQRT(E(K)/AA(I))/143.92
!               AK = radius*sqrt(Energy/A_mass)/H_bar
!
                CALL PEPS(AK,L,DUMMY,VL)
!
                RK=AK*R(L,I)/RC(I)
!
                CALL PEPS(RK,L,PL,DUMMY)
!
!                 CALCULATE SUM OVER ALL POSSIBLE COMPOUND SPINS
                S3=0.
                J=0
                J2MN=J2N(L,I)
                J2MX=J2X(L,I)
                DO 18 J2=J2MN,J2MX,2
                    J=J+1

!                     FIND NUMBER OF POSSIBLE CHANNEL SPINS FOR GIVEN J AND L
                    EJL=1.
                    IF(J2.LT.J2X(L,I).AND.J2.GT.J2N(L,I).OR.2*L.EQ.I2+2.AND.L.NE.1.0.AND.J2.EQ.J2N(L,I)) EJL=2.

!                     REDUCED WIDTHS FOR ELASTIC AND INELASTIC SCATTERING
                    GNR(J,L,I)  =S(L,I) *EJL*DJL(J,L,I)
                    GNRIN(J,L,I)=SI(L,I)*EJL*DJL(J,L,I)

!                     NEUTRON WIDTHS FOR ELASTIC AND INELASTIC SCATTERING
                    SQ=SQRT(E(K)*1000.)*VL
                    GN  =GNR  (J,L,I)*SQ
                    GNIN=GNRIN(J,L,I)*SQ

!                     FIND PSI-FUNCTION GIVING PORTER-THOMAS AVERAGE
                    ETA=SQRT((GG(L,I)*EDGG+GNIN)/(2.*GN*EDD))
!
                    CALL PFCN(0.,ETA,U,V,KZ)
!
                    PSI0=1.77245  *ETA*U

!                     FORM J-SUM (COMPOUND SPINS)
                    S3=S3+G(J,L,I)**2*(1.-PSI0)
18              CONTINUE

!                 CONTRIBUTION TO EFFECTIVE CAPTURE CROSS SECTION
                SGL(I,K,L)=4.09E3*AA(I)*AB(I)*GG(L,I)*EDGG*S3/E(K)/D(L,I)/EDD/SUMGJ(L,I)*EFF(L,I)

!                 PURE CAPTURE CROSS SECTION

                SC(K)=SC(K)+SGL(I,K,L)/EFF(L,I)

!                 CONTRIBUTION TO POTENTIAL SCATTERING CROSS SECTION

                SPL(I,K,L)=2.605E3/E(K)*AA(I)*AB(I)*(2.*FLOAT(L)-1.)*SIN(PL)**2

!                 CONTRIBUTION TO TOTAL CROSS SECTION

                STL(I,K,L)=SPL(I,K,L)+4.09E6/SQRT(1000.*E(K))*AA(I)*AB(I)*S(L,I)*(2.*FLOAT(L)-1.)*VL*COS(2.*PL)

!                 FORM L-SUMS (PARTIAL WAVES)

                SGI(I,K)=SGI(I,K)+SGL(I,K,L)
                SPI(I,K)=SPI(I,K)+SPL(I,K,L)
                STI(I,K)=STI(I,K)+STL(I,K,L)
19          CONTINUE

!             FORM I-SUMS (ISOTOPES)
!             AVERAGE EFFECTIVE CAPTURE CROSS SECTION, K-TH ENERGY

            SG(K)=SG(K)+SGI(I,K)

!             AVERAGE TOTAL CROSS SECTION, K-TH ENERGY

            ST(K)=ST(K)+STI(I,K)

!             POTENTIAL SCATTERING CROSS SECTION, K-TH ENERGY

            SP(K)=SP(K)+SPI(I,K)
20      CONTINUE
21  CONTINUE
    WRITE(8,109)
109 FORMAT(//,//, &
    ' A N A L Y T I C A L   R E S U L T S '/ &
    ' =================================== '// &
    ' AVERAGE CROSS SECTIONS                               ', &
    'ISOTOPIC CONTRIBUTIONS                    PARTIAL-WAVE', &
    ' CONTRIBUTIONS      '// &
    ' NEUTRON   TOTAL    POTENTIAL     PURE    EFFECTIVE   ', &
    'NUCLEON  TOTAL    POTENTIAL  EFFECTIVE    L    TOTAL  ', &
    '  POTENTIAL  EFFECTIVE'/ &
    '  ENERGY           SCATTERING    CAPTURE   CAPTURE    ', &
    'NUMBER           SCATTERING   CAPTURE                 ', &
    ' SCATTERING   CAPTURE'/ &
    '  (KEV)    (BARN)     (BARN)     (BARN)    (BARN)     ', &
    '         (BARN)     (BARN)    (BARN)           (BARN) ', &
    '    (BARN)    (BARN) '/)
! Loop over energies
    DO 25 K=1,NE
        WRITE(8,110)E(K),ST(K),SP(K),SC(K),SG(K),A(1),STI(1,K),SPI(1,K),SGI(1,K),STL(1,K,1),SPL(1,K,1),SGL(1,K,1)
       
110     FORMAT(/,0PF7.1,1P,4E11.3,0PF7.0,1P,3E11.3,4H   0,1P,3E11.3)
! Loop over isotopes
        DO 24 I=1,NI
            WRITE(11, 212) E(K),STL(I,K,1),SPL(I,K,1),SGL(I,K,1), &
            STL(I,K,2),SPL(I,K,2),SGL(I,K,2),STL(I,K,3),SPL(I,K,3), &
            SGL(I,K,3)
            IF(I.EQ.1)GO TO 22
            WRITE(8,111) A(I),STI(I,K),SPI(I,K),SGI(I,K),STL(I,K,1),SPL(I,K,1),SGL(I,K,1)
111         FORMAT(/,051X           ,0PF7.0,1P,3E11.3,4H   0,1P,3E11.3)
22          IF(NL(I).EQ.1)GO TO 24
            LX=NL(I)
!           Loop over L-values
            DO 23 L=2,LX
                LL=L-1
                WRITE(8,112)LL,STL(I,K,L),SPL(I,K,L),SGL(I,K,L)
112             FORMAT(94X,I1,1P,3E11.3)
212             FORMAT(0PF7.1,9E11.3)
23          CONTINUE
24      CONTINUE
25  CONTINUE
!             END OF ANALYTICAL CALCULATION
!
    print *, "Done!"
!
    IF(ZH(1).LT.1.)GO TO 1
!             SECOND PART
!             MONTE CARLO CALCULATION OF CROSS SECTIONS AND OF
!             PROBABILITY FOR DETECTED CAPTURE
!
    print *, "Beginning Monte Carlo simulations..."
!
    IF(XN(1).GT.0.0.AND.RX(1).GT.0.0.AND.RN(1).EQ.0.0)ITYPO=1 ! -- Circular capture sample geometry
    IF(XN(1).EQ.0.0.AND.RX(1).GT.0.0.AND.RN(1).GT.0.0)ITYPO=2 ! -- Spherical shell geometry
    IF(XN(1).GT.0.0.AND.RX(1).EQ.0.0.AND.RN(1).EQ.0.0)ITYPO=3 ! -- Transmission geometry
    IF(XN(1).GT.0.0.AND.RX(1).GT.0.0.AND.RN(1).GT.0.0)ITYPO=4 ! -- Self-indication geometry
    IF(ITYPO.EQ.1)WRITE(8,113)
    IF(ITYPO.EQ.2)WRITE(8,114)
    IF(ITYPO.EQ.3)WRITE(8,115)
    IF(ITYPO.EQ.4)WRITE(8,116)
113 FORMAT(//,//, &
    ' M O N T E   C A R L O   R E S U L T S '/ &
    ' ===================================== '// &
    ' CAPTURE DATA, CIRCULAR DISC SAMPLE    '// &
    ' SAMPLE    NEUTRON  AV. EFF.  AVERAGE   FIRST-COLL.  P', &
    'ROB. FOR   SELF-SHIELD.  AVERAGE     NUMBER OF  NUMBER', &
    ' OF  SAMPLE  AV.   '/ &
    ' THICK-    ENERGY   CAPTURE   TOTAL     CAPTURE      C', &
    'APTURE     CORRECTION    NUMBER OF   HISTORIES  RESONA', &
    'NCE  RADIUS  TRANS-'/ &
    ' NESS               CROSS     CROSS     PROBAB./     A', &
    'FTER 1     FACTOR/       COLLISIONS             PAIRS ', &
    '             MISS./'/ &
    '                    SECTION/  SECTION/  UNCERT.      C', &
    'OLLISION/  UNCERT.                                    ', &
    '             UNCERT.'/ &
    '                    UNCERT.   UNCERT.                U', &
    'NCERT. '/ &
    ' (NUC./B)   (KEV)   (B)/(B)   (B)/(B)                 ', &
    '             /(PERCENT)                               ', &
    '   (NUC./B)' /)
114 FORMAT(//,//, &
    ' M O N T E   C A R L O   R E S U L T S '/ &
    ' ===================================== '// &
    ' CAPTURE DATA, SPHERICAL SHELL'// &
    ' SHELL     NEUTRON  AV. EFF.  AVERAGE   FIRST-COLL.  P', &
    'ROB. FOR   SHELL         AVERAGE     NUMBER OF  NUMBER', &
    ' OF  SAMPLE  AV.    '/ &
    ' THICK-    ENERGY   CAPTURE   TOTAL     CAPTURE      C', &
    'APTURE     TRANS-        NUMBER OF   HISTORIES  RESONA', &
    'NCE  RADIUS  TRANS- '/ &
    ' NESS               CROSS     CROSS     PROBAB./     A', &
    'FTER 1     MISSION/      COLLISIONS             PAIRS ', &
    '             MISS./ '/ &
    '                    SECTION/  SECTION/  UNCERT.      C', &
    'OLLISION/  UNCERT.                                    ', &
    '             UNCERT.'/ &
    '                    UNCERT.   UNCERT.                U', &
    'NCERT. '/ &
    ' (NUC./B)   (KEV)   (B)/(B)   (B)/(B)                 ', &
    '                           (NUC./B)' /)
115 FORMAT(//,//, &
    ' M O N T E   C A R L O   R E S U L T S '/ &
    ' ===================================== '// &
    ' TRANSMISSION DATA '// &
    ' SAMPLE    NEUTRON  AV. EFF.  AVERAGE   AVERAGE   SELF', &
    '-SHIELD.  NUMBER OF  NUMBER OF '/ &
    ' THICK-    ENERGY   CAPTURE   TOTAL     TRANS.    CORR', &
    'ECTION    HISTORIES  RESONANCE '/ &
    ' NESS               CROSS     CROSS     MISSION/  FACT', &
    'OR/                  PAIRS     '/ &
    '                    SECTION/  SECTION/  UNCERT.   UNCE', &
    'RT. '/ &
    '                    UNCERT.   UNCERT.                 ', &
    '    '/ &
    ' (NUC./B)   (KEV)   (B)/(B)   (B)/(B)               /(', &
    'PERCENT) '/)
116 FORMAT(//,//, &
    ' M O N T E   C A R L O   R E S U L T S '/ &
    ' ===================================== '// &
    ' SELF-INDICATION DATA, DISC SAMPLES '// &
    ' CAPTURE   NEUTRON  AV. EFF.  AVERAGE   FIRST-COLL.  P', &
    'ROB. FOR   SELF-SHIELD.  AVERAGE     NUMBER OF  NUMBER', &
    ' OF  SAMPLE  FIRST '/ &
    ' SAMPLE    ENERGY   CAPTURE   TOTAL     CAPTURE      C', &
    'APTURE     CORRECTION    NUMBER OF   HISTORIES  RESONA', &
    'NCE  RADIUS  SAMPLE'/ &
    ' THICK.             CROSS     CROSS     PROBAB./     A', &
    'FTER 1     FACTOR/       COLLISIONS             PAIRS ', &
    '             THICK.'/ &
    '                    SECTION/  SECTION/  UNCERT.      C', &
    'OLLISION/  UNCERT. '/ &
    '                    UNCERT.   UNCERT.                U', &
    'NCERT. '/ &
    ' (NUC./B)   (KEV)   (B)/(B)   (B)/(B)                 ', &
    '           /(PERCENT)                                 ', &
    '    (NUC./B) (NUC./B)'/)
!             DETERMINATION OF EXTREME VALUES OF L FOR EACH PARITY AND J
    DO 36 I=1,NI
        LX=NL(I)
!             TWICE EXTREME COMPOUND SPIN QUANTUM NUMBERS
        I2=2.*SPIN(I)
        JX2(I)=I2+2*LX-1
        JN2(I)=I2-2*LX+1
        IF(I2.GT.6)GO TO 26
        IF(MOD(I2,2).EQ.0)JN2(I)=1
        IF(MOD(I2,2).EQ.1)JN2(I)=0
!             JMX(I): NUMBER OF POSSIBLE J VALUES
26      JMX(I)=(JX2(I)-JN2(I))/2+1
        JX=JMX(I)
!             J: COUNTER FOR COMPOUND SPIN IN ASCENDING ORDER
        DO 35 J=1,JX
            J2=JN2(I)+2*J-2
!             M: PARITY LABEL, M=1 FOR PARITY OF TARGET GROUND STATE
            M=1
!             EXTREMAL L VALUES FOR GIVEN PARITY, J, AND ISOTOPE
27          IF(LX.GT.M+1)GO TO 30
            IF(J2N(M,I).LE.J2.AND.J2X(M,I).GE.J2)GO TO 29
28          LJX(M,J,I)=-1
            LJN(M,J,I)=-1
            GO TO 34
29          LJX(M,J,I)=M
            LJN(M,J,I)=M
            GO TO 34
30          IF(J2N(M  ,I).LE.J2.AND.J2X(M  ,I).GE.J2)GO TO 31
            IF(J2N(M+2,I).LE.J2.AND.J2X(M+2,I).GE.J2)GO TO 33
            GO TO 28
31          IF(J2N(M+2,I).LE.J2.AND.J2X(M+2,I).GE.J2)GO TO 32
            GO TO 29
32          LJX(M,J,I)=M+2
            LJN(M,J,I)=M
            GO TO 34
33          LJX(M,J,I)=M+2
            LJN(M,J,I)=M+2
34          IF(M.GE.2)GO TO 35
            M=M+1
            IF(M.LE.LX)GO TO 27
            LJX(M,J,I)=-1
            LJN(M,J,I)=-1
35      CONTINUE
36  CONTINUE
!             BEGIN N-LOOP (SAMPLE THICKNESSES)
37  DO 46 NQ=1,NN
        N=NQ
!             BEGIN K-LOOP (ENERGIES)
        DO 45 KQ=1,NE
            PRINT*, "Running energy bin ", KQ, " of ", NE
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
38          CONTINUE
            GO TO (39,40,42,43),ITYPE
!           39 xx and 43 off to shut off multiple scattering (39, 40, 42, 43)
!           CYLINDRICAL SAMPLE CAPTURE
39          CALL MUSC(COMM,ZL,SGI,SGL,AP,UM,A,AB,BE,PE,TEFF,SPIN,NL,AA,RC,GG, &
                D,S,SI,R,EFF,J2X,J2N,SUMGJ,XN,E,SG,SP,G,GNR,GNRIN,DJL,SSCF,DSSCF, &
                N,K,I,L,J,NN,NE,NI,LX,RN,RX,ZH,PO,PS,DPO,DPS,ZP,SGM,STM,DSGM,DSTM, &
                TMC,DTMC,TAU,DTAU,SNC,XNSS,ITYPE,IHIST,LJX,LJN,JX2,JN2,JMX, RB)
!
            GO TO 41
!
!           SPHERICAL SHELL CAPTURE
40          CALL MUSS(COMM,ZL,SGI,SGL,AP,UM,A,AB,BE,PE,TEFF,SPIN,NL,AA,RC,GG, &
                D,S,SI,R,EFF,J2X,J2N,SUMGJ,XN,E,SG,SP,G,GNR,GNRIN,DJL,SSCF,DSSCF, &
                N,K,I,L,J,NN,NE,NI,LX,RN,RX,ZH,PO,PS,DPO,DPS,ZP,SGM,STM,DSGM,DSTM, &
                TMC,DTMC,TAU,DTAU,SNC,XNSS,ITYPE,IHIST,LJX,LJN,JX2,JN2,JMX, RB)
!
            XN(N)=RX(N)-RN(N)
41          IF(K.EQ.1) &
                WRITE(8,117)XN(N),E(K),SGM,STM,PO,PS,SSCF,SNC, ZH(K),ZP,RX(N), &
                TMC  ,DSGM,DSTM,DPO,DPS,DSSCF,DTMC
            IF(K.GT.1) &
                WRITE(8,118)      E(K),SGM,STM,PO,PS,SSCF,SNC, ZH(K),ZP,RX(N), &
                TMC  ,DSGM,DSTM,DPO,DPS,DSSCF,DTMC
                WRITE(12,312)E(K),STM,DSTM,SGM,DSGM,SSCF,DSSCF
            GO TO 44
!
!           TRANSMISSION
42          CALL MOCT(COMM,ZL,SGI,SGL,AP,UM,A,AB,BE,PE,TEFF,SPIN,NL,AA,RC,GG, &
                D,S,SI,R,EFF,J2X,J2N,SUMGJ,XN,E,SG,SP,G,GNR,GNRIN,DJL,SSCF,DSSCF, &
                N,K,I,L,J,NN,NE,NI,LX,RN,RX,ZH,PO,PS,DPO,DPS,ZP,SGM,STM,DSGM,DSTM, &
                TMC,DTMC,TAU,DTAU,SNC,XNSS,ITYPE,IHIST,LJX,LJN,JX2,JN2,JMX, RB)
!
            IF(K.EQ.1) &
                WRITE(8,119)XN(N),E(K),SGM,STM,TMC,TAU,ZH(K),ZP,DSGM,DSTM,DTMC,DTAU
            IF(K.GT.1) &
                WRITE(8,120)      E(K),SGM,STM,TMC,TAU,ZH(K),ZP,DSGM,DSTM,DTMC,DTAU
            GO TO 44
!
!           SELF-INDICATION
43          CALL MUSC(COMM,ZL,SGI,SGL,AP,UM,A,AB,BE,PE,TEFF,SPIN,NL,AA,RC,GG, &
                D,S,SI,R,EFF,J2X,J2N,SUMGJ,XN,E,SG,SP,G,GNR,GNRIN,DJL,SSCF,DSSCF, &
                N,K,I,L,J,NN,NE,NI,LX,RN,RX,ZH,PO,PS,DPO,DPS,ZP,SGM,STM,DSGM,DSTM, &
                TMC,DTMC,TAU,DTAU,SNC,XNSS,ITYPE,IHIST,LJX,LJN,JX2,JN2,JMX, RB)
!
            IF(K.EQ.1) &
                WRITE(8,121)XN(N),E(K),SGM,STM,PO,PS,SSCF,SNC, ZH(K),ZP,RX(N), &
                RN(N),DSGM,DSTM,DPO,DPS,DSSCF
            IF(K.GT.1) &
                WRITE(8,122)      E(K),SGM,STM,PO,PS,SSCF,SNC, ZH(K),ZP,RX(N), &
                RN(N),DSGM,DSTM,DPO,DPS,DSSCF
117         FORMAT(/ &
            1PE10.3,0PF7.3,1PE11.3,2E10.3,3E13.3,0P,2F10.0,1X,1P,2E10.3/ &
                     17X,1PE11.3,2E10.3,2E13.3,44X        ,1P E10.3/)
118         FORMAT(0PF17.3,1PE11.3,2E10.3,3E13.3,0P,2F10.0,1X,1P,2E10.3/ &
                     17X,1PE11.3,2E10.3,2E13.3,44X        ,1P E10.3/)
119         FORMAT(/ &
            1PE10.3,0PF7.3,1PE11.3,2E10.3,E13.4,0P,2F11.0/ &
                     17X,1PE11.3,2E10.3,E13.4/)
120         FORMAT(0PF17.3,1PE11.3,2E10.3,E13.4,0P,2F11.0/ &
                     17X,1PE11.3,2E10.3,E13.4/)
121         FORMAT(/ &
            1PE10.3,0PF7.3,1PE11.3,2E10.3,3E13.3,0P,2F10.0,1X,1P,2E10.3/ &
                     17X,1PE11.3,2E10.3,2E13.3/)
122         FORMAT(0PF17.3,1PE11.3,2E10.3,3E13.3,0P,2F10.0,1X,1P,2E10.3/ &
                     17X,1PE11.3,2E10.3,2E13.3/)
312         FORMAT(0PF10.3,6E11.3)
44          CONTINUE
45      CONTINUE
46  CONTINUE
    GO TO 1
999 STOP  
end program sesh
!
!--------------------------- END MAIN ---------------------------------
!
!end program sesh



