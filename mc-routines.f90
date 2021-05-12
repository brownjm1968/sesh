module mc_routines
    contains

subroutine musc(COMM,ZL,SGI,SGL,AP,UM,A,AB,BE,PE,TEFF,SPIN,NL,AA,RC,GG,D,S,SI, &
                R,EFF,J2X,J2N,SUMGJ,XN,E,SG,SP,G,GNR,GNRIN,DJL,SSCF,DSSCF,N,K, &
                I,L,J,NN,NE,NI,LX,RN,RX,ZH,PO,PS,DPO,DPS,ZP,SGM,STM,DSGM,DSTM, &
                TMC,DTMC,TAU,DTAU,SNC,XNSS,ITYPE,IHIST,LJX,LJN,JX2,JN2,JMX, RB)

    implicit none
!             MUSC  YIELDS THE MULTIPLE SCATTERING CORRECTION FOR A
!             CYLINDRICAL CAPTURE SAMPLE (MONTE CARLO)
!                XN(N)      CAPTURE SAMPLE THICKNESS
!                RX(N)      CAPTURE SAMPLE RADIUS
!                RB(N)      NEUTRON BEAM RADIUS
!                RN(N)      FILTER  SAMPLE THICKNESS
    
    real(4), allocatable, dimension(:) :: COMM,ZL,AP,UM,A,SS,STT,AB,BE,PE,TEFF, &
                              SPIN,AA,RC,XN,E,SG,SP,RN,RX,ZH,RB,ZST,ZSG,ZT,ZF0, &
                              GNL,PSM,POM,EE,STMC,SGMC,TM
    integer, allocatable, dimension(:) :: NL,JX2,JN2,JMX
    real(4), allocatable, dimension(:,:) :: SGI,GG,D,S,SI,R,EFF,SUMGJ,DOP
    integer, allocatable, dimension(:,:) :: J2X,J2N
    real(4), allocatable, dimension(:,:,:) :: SGL,G,GNR,GNRIN,DJL,SINE,COSE,PLE,GGE
    integer, allocatable, dimension(:,:,:) :: LJX,LJN
    real(4), allocatable, dimension(:,:,:,:) :: DE,GNE,GNINE

    real(4) :: arbitrary,AK,B1,B11,B12,B2,B23,B3,BBOTH,C0,CG,CT,CTH,CTHC,DD,  &
               DEN,DPO,DPS,DSGM,DSSCF,DSTM,DTAU,DTMC,DUMMY,DYN,EDD,EDGG,ETA,  &
               EXN1,FE,FP,GNINS,GNS,GT,H,HNEG,HPOS,O,PHI,PO,POMG,PS,PSMG,RHO, &
               RK,SAM,SC,SCC,SEM,SGM,SIM,SN,SNC,SOM,SQ,SSCF,ST,STH,STM,SUM,   &
               SUMF0,SUMSG,SUMST,SUMT,SYM,T,TAU,TF,TMC,TP,TX,U,UU,V,VL,VV,W,  &
               WG,WI,WN,X,XI,XNSS,Y,Z,ZP,ZSIR
    integer :: I,IHIST,ITYPE,IZ,J,JL,JX,K,KZ,L,L1,L2,LP,LX,M,MCHD,MDIV,MM,MP,N, &
               NC,NE,NH,NI,NN,NP,NS

    allocate(SS(100),STT(100),ZST(101),ZSG(101),ZT(101),ZF0(101),GNL(8),      &
             COSE(10,4,10),SINE(10,4,10),DE(10,8,4,10),GNE(10,8,4,10),        &
             GNINE(10,8,4,10),PLE(10,4,10),GGE(10,4,10),DOP(10,10),PSM(1000), &
             POM(1000),EE(10),STMC(1000),SGMC(1000),TM(1000))

    
    IF(IHIST.NE.1)GO TO 51
    do IZ=1,101
        ZST(IZ)=0.
        ZSG(IZ)=0.
        ZT (IZ)=0.
        ZF0(IZ)=0.
    end do  
51  CONTINUE
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
!   A1=(1.-EXP(-STT(K)*XN(N)))/STT(K)
!   A2=A1*SG(K)
!   IF(ITYPE.EQ.4)A1=A1*EXP(-STT(K)*RN(N))
    FE=1.-2./(A(1)*AA(1))
    do NC=1,10
        IF(NC.EQ.1)EE(NC)=E(K)
        IF(NC.GT.1)EE(NC)=EE(NC-1)*FE
        do I=1,NI
    ! ---- Doppler? -------------------------------------------------------
            DOP(NC,I)=SQRT(0.72531*A(I)/TEFF(I)/EE(NC))
    !
             CALL ENDEP(EE(NC),A(I),BE(I),PE(I),AP(I),UM(I),EDGG,EDD)
            
            LX=NL(I)
            do L=1,LX
                GGE(NC,L,I)=GG(L,I)*EDGG
                AK=RC(I)*SQRT(EE(NC)/AA(I))/143.92
        !
        ! ---- penetrabilities, etc..... --------------------------------------
                 CALL PEPS(AK,L,DUMMY,VL)
        !
                RK=AK*R(L,I)/RC(I)
        !
                 CALL PEPS(RK,L,PLE(NC,L,I),DUMMY)
        !
                SINE(NC,L,I)=SIN(2.*PLE(NC,L,I))
                COSE(NC,L,I)=COS(2.*PLE(NC,L,I))
                JX=(J2X(L,I)-J2N(L,I))/2+1
                do JL=1,JX
                    J=(J2N(L,I)-JN2(I))/2+JL
                    DE(NC,J,L,I)=DJL(JL,L,I)*EDD 
                    SQ=SQRT(EE(NC)*1000.)*VL*EDD
                    GNE  (NC,J,L,I)=GNR  (JL,L,I)*SQ
                    GNINE(NC,J,L,I)=GNRIN(JL,L,I)*SQ
                end do
            end do
        end do
    end do
      
!             BEGIN LOOP OF NEUTRON HISTORIES
    NH=ZH(K)
    NP=ZP
    do M=1,NH
!       PRINT HISTORY COUNTER
        IF(NH.LE.100)MDIV=10
        IF(NH.GT.100 .AND. NH.LE.1000)MDIV=100
        IF(NH.GT.1000)MDIV=1000
        MCHD=MOD(M,MDIV)
        IF(MCHD.EQ.0)THEN
        PRINT*, " History ", M, " of ", NH
        ENDIF
        
! *** INITIALIZE MONTE CARLO RUN ***
!
        U=0.                            ! INITIAL DIRECTION VECTOR IS NORMAL TO THE SAMPLE SURFACE
        V=0.
        W=1.
        IF(RB(N).GT.RX(N))RB(N)=RX(N)   ! IF THE BEAM IS LARGER THAN THE SAMPLE, SET THE BEAM SIZE EQUAL TO THE SAMPLE SIZE
        IF(RB(N).LE.(0.0))RB(N)=RX(N)   ! IF THE BEAM SIZE IS SET NEGATIVE OR EQUAL TO ZERO, SET THE BEAM SIZE EQUAL TO THE SAMPLE SIZE 
        arbitrary = 0.0d0
        X=RB(N)*SQRT(RANDOM(arbitrary))         ! SAMPLE ALONG THE BEAM RADIUS
        Y=RB(N)*SQRT(RANDOM(arbitrary))
        Z=0.
        TX=XN(N)
        NC=0
!
!
!             BEGIN NEUTRON HISTORY LOOP
! 
2       IF(NC.LT.10)NC=NC+1
!
!             FIND CROSS SECTIONS FOR COLLISION
        SN=0.
        SC=0.
        ST=SP(K)
        DO 3 I=1,NI
            JX=JMX(I)
            DO 83 J=1,JX
!                 LP=1: PARITY OF TARGET GROUND STATE
!                 LP=2: OPPOSITE PARITY
                DO 93 LP=1,2
                    L1=LJN(LP,J,I)
                    L2=LJX(LP,J,I)
                    IF(L1.LT.0)GO TO 93
                    JL=J-(J2N(L1,I)-JN2(I))/2
!                      4.617E3=2.605E3*1.772454
                    C0=4.617E3/EE(NC)*AA(I)*AB(I)*G(JL,L1,I)
!                           MP: COUNTER FOR RESONANCE PAIRS
                    MP=1
!                       FIND CENTRAL LEVEL INTERVAL AND ENERGY IN IT
!
                    CALL SPACE(DE(NC,J,L1,I),DD)
!
                    H=DD*RANDOM(arbitrary)
!                     FIND REDUCED NEUTRON WIDTH
17                  GNS=0.
                    GNINS=0.
                    DO L = L1,L2,2
!
                        CALL PORTER(GNE(NC,J,L,I),GNL(L))
                        
                        GNS=GNS+GNL(L)
                        GNINS=GNINS+GNINE(NC,J,L,I)
                    end do
                    GT=GNS+GNINS+GGE(NC,L1,I)
                    ETA=GT*DOP(NC,I)
                    XI=2.*ETA*H/GT
                    CT=C0*ETA/GT
                    CG=CT*GNS*GGE(NC,L1,I)/GT
!
                    CALL PFCN(XI,ETA,UU,VV,KZ)
!
                    SCC=CG*UU
                    SN=SN+SCC
!                     EFFECTIVE CAPTURE CROSS SECTION
                    SC=SC+SCC*EFF(L1,I)
!                     TOTAL CROSS SECTION
                    do L=L1,L2,2
                        ST=ST+CT*GNL(L)*(UU*COSE(NC,L,I)+VV*SINE(NC,L,I))
                    end do
                    IF(MP.EQ.1)GO TO 18
!                         CHECK IF SECOND MEMBER OF PAIR IS ALREADY INCLUDED
                    IF(H.LT.0.)GO TO 16
!
                    CALL WIGNER(DE(NC,J,L1,I),DD)
!
                    HNEG=HNEG-DD
                    H=HNEG
                    GO TO 17
!                         CHECK IF SECOND MEMBER OF PAIR IS ALREADY INCLUDED
18                  IF(H.LT.0.)GO TO 16
                    HPOS=H
                    HNEG=H-DD
                    H=HNEG
                    GO TO 17
!                         INCLUDE ANOTHER PAIR OF RESONANCES IF REQUIRED
16                  IF(MP.EQ.NP)GO TO 93
                    MP=MP+1
!
                    CALL WIGNER(DE(NC,J,L1,I),DD)
!
                    HPOS=HPOS+DD
                    H=HPOS
                    GO TO 17
93              CONTINUE
83          CONTINUE
3       CONTINUE
        SN=ST-SN
        T=EXP(-TX*ST)
!             INTERACTING FRACTION
        WI=1.-T
        WG=WI*SC
!               SURVIVING FRACTION
        WN=WI*SN/ST
        IF(NC.GT.1)GO TO 6
        TF=EXP(-RN(N)*ST)
        WG=WG*TF
        WN=WN*TF
!               CALCULATE HISTOGRAM VALUES IF REQUIRED
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
44      CONTINUE
        POMG=WG/ST
        SOM=SOM+POMG
        PSMG=0.
!
        CALL AVERT(TM(NS),T,SGMC(NS),SC,STMC(NS),ST,NS,1)
!
        GO TO 7
6       WG=WG/ST
        SUM=SUM+WG
        PSMG=PSMG+WG
!             PLAY RUSSIAN ROULETTE WITH SURVIVAL PROBABILITY WN
7       IF(RANDOM(arbitrary).GT.WN)GO TO 19
!             FREE PATH SAMPLING
        FP=-ALOG(1.-WI*RANDOM(arbitrary))/ST
!             COORDINATES OF COLLISION
        X=X+U*FP
        Y=Y+V*FP
        Z=Z+W*FP
!               SCATTERING ANGLES
        PHI=6.28318 *RANDOM(arbitrary)
        CTHC=2.*RANDOM(arbitrary)-1.
        CTH=(1.+A(1)*CTHC)/SQRT(1.+A(1)**2+2.*A(1)*CTHC)
        STH=SQRT(1.-CTH**2)
!               NEW DIRECTION COSINES
        RHO=SQRT(1.-W**2)
        IF(RHO.GT.0.001)GO TO 8
        U=W*STH*COS(PHI)
        V=W*STH*SIN(PHI)
        W=W*CTH
        GO TO 9
8       O=CTH*U+STH*(U*W*COS(PHI)-V*SIN(PHI))/RHO
        V=CTH*V+STH*(V*W*COS(PHI)+U*SIN(PHI))/RHO
        W=CTH*W-STH*COS(PHI)*RHO
        U=O
!             PATH LENGTH TO SAMPLE SURFACE
9       B1=U**2+V**2
        B2=U*X+V*Y
        B3=RX(N)**2-X**2-Y**2
        IF(B2.LT.0.)GO TO 10
        TX=B3/(SQRT(B2**2+B1*B3)+B2)
        GO TO 11
10      TX=(SQRT(B2**2+B1*B3)-B2)/B1
11      IF(W.LT.-0.000001)GO TO 12
        IF(W.GT.+0.000001)GO TO 13
        GO TO 2
12      TP=-Z/W
        GO TO 14
13      TP=(XN(N)-Z)/W
14      IF(TP.LT.TX)TX=TP
        GO TO 2
19      SNC=SNC+FLOAT(NC)
!
        CALL AVERP(POM(MM),POMG,PSM(MM),PSMG,MM,1)
end do
!
!             END OF HISTORY LOOP
!             MONTE CARLO RESULTS
!
    CALL AVERP(POM(MM),POMG,PSM(MM),PSMG,MM,2)
    CALL AVERT(TM(NS),T,SGMC(NS),SC,STMC(NS),ST,NS,2)
!
    TMC=0.
    SGM=0.
    STM=0.
    do M=1,NS
        TMC=TMC+TM(M)
        SGM=SGM+SGMC(M)
        STM=STM+STMC(M)
    end do
    TMC=TMC/FLOAT(NS)
    SGM=SGM/FLOAT(NS)
    STM=STM/FLOAT(NS)
    PO =SOM/ZH(K)
    PS =SUM/ZH(K)
    SNC=SNC/ZH(K)
    ZSIR=ZSIR/(ZH(K)*SGM)
!             STATISTICAL UNCERTAINTIES
    do M=1,NS
        DTMC=DTMC+(TM(M)-TMC)**2
        SIM=SIM+(SGMC(M)-SGM)**2
        SYM=SYM+(STMC(M)-STM)**2
    end do
    B11=0.
    B12=0.
    B23=0.
    BBOTH=0.
    do M=1,MM
        BBOTH=BBOTH+(POM(M)+PSM(M)-PO-PS)*(SGMC(M)-SGM)
        IF(ITYPE.EQ.4)B12=B12+(POM(M)+PSM(M)-PO-PS)*(STMC(M)-STM)
        IF(ITYPE.EQ.4)B23=B23+(STMC(M)-STM)*(SGMC(M)-SGM)
                      B11=B11+(POM(M)-PO)*(PSM(M)-PS)
        SAM=SAM+(POM(M)-PO)**2
        SEM=SEM+(PSM(M)-PS)**2
    end do
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
!             SELF-SHIELDING CORRECTION FACTOR AND ITS PERCENT ERROR
! set PO = 0 to remove multiple scattering contribution
!      PO=0
    SSCF=(PO+PS)/(XN(N)*SGM)
    IF(ITYPE.EQ.4)EXN1=EXP(RN(N)*STM)
    IF(ITYPE.EQ.4)SSCF=SSCF*EXN1
    IF(ITYPE.EQ.1)DSSCF=SQRT(B11/(XN(N)*SGM)**2-2.*BBOTH*SSCF/(XN(N)*SGM**2)+(SSCF/SGM)**2*DSGM**2)/SSCF*100.0
    IF(ITYPE.EQ.4)DSSCF=SQRT((EXN1/(XN(N)*SGM))**2*B11+2.*(EXN1/(XN(N)*SGM))* &
                  RN(N)*SSCF*B12-2.*(EXN1/(XN(N)*SGM))*SSCF/SGM*BBOTH+(RN(N)* &
                  SSCF*DSTM)**2-2.*RN(N)*SSCF**2/SGM*B23+(SSCF/SGM*DSGM)**2)/SSCF*100.0
!           WRITE HISTOGRAM VALUES IF REQUIRED
    IF(IHIST.NE.1)GO TO 53
    SUMST=0.
    SUMSG=0.
    SUMT =0.
    SUMF0=0.
    do IZ=1,101
        SUMST=SUMST+ZST(IZ)
        SUMSG=SUMSG+ZSG(IZ)
        SUMT =SUMT + ZT(IZ)
        SUMF0=SUMF0+ZST(IZ)
    end do 
    do IZ=1,101
        ZST(IZ)=ZST(IZ)/SUMST
        ZSG(IZ)=ZSG(IZ)/SUMSG
        ZT(IZ)=ZT(IZ)/SUMT
        ZF0(IZ)=ZF0(IZ)/SUMF0
    end do
    WRITE(8,198)ZSIR
198 FORMAT(/1PE15.3/)
    WRITE(8,199)(IZ,ZST(IZ),ZSG(IZ),ZT(IZ),ZF0(IZ),IZ=1,101)
199 FORMAT(1H1,///,' NO.      RELATIVE FREQUENCIES OF                 ST     SG     T      F0       ',(I4,0PF8.4,3F7.4))
53  CONTINUE
    RETURN
end subroutine musc
!
subroutine muss(COMM,ZL,SGI,SGL,AP,UM,A,AB,BE,PE,TEFF,SPIN,NL,AA,RC,GG,D,S,SI, &
                R,EFF,J2X,J2N,SUMGJ,XN,E,SG,SP,G,GNR,GNRIN,DJL,SSCF,DSSCF,N,K, &
                I,L,J,NN,NE,NI,LX,RN,RX,ZH,PO,PS,DPO,DPS,ZP,SGM,STM,DSGM,DSTM, &
                TMC,DTMC,TAU,DTAU,SNC,XNSS,ITYPE,IHIST,LJX,LJN,JX2,JN2,JMX, RB)

    implicit none
!             MUSS YIELDS THE MULTIPLE SCATTERING CORRECTION FOR A
!             SPHERICAL SHELL (MONTE CARLO)

    real(4), allocatable, dimension(:) :: COMM,ZL,AP,UM,A,SS,STT,AB,BE,PE,TEFF, &
                              SPIN,AA,RC,XN,E,SG,SP,RN,RX,ZH,RB,ZST,ZSG,ZT,ZF0, &
                              GNL,EE,STMC,SGMC,TM,PSM,POM
    integer, allocatable, dimension(:) :: NL,JX2,JN2,JMX
    real(4), allocatable, dimension(:,:) :: SGI,GG,D,S,SI,R,EFF,SUMGJ,DOP
    integer, allocatable, dimension(:,:) :: J2X,J2N
    real(4), allocatable, dimension(:,:,:) :: SGL,G,GNR,GNRIN,DJL,SINE,COSE,PLE,GGE
    integer, allocatable, dimension(:,:,:) :: LJX,LJN
    real(4), allocatable, dimension(:,:,:,:) :: DE,GNE,GNINE

    real(4) :: arbitrary,A1,A2,AK,B1,B11,B12,B2,B23,B3,BBOTH,C0,CG,CT,CTH,CTHC, &
               DD,DEN,DPO,DPS,DSGM,DSSCF,DSTM,DTAU,DTMC,DUMMY,DYN,EDD,EDGG,ETA, &
               EXN1,FE,FP,GNINS,GNS,GT,H,HNEG,HPOS,O,PHI,PO,POMG,PS,PSMG,RHO,   &
               RK,SAM,SC,SCC,SEM,SGM,SIM,SN,SNC,SOM,SQ,SSCF,ST,STH,STM,SUM,     &
               SUMF0,SUMSG,SUMST,SUMT,SYM,T,TAU,TF,TMC,TP,TX,U,UU,V,VL,VV,W,WG, &
               WI,WN,X,XI,XNSS,Y,Z,ZP,ZSIR,Q
    integer :: I,IHIST,ITYPE,IZ,J,JL,JX,K,KZ,L,L1,L2,LP,LX,M,MCHD,MDIV,MM,MP,N, &
               NC,NE,NH,NI,NN,NP,NS

    allocate(SS(100),STT(100),ZST(101),ZSG(101),ZT(101),ZF0(101),GNL(8),        &
             COSE(10,4,10),SINE(10,4,10),EE(10),STMC(1000),SGMC(1000),TM(1000), &
             PSM(1000),POM(1000),PLE(10,4,10),GGE(10,4,10),DOP(10,10),          &
             DE(10,8,4,10),GNE(10,8,4,10),GNINE(10,8,4,10))
    
    IF(IHIST.NE.1)GO TO 51
    DO 50 IZ=1,101
    ZST(IZ)=0.
    ZSG(IZ)=0.
    ZT (IZ)=0.
50  ZF0(IZ)=0.
51  CONTINUE
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
!
        CALL ENDEP(EE(NC),A(I),BE(I),PE(I),AP(I),UM(I),EDGG,EDD)
!
        GGE(NC,L,I)=GG(L,I)*EDGG
        RK=RC(I)*SQRT(EE(NC)/AA(I))/143.92
!
        CALL PEPS(RK,L,PLE(NC,L,I),VL)
!
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
20  CONTINUE
!             BEGIN LOOP OF NEUTRON HISTORIES
    NH=ZH(K)
    NP=ZP
    DO 1 M=1,NH
!       PRINT HISTORY COUNTER
        IF(NH.LE.100)MDIV=10
        IF(NH.GT.100 .AND. NH.LE.1000)MDIV=100
        IF(NH.GT.1000)MDIV=1000
        MCHD=MOD(M,MDIV)
        IF(MCHD.EQ.0)THEN
        PRINT*, " History ", M, " of ", NH
        ENDIF
!               INITIALIZE
        Q=0.
        CTH=-1.
        STH=0.
        Z=RN(N)
        TX=XNSS
        W=1.
        NC=0
!            BEGIN OF COLLISION LOOP
2       IF(NC.LT.10)NC=NC+1
!            FIND CROSS SECTIONS FOR COLLISION
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
!                    4.617E3=2.605E3*1.772454
                    JL=J-(J2N(L1,I)-JN2(I))/2
                    C0=4.617E3/EE(NC)*AA(I)*AB(I)*G(JL,L1,I)
                    MP=1
!                    FIND CENTRAL LEVEL INTERVAL AND ENERGY IN IT
!
                    CALL SPACE(DE(NC,J,L1,I),DD)
!
                    H=DD*RANDOM(arbitrary)
17                  GNS=0.
                    GNINS=0.
                    DO 40 L=L1,L2,2
!                     FIND REDUCED NEUTRON WIDTH
!
                    CALL PORTER(GNE(NC,J,L,I),GNL(L))
!
                    GNS=GNS+GNL(L)
40                  GNINS=GNINS+GNINE(NC,J,L,I)
                    GT=GNS+GNINS+GGE(NC,L1,I)
                    ETA=GT*DOP(NC,I)
                    XI=2.*ETA*H/GT
                    CT=C0*ETA/GT
                    CG=CT*GNS*GGE(NC,L1,I)/GT
!
                    CALL PFCN(XI,ETA,UU,VV,KZ)
!
!                     CAPTURE CROSS SECTION
                    SCC=CG*UU
                    SN=SN+SCC
!                     EFFECTIVE CAPTURE CROSS SECTION
                    SC=SC+SCC*EFF(L1,I)
!                     TOTAL CROSS SECTION
                    DO 41 L=L1,L2,2
41                  ST=ST+CT*GNL(L)*(UU*COSE(NC,L,I)+VV*SINE(NC,L,I))
                    IF(MP.EQ.1)GO TO 18
!                     CHECK IF SECOND MEMBER OF PAIR IS ALREADY INCLUDED
                    IF(H.LT.0.)GO TO 16
!
                    CALL WIGNER(DE(NC,J,L1,I),DD)
!
                    HNEG=HNEG-DD
                    H=HNEG
                    GO TO 17
!                     CHECK IF SECOND MEMBER OF PAIR IS ALREADY INCLUDED
18                  IF(H.LT.0.)GO TO 16
                    HPOS=H
                    HNEG=H-DD
                    H=HNEG
                    GO TO 17
!                     INCLUDE ANOTHER PAIR OF RESONANCES IF REQUIRED
16                  IF(MP.EQ.NP)GO TO 93
                    MP=MP+1
!
                    CALL WIGNER(DE(NC,J,L1,I),DD)
!
                    HPOS=HPOS+DD
                    H=HPOS
                    GO TO 17
93              CONTINUE
83          CONTINUE
3       CONTINUE
        SN=ST-SN
        T=EXP(-TX*ST)
!             INTERACTING FRACTION
        WI=1.-T
!             CAPTURED AND DETECTED FRACTION
        WG=WI*SC
        IF(ITYPE.EQ.4)WG=WG*EXP(-RN(N)*ST)
!             SURVIVING FRACTION
        WN=WI*SN/ST
        IF(NC.GT.1)GO TO 6
!             CALCULATE HISTOGRAM VALUES IF REQUIRED
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
52      CONTINUE
        POMG=WG*(1./ST-A1/WI)+A2
        SOM=SOM+POMG
        PSMG=0.
!
        CALL AVERT(TM(NS),T,SGMC(NS),SC,STMC(NS),ST,NS,1)
!
        GO TO 7
6       WG=WG/ST
        PSMG=PSMG+WG
        SUM=SUM+WG
!             PLAY RUSSIAN ROULETTE WITH SURVIVAL PROBABILITY WN
7       IF(RANDOM(arbitrary).GT.WN)GO TO 19
!             CALCULATE RADIAL COORDINATE (Z) OF NC-TH COLLISION
        T=-ALOG(1.-WI*RANDOM(arbitrary))/ST
        IF(T.GT.(Z*CTH-Q))  T=T+2.*Q
        Z=SQRT(Z**2+  T**2-2.*Z*  T*CTH)
!             ANGLE WITH RESPECT TO RADIUS AFTER COLLISION
        CTH=2.*RANDOM(arbitrary)-1.
        STH=SQRT(1.-CTH**2)
!             MAXIMUM PATH LENGTH IN MATTER
        Q=0.
        IF(CTH.GT.SQRT(1.-(RN(N)/Z)**2))Q=SQRT(RN(N)**2-(Z*STH)**2)
        TX=Z*CTH+SQRT(RX(N)**2-(Z*STH)**2)-2.*Q
        GO TO 2
19      SNC=SNC+FLOAT(NC)
!
1   CALL AVERP(POM(MM),POMG,PSM(MM),PSMG,MM,1)
!
!             END OF HISTORY LOOP
!             MONTE CARLO RESULTS
!
    CALL AVERP(POM(MM),POMG,PSM(MM),PSMG,MM,2)
    CALL AVERT(TM(NS),T,SGMC(NS),SC,STMC(NS),ST,NS,2)
!
    TMC=0.
    SGM=0.
    STM=0.
    DO 22 M=1,NS
    TMC=TMC+TM(M)
    SGM   =SGM   +SGMC(M)
22  STM   =STM   +STMC(M)
    SGM   =SGM   /FLOAT(NS)
    TMC=TMC/FLOAT(NS)
    STM   =STM   /FLOAT(NS)
    PO=SOM/ZH(K)
    PS=SUM/ZH(K)
    SNC=SNC/ZH(K)
!             STATISTICAL UNCERTAINTIES
    DO 23 M=1,NS
    DTMC=DTMC+(TM(M)-TMC)**2
    SIM=SIM+(SGMC(M)-SGM   )**2
23  SYM=SYM+(STMC(M)-STM   )**2
    B11=0.
    B12=0.
    B23=0.
    BBOTH=0.
    DO 15 M=1,MM
    BBOTH=BBOTH+(POM(M)+PSM(M)-PO-PS)*(SGMC(M)-SGM)
                  B11=B11+(POM(M)-PO)*(PSM(M)-PS)
    SAM=SAM+(POM(M)-PO     )**2
15  SEM=SEM+(PSM(M)-PS     )**2
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
!             SELF-SHIELDING CORRECTION FACTOR AND ITS PERCENT ERROR
    SSCF=1.-PO-PS
    DSSCF=SQRT(DPO**2+DPS**2)
!             WRITE HISTOGRAM VALUES IF REQUIRED
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
54  CONTINUE
    DO 55 IZ=1,101
    ZST(IZ)=ZST(IZ)/SUMST
    ZSG(IZ)=ZSG(IZ)/SUMSG
    ZT(IZ)=ZT(IZ)/SUMT
    ZF0(IZ)=ZF0(IZ)/SUMF0
55  CONTINUE
    WRITE(8,199)(IZ,ZST(IZ),ZSG(IZ),ZT(IZ),ZF0(IZ),IZ=1,101)
199 FORMAT(1H1,///,' NO.    TALLY FOR HISTOGRAM OF                  ST     SG     T      F0     ',(I4,0PF8.4,3F7.4))
53  CONTINUE
    RETURN
end subroutine muss
!
subroutine moct(COMM,ZL,SGI,SGL,AP,UM,A,AB,BE,PE,TEFF,SPIN,NL,AA,RC,GG,D,S,SI, &
                R,EFF,J2X,J2N,SUMGJ,XN,E,SG,SP,G,GNR,GNRIN,DJL,SSCF,DSSCF,N,K, &
                I,L,J,NN,NE,NI,LX,RN,RX,ZH,PO,PS,DPO,DPS,ZP,SGM,STM,DSGM,DSTM, &
                TMC,DTMC,TAU,DTAU,SNC,XNSS,ITYPE,IHIST,LJX,LJN,JX2,JN2,JMX, RB)

    implicit none
!             MOCT  YIELDS THE SELF-SHIELDING CORRECTION FACTOR AND
!             TRANSMISSION FOR A CYLINDRICAL SAMPLE

    real(4), allocatable, dimension(:) :: COMM,ZL,AP,UM,A,SS,STT,AB,BE,PE,TEFF, &
                              SPIN,AA,RC,XN,E,SG,SP,RN,RX,ZH,RB,ZST,ZSG,ZT,ZF0, &
                              GNL,TMCG,TN,SGMC,STMC,DOP
    integer, allocatable, dimension(:) :: NL,JX2,JN2,JMX
    real(4), allocatable, dimension(:,:) :: SGI,GG,D,S,SI,R,EFF,SUMGJ,PL,SINE,  &
                                            COSE,GGE,VL,EDGG,EDD
    integer, allocatable, dimension(:,:) :: J2X,J2N
    real(4), allocatable, dimension(:,:,:) :: SGL,G,GNR,GNRIN,DJL,GN,GNIN,DE
    integer, allocatable, dimension(:,:,:) :: LJX,LJN

    real(4) :: arbitrary,AK,AVT,AVTS,B1,B11,B12,B2,B23,B3,BBOTH,C0,CG,CT,CTH, &
               CTHC,DD, DEN,DPO,DPS,DSGM,DSSCF,DSTM,DT,DTS,DTS2,DTAU,DTMC,    &
               DUMMY,DYN,EE,ETA,EXN1,FE,FP,GNE,GNINE,GNINS,GNS,GT,H,HNEG,     &
               HPOS,O,PHI,PLE,PO,POM,POMG,PS,PSM,PSMG,RHO,RK,SAM,SC,SCC,SEM,  &
               SGM,SIM,SN,SNC,SOM,SQ,SSCF,ST,STH,STM,SUM,SUMF0,SUMSG,SUMST,   &
               SUMT,SYM,T,TAU,TF,TM,TMC,TP,TX,U,UU,V,VV,W,WG,WI,WN,X,XI,XNSS, &
               Y,Z,ZP,ZSIR,TNUM,TS,TSSM,TSUM,A1,A2,Q
    integer :: I,IHIST,IT,ITYPE,IZ,J,JL,JX,K,KZ,L,L1,L2,LP,LX,M,MCHD,MDIV,MM, &
               MP,N,NC,NE,NH,NHIST,NI,NN,NP,NS

    allocate(ZST(101),ZSG(101),ZT(101),TN(101),GNL(8),COSE(4,10),SINE(4,10), &
             GN(8,4,10),GNIN(8,4,10),PL(4,10),TMCG(1000),SGMC(1000),         &
             STMC(1000),DOP(10),GGE(4,10),VL(4,10),EDGG(10,100),EDD(10,100), &
             DE(8,4,10))

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
2   ZT(IT)=0.
    DO 20 I=1,NI
    DOP(I)=SQRT(0.72531*A(I)/TEFF(I)/E(K))
!
        CALL ENDEP(E(K),A(I),BE(I),PE(I),AP(I),UM(I),EDGG(I,K),EDD(I,K))
!
    LX=NL(I)
    DO 20 L=1,LX
    GGE(L,I)=GG(L,I)*EDGG(I,K)
    AK=RC(I)*SQRT(E(K)/AA(I))/143.92
!
        CALL PEPS(AK,L,DUMMY,VL(L,I))
!
    RK=AK*R(L,I)/RC(I)
!
        CALL PEPS(RK,L,PL(L,I),DUMMY)
!
    SINE(L,I)=SIN(2.*PL(L,I))
    COSE(L,I)=COS(2.*PL(L,I))
    JX=(J2X(L,I)-J2N(L,I))/2+1
    DO 20 JL=1,JX
    J=(J2N(L,I)-JN2(I))/2+JL
    DE(J,L,I)=DJL(JL,L,I)*EDD(I,K)
    GN  (J,L,I)=GNR  (JL,L,I)*SQRT(E(K)*1000.)*VL(L,I)*EDD(I,K)
    GNIN(J,L,I)=GNRIN(JL,L,I)*SQRT(E(K)*1000.)*VL(L,I)*EDD(I,K)
20  CONTINUE
!             BEGIN LOOP OF NEUTRON HISTORIES
    NH=ZH(K)
    NP=ZP
    DO 1 M=1,NH
    NHIST=0
    NHIST=NHIST+1
!   PRINT HISTORY COUNTER
    IF(NH.LE.100)MDIV=10
    IF(NH.GT.100 .AND. NH.LE.1000)MDIV=100
    IF(NH.GT.1000)MDIV=1000
    MCHD=MOD(M,MDIV)
    IF(MCHD.EQ.0)THEN
    PRINT*, " History ", M, " of ", NH
    ENDIF
!           FIND CROSS SECTIONS FOR COLLISION
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
!           4.617E3=2.605E3*1.772454
    JL=J-(J2N(L1,I)-JN2(I))/2
    C0=4.617E3/E(K)*AA(I)*AB(I)*G(JL,L1,I)
    MP=1
!            FIND CENTRAL LEVEL INTERVAL AND ENERGY IN IT
!
        CALL SPACE(DE(J,L1,I),DD)
!
    H=DD*RANDOM(arbitrary)
17  GNS=0.
    GNINS=0.
    DO 40 L=L1,L2,2
!             FIND REDUCED NEUTRON WIDTH
!
        CALL PORTER(GN(J,L,I),GNL(L))
!
    GNS=GNS+GNL(L)
40  GNINS=GNINS+GNIN(J,L,I)
    GT=GNS+GNINS+GGE(L1,I)
    ETA=GT*DOP(I)
    XI=2.*ETA*H/GT
    CT=C0*ETA/GT
    CG=CT*GNS*GGE(L1,I)/GT
!
        CALL PFCN(XI,ETA,UU,VV,KZ)
!
!           CAPTURE CROSS SECTION
    SCC=CG*UU
    SN=SN+SCC
!           EFFECTIVE CAPTURE CROSS SECTION
    SC=SC+SCC*EFF(L1,I)
!             TOTAL CROSS SECTION
    DO 41 L=L1,L2,2
41  ST=ST+CT*GNL(L)*(UU*COSE(L,I)+VV*SINE(L,I))
    IF(MP.EQ.1)GO TO 18
!             CHECK IF SECOND MEMBER OF PAIR IS ALREADY INCLUDED
    IF(H.LT.0.)GO TO 16
!
        CALL WIGNER(DE(J,L1,I),DD)
!
    HNEG=HNEG-DD
    H=HNEG
    GO TO 17
!            CHECK IF SECOND MEMBER OF PAIR IS ALREADY INCLUDED
18  IF(H.LT.0.)GO TO 16
    HPOS=H
    HNEG=H-DD
    H=HNEG
    GO TO 17
!             INCLUDE ANOTHER PAIR OF RESONANCES IF REQUIRED
16  IF(MP.EQ.NP)GO TO 93
    MP=MP+1
!
        CALL WIGNER(DE(J,L1,I),DD)
!
    HPOS=HPOS+DD
    H=HPOS
    GO TO 17
93  CONTINUE
83  CONTINUE
3   CONTINUE
    SN=ST-SN
    T=EXP(-XN(N)*ST)
    TS=T*ST
    write(99,'(2E12.5)')ST, T
    IF(IHIST.NE.1)GO TO 52
    IZ=20.*ST/STT(K)+1.
    IF(IZ.LT.1)  IZ=1
    IF(IZ.GT.101)IZ=101
    ZST(IZ)=ZST(IZ)+1.
    IZ=20.*SC/SG(K)+1.
    IF(IZ.LT.1)  IZ=1
    IF(IZ.GT.101)IZ=101
    ZSG(IZ)=ZSG(IZ)+1.
52  IZ=2.*T/DT+1.
    IF(IZ.LT.1  )IZ=1
    IF(IZ.GT.101)IZ=101
    ZT(IZ)=ZT(IZ)+1.
    TNUM=TNUM+1.
    TSUM=TSUM+T
    TSSM=TSSM+TS
!
        CALL AVERT(TMCG(MM),T,SGMC(MM),SC,STMC(MM),ST,MM,1)
!
1   CONTINUE
!             END OF HISTORY LOOP
    AVT=TSUM/TNUM
    AVTS=TSSM/TNUM
!             MONTE CARLO RESULTS
!
        CALL AVERT(TMCG(MM),T,SGMC(MM),SC,STMC(MM),ST,MM,2)
!
    TMC=0.
    SGM=0.
    STM=0.
    DO 22 M=1,MM
    SGM=SGM+SGMC(M)
    STM=STM+STMC(M)
22  TMC=TMC+TMCG(M )
    TMC=TMC/FLOAT(MM)
    SGM=SGM/FLOAT(MM)
    STM=STM/FLOAT(MM)
!           STATISTICAL UNCERTAINTY
    DTMC=0.
    DTS=0.
    DO 23 M=1,MM
    SIM=SIM+(SGMC(M)-SGM)**2
    SYM=SYM+(STMC(M)-STM)**2
    DTS=DTS+(STMC(M)-STM)*(TMCG(M)-TMC)
23  DTMC=DTMC+(TMCG(M)-TMC)**2
    DEN=FLOAT(MM*(MM-1))
    DSGM=SQRT(SIM/DEN)
    DSTM=SQRT(SYM/DEN)
    DTMC=SQRT(DTMC/DEN)
    DTS2=DTS/DEN
    TAU=EXP(XN(N)*STM)*TMC
    DTAU=(XN(N)*DSTM)**2+2.*XN(N)*DTS2/TMC+(DTMC/TMC)**2
    IF(DTAU.LT.0.0)DTAU=0.0
    DTAU=100.*SQRT(DTAU)
!           WRITE HISTOGRAM VALUES IF REQUIRED
    IF(IHIST.NE.1)GO TO 53
    SUMST=0.
    SUMSG=0.
    SUMT =0.
    DO 54 IZ=1,101
    SUMST=SUMST+ZST(IZ)
    SUMSG=SUMSG+ZSG(IZ)
    SUMT =SUMT + ZT(IZ)
54  CONTINUE
    DO 55 IZ=1,101
    ZST(IZ)=ZST(IZ)/SUMST
    ZSG(IZ)=ZSG(IZ)/SUMSG
    ZT(IZ)=ZT(IZ)/SUMT
55 CONTINUE
    WRITE(8,199)(IZ,ZST(IZ),ZSG(IZ),ZT(IZ),IZ=1,101)
199 FORMAT(1H1,///,' NO.    TALLY FOR HISTOGRAM OF                ST     SG     T           ',(I4,0PF8.4,2F7.4))
!
!C     CALL PLOT(TN,ZT,100,3,3,1,1,1,1,100.,0.,0.125,100.,0.,0.125,'DISTR
!C   1IBUTION OF TRANSMISSION VALUES',1)
53  CONTINUE
    RETURN
end subroutine moct

!
subroutine peps(RK,L,PL,VL)
!             PEPS CALCULATES CENTRIFUGAL-BARRIER PENETRABILITIES, VL,
!             AND HARD-SPHERE PHASE SHIFTS, PL.
    
    implicit none

    real(4), intent(in) :: RK
    integer, intent(in) :: L
    real(4), intent(out) :: PL
    real(4) :: VL

    GO TO (1,2,3,4),L
1   PL=RK
    VL=1.
    RETURN
2   PL=RK-ATAN(RK)
    VL=RK**2/(RK**2+1.)
    RETURN
3   PL=RK-ATAN(3.*RK/(3.-RK**2))
    VL=RK**4/(RK**4+3.*RK**2+9.)
    RETURN
4   PL=RK-ATAN(RK*(15.-RK**2)/3./(5.-2.*RK**2))
    VL=RK**6/(RK**6+6.*RK**4+45.*RK**2+225.)
    RETURN
end subroutine peps
!
subroutine endep(E,A,B,P,AP,UM,EDGG,EDD)
!             ENDEP CALCULATES THE FACTORS WHICH DESCRIBE THE ENERGY
!             DEPENDENCE OF LEVEL SPACINGS AND RADIATION WIDTHS

    implicit none

    real(4), intent(in) :: A,B,P,AP,UM
    real(4), intent(inout) :: E

    real(4) :: Y1,Y2,A1,E1,E2,ER,F1,F2,GR,SUM1,SUM2,U,U1X,U2X,UB,U1,U2,EDGG,EDD
    integer :: N

    DIMENSION Y1(21),Y2(21)
!    PRINT*, E, A, B, P, AP, UM, EDGG, EDD
    E=0.001*E
    A1=A+1.
!           LEVEL DENSITY FACTOR
    U =B+E-P
    UB=B-P
    EDD=RHO(UB,A,P,AP,UM)/RHO(U,A,P,AP,UM)
!           RADIATION WIDTH FACTOR
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
1   CONTINUE
!
    CALL SIMP(Y1,1,21,SUM1)
    CALL SIMP(Y2,1,21,SUM2)
!
    EDGG=EDD*SUM1/SUM2
    E=1000.*E
    RETURN
end subroutine endep
!
function rho(U,A,PE,AP,UM)
!
!    GILBERT-CAMERON COMPOSITE LEVEL DENSITY
!    (ALL SPINS AND PARITIES)
!
!    Parameters: 
!         U:   EXCITATION ENERGY (MEV)
!         A:   ATOMIC WEIGHT OF TARGET NUCLEUS (AMU)
!         PE:  PAIRING CORRECTION (MEV)
!         AP:  A-PARAMETER (1/MEV)
!         UM:  MATCHING ENERGY (MEV)
!
    implicit none

    real(4), intent(in) :: U,A,PE,AP,UM

    real(4) :: rho,VARJ2
    real(4) :: AI,AU4,T,UEFF

    AI=FLOAT(INT(A+1.5))
    IF(U.GT.UM)GO TO 1
!   LOW ENERGIES: CONSTANT-TEMPERATURE FORMULA
    UEFF=UM-PE
    AU4=AP*UEFF*4.
    VARJ2=0.1460*SQRT(AU4)*CBRT(INT(AI)*INT(AI))
    T=1./(SQRT(AP/UEFF)-1.5/UEFF)
    RHO=(AP/3.)*EXP(SQRT(AU4))/AU4**1.25*SQRT(2./VARJ2)
    RHO=RHO*EXP((U-UM)/T)
    RETURN
!   HIGH ENERGIES: FERMI-GAS FORMULA
1   AU4=AP*(U-PE)*4.
    VARJ2=0.1460*SQRT(AU4)*CBRT(INT(AI)*INT(AI))
    RHO=(AP/3.)*EXP(SQRT(AU4))/AU4**1.25*SQRT(2./VARJ2)
!    PRINT*, RHO
    RETURN
end function rho
!
subroutine aver(SGMC,SC,STMC,ST,NS,I)

    implicit none

    real(4), intent(inout) :: SGMC, STMC
    real(4), intent(in) :: SC, ST
    integer, intent(inout) :: NS
    integer, intent(in) :: I

    real(4) :: corr,W
    integer :: N,J
    
    DATA N,J/100,0/
!        N   GROUP SIZE
!        I=1 SUPPLY DATA TO BE GROUPED
!        I=2 FINISH LAST GROUP,WHICH MAY CONTAIN LESS THAN N MEMBERS
    W=1./FLOAT(N)
    GO TO (1,3),I
1   IF(J.GT.0)GO TO 2
    SGMC=0.
    STMC=0.
    J=1
2   SGMC=SGMC+SC
    STMC=STMC+ST
    J=J+1
    IF(J.LE.N)RETURN
    J=0
    SGMC=SGMC*W
    STMC=STMC*W
    NS=NS+1
    RETURN
3   IF(J.EQ.0)GO TO 4
    CORR=1.0/FLOAT(J-1)
    SGMC=SGMC*CORR
    STMC=STMC*CORR
    J=0
    RETURN
4   NS=NS-1
5   RETURN
end subroutine aver
!
subroutine averp(PZ,PZG,PS,PSG,MM,I)

    implicit none

    real(4), intent(inout) :: PZ,PS
    real(4), intent(in) :: PZG,PSG
    integer, intent(inout) :: MM
    integer, intent(in) :: I

    real(4) :: W,corr
    integer :: N,J

    DATA N,J/100,0/
!        N   GROUP SIZE
!        I=1 SUPPLY DATA TO BE GROUPED
!        I=2 FINISH LAST GROUP,WHICH MAY CONTAIN LESS THAN N MEMBERS
    W=1./FLOAT(N)
    GO TO (1,3),I
1   IF(J.GT.0)GO TO 2
    PZ=0.
    PS=0.
    J=1
2   PZ=PZ+PZG
    PS=PS+PSG
    J=J+1
    IF(J.LE.N)RETURN
    J=0
    PZ=PZ*W
    PS=PS*W
    MM=MM+1
    RETURN
3   IF(J.EQ.0)GO TO 4
    CORR=1.0/FLOAT(J-1)
    PZ=PZ*CORR
    PS=PS*CORR
    J=0
    RETURN
4   MM=MM-1
5   RETURN
end subroutine averp
!
subroutine avert(TM,T,SGM,SC,STM,ST,MM,I)

    implicit none

    real(4), intent(in) :: T,SC,ST
    integer, intent(in) :: I
    real(4), intent(inout) :: TM
    integer, intent(inout) :: MM

    real(4) :: corr,W,SGM,STM
    integer :: N,J

    DATA N,J/100,0/
!        N   GROUP SIZE
!        I=1 SUPPLY DATA TO BE GROUPED
!        I=2 FINISH LAST GROUP,WHICH MAY CONTAIN LESS THAN N MEMBERS
    W=1./FLOAT(N)
    GO TO (1,3),I
1   IF(J.GT.0)GO TO 2
    TM=0.
    SGM=0.
    STM=0.
    J=1
2   TM=TM+T
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
3   IF(J.EQ.0)GO TO 4
    CORR=1.0/FLOAT(J-1)
    TM=TM*CORR
    SGM=SGM*CORR
    STM=STM*CORR
    J=0
    RETURN
4   MM=MM-1
5   RETURN
end subroutine avert
!
! ----------- Fadeeva function ---------------------------------------
! - If Y is less than 0.001 call the more efficient PFCN 
subroutine pfcn(X,Y,U,V,K)

    implicit none

    real(4), intent(in) :: X,Y
    real(4), intent(inout) :: U,V
    integer, intent(inout) :: K

    real(4) :: AM,CXP2,CXY2,D,DM,DP,EMN,EYP2,EYX,F,P,PI,PI2,SU,SV,SXP2,SXY2, &
               XNM,XNP,XP2,XY,XY2,Y2
    integer :: I,IS,N

    DATA PI,PI2,N,IS/3.14159 ,6.28318 ,5,0/
    DIMENSION EMN(10)
    IF(Y.GT.0.01)GO TO 40
    IF(IS.NE.0)GO TO 5
    IS=1
    DO 2 I=2,N
2   EMN(I)=EXP(-FLOAT(I-1)**2)/PI
5   CONTINUE
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
10  CONTINUE
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
20  U=SU
    V=SV
    RETURN
30  F=1.77245  *Y
    U=1./((1.+XY**2)*F)
    V=XY*U
    RETURN
!
40  CALL PFCNP(X,Y,U,V,K)
!
    RETURN
end subroutine pfcn
!
! ----------- Fadeeva function ---------------------------------------
! - If Y is less than 0.001 this is more efficient than PFCN
subroutine pfcnp(X,Y,U,V,L)

    implicit none

    real(4), intent(in) :: X,Y
    real(4), intent(inout) :: U,V
    integer, intent(inout) :: L

    real(4) :: c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18, &
               c19,c20,c21,CO,T,T1,W283,W287,Z
    
    integer :: I,II,J,K,M
!
!   PFCN YIELDS REAL AND IMAGINARY PART OF THE COMPLEX
!   PROBABILITY INTEGRAL
!
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
8   IF(C6.GE.0.0)GO TO 14
    ASSIGN 245 TO I
    GO TO 20
11  ASSIGN 257 TO I
    GO TO 46
14  ASSIGN 255 TO I
    GO TO 46
20  Z=C6*C6-C5*C5
    CO=EXP(Z)
    C7=CO+CO
    CO=C5*C6
    C9=CO+CO
    C8=-C7*SIN(C9)
    C7=C7*COS(C9)
46  C5=ABS(C5)
    C6=ABS(C6)
    IF(C5.GE.6.0)GO TO 219
50  IF(C6.LE.0.5)GO TO 65
    IF(C6.LE.3.0)GO TO 61
    IF(C6.GT.6.0)GO TO 219
    C9=0.5
    GO TO 73
61  IF(C6.LE.1.5)GO TO 71
    C9=0.25
    GO TO 73
65  C10=C6
    C6=0.5
    ASSIGN 128 TO J
71  C9=0.09375
73  C11=0.0
    C17=0.0
    C18=0.0
    ASSIGN 123 TO K
79  C21=C5-C11
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
end subroutine pfcnp
! ---------------------------------------------------------------------
!
subroutine wigner(DAV,D)
    implicit none
!             WIGNER YIELDS LEVEL SPACINGS WITH THE WIGNER DISTRIBUTION.
!             DAV IS THE AVERAGE LEVEL SPACING.
    real(4), intent(in) :: DAV
    real(4), intent(out) :: D
    real(4)::www
    real(4)::arbitrary

1   WWW=RANDOM(arbitrary)
    IF(WWW.EQ.0.) GO TO 1
! ---- 2/sqrt(pi) = 1.12837
    D=1.12837  *DAV*SQRT(-ALOG(WWW))
    RETURN
end subroutine wigner
!
subroutine porter(GNAV,GN)
    implicit none
!             PORTER YIELDS NEUTRON WIDTHS WITH THE PORTER-THOMAS
!             DISTRIBUTION. GNAV IS THE AVERAGE WIDTH.
!             V. NEUMANN SCHEME
    real(4), intent(in) :: GnAv
    real(4), intent(out) :: Gn
    real(4) :: arbitrary,r,rr,x

1   R=RANDOM(arbitrary)*0.929
    RR=1./(1.-R)**2
    X=R**2*RR
    IF(X.GT.170.) GO TO 1
    IF(RANDOM(arbitrary).GT.0.560*EXP(-X)*RR)GO TO 1
    GN=2.*GNAV*X
    RETURN
end subroutine porter
!
subroutine space(DAV,D)
    implicit none
    real(4), intent(in) :: dav
    real(4), intent(out) :: d
    real(4)::arbitrary,r,x

!             SPACE YIELDS LEVEL SPACINGS WITH A D-WEIGHTED WIGNER
!             DISTRIBUTION. DAV IS THE AVERAGE OF THE WIGNER DISTRIBU-
!             TION. V. NEUMANN SCHEME
! ---------------------------------------------------------------------
! ---- Sample a simpler function that envelopes the porter distribution.
! ---- If the sample is within the Porter, keep. If not re-sample
! ---------------------------------------------------------------------
1   R=RANDOM(arbitrary)
    X=(R/(1.-R))**0.3333
    IF(X.GT.13.) GO TO 1
    IF(RANDOM(arbitrary).GT.0.494*EXP(-X**2)/(1.-R)**2)GO TO 1
! ---- 2/sqrt(pi) = 1.12837
    D=1.12837  *DAV*X
    RETURN
end subroutine space

function random(arbitrary)
    implicit none
    real(4), intent(in) :: arbitrary
    real(4) :: random

    CALL random_number(random) ! built-in Fortran
    RETURN
end function random
!
subroutine simp(Y,M,N,SUM)
    implicit none
!
!         SIMPSON'S RULE
!
    real(4), intent(in) :: y
    real(4), intent(out) :: sum
    integer(4), intent(in) :: m,n
    integer(4) :: k,k1,k2
    DIMENSION Y(21)
    SUM=0.
    K1=M+1
    K2=N-1
    do 1 K=K1,K2,2
1       SUM=SUM+Y(K-1)+4.*Y(K)+Y(K+1)
        SUM=SUM/3.
    RETURN
end subroutine
!
function cbrt(MANT)
    implicit none

    integer(4), intent(in) :: mant
    real(4) :: expo,cbrt

    EXPO=1./3.
    CBRT=MANT**EXPO
    RETURN
end function cbrt
!
!
!
end module mc_routines
