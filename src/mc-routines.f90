module mc_routines
    use physics
    contains

subroutine musc(AP,UM,A,abundance,BE,PE,eff_temp,NL,AA,chan_rad,GG, &
                R,EFF,J2X,J2N,XN,E,SG,SP,G,GNR,GNRIN,DJL,SSCF,DSSCF,N,K, &
                I,L,J,NI,LX,RN,RX,ZH,PO,PS,DPO,DPS,numResPairs,SGM,STM,DSGM,DSTM, &
                TMC,DTMC,SNC,ITYPE,IHIST,LJX,LJN,JN2,JMX, RB)

    implicit none
!             MUSC  YIELDS THE MULTIPLE SCATTERING CORRECTION FOR A
!             CYLINDRICAL CAPTURE SAMPLE (MONTE CARLO)
!                XN(N)      CAPTURE SAMPLE THICKNESS
!                RX(N)      CAPTURE SAMPLE RADIUS
!                RB(N)      NEUTRON BEAM RADIUS
!                RN(N)      FILTER  SAMPLE THICKNESS
    
    real(8), allocatable, dimension(:) :: AP,UM,A,SS,STT,abundance,BE,PE,eff_temp, &
                              AA,chan_rad,XN,E,SG,SP,RN,RX,ZH,RB,ZST,ZSG,ZT,ZF0, &
                              GNL,PSM,POM,EE,STMC,SGMC,TM
    integer, allocatable, dimension(:) :: NL,JN2,JMX
    real(8), allocatable, dimension(:,:) :: GG,R,EFF,DOP
    integer, allocatable, dimension(:,:) :: J2X,J2N
    real(8), allocatable, dimension(:,:,:) :: G,GNR,GNRIN,DJL,SINE,COSE,PLE,GGE
    integer, allocatable, dimension(:,:,:) :: LJX,LJN
    real(8), allocatable, dimension(:,:,:,:) :: DE,GNE,GNINE

    real(8) :: AK,B1,B11,B12,B2,B23,B3,BBOTH,C0,CG,CT,CTH,CTHC,DD,  &
               DEN,DPO,DPS,DSGM,DSSCF,DSTM,DTMC,DUMMY,DYN,EDD,EDGG,ETA,  &
               EXN1,FE,FP,GNINS,GNS,GT,H,HNEG,HPOS,O,PHI,PO,POMG,PS,PSMG,RHO, &
               RK,SAM,SC,SCC,SEM,SGM,SIM,SN,SNC,SOM,SQ,SSCF,ST,STH,STM,SUM,   &
               SUMF0,SUMSG,SUMST,SUMT,SYM,T,TF,TMC,TP,TX,U,UU,V,VL,VV,W,  &
               WG,WI,WN,X,XI,Y,Z,ZSIR
    integer :: I,IHIST,ITYPE,IZ,J,JL,JX,K,KZ,L,L1,L2,LP,LX,M,MCHD,MDIV,MM,MP,N, &
               NC,NH,NI,NP,NS,numResPairs

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
            DOP(NC,I)=SQRT(0.72531*A(I)/eff_temp(I)/EE(NC))
    !
             CALL ENDEP(EE(NC),A(I),BE(I),PE(I),AP(I),UM(I),EDGG,EDD)
            
            LX=NL(I)
            do L=1,LX
                GGE(NC,L,I)=GG(L,I)*EDGG
                AK=chan_rad(I)*SQRT(EE(NC)/AA(I))/143.92
        !
        ! ---- penetrabilities, etc..... --------------------------------------
                 CALL PEPS(AK,L,DUMMY,VL)
        !
                RK=AK*R(L,I)/chan_rad(I)
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
    NH=INT(ZH(K))
    NP=numResPairs
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
        X=RB(N)*SQRT(random())         ! SAMPLE ALONG THE BEAM RADIUS
        Y=RB(N)*SQRT(random())
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
                    C0=4.617E3/EE(NC)*AA(I)*abundance(I)*G(JL,L1,I)
!                           MP: COUNTER FOR RESONANCE PAIRS
                    MP=1
!                       FIND CENTRAL LEVEL INTERVAL AND ENERGY IN IT
!
                    CALL SPACE(DE(NC,J,L1,I),DD)
!
                    H=DD*random()
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
        IZ=int( 20.*ST/STT(K)+1. )
        IF(IZ.LT.1)  IZ=1
        IF(IZ.GT.101)IZ=101
        ZST(IZ)=ZST(IZ)+1.
        IZ=int( 20.*SC/SG(K)+1. )
        IF(IZ.LT.1)  IZ=1
        IF(IZ.GT.101)IZ=101
        ZSG(IZ)=ZSG(IZ)+1.
        IZ=int( 100.*T+1. )
        IF(IZ.LT.1)  IZ=1
        IF(IZ.GT.101)IZ=101
        ZT(IZ)=ZT(IZ)+1.
        IZ=int( 20.*(WG/ST)/(XN(N)*SG(K)*EXP(-RN(N)*STT(K))) )
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
7       IF(random().GT.WN)GO TO 19
!             FREE PATH SAMPLING
        FP=-LOG(1.-WI*random())/ST
!             COORDINATES OF COLLISION
        X=X+U*FP
        Y=Y+V*FP
        Z=Z+W*FP
!               SCATTERING ANGLES
        PHI=6.28318 *random()
        CTHC=2.*random()-1.
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
subroutine muss(AP,UM,A,abundance,BE,PE,eff_temp,NL,AA,chan_rad,GG, &
                R,EFF,J2X,J2N,E,SG,SP,G,GNR,GNRIN,DJL,SSCF,DSSCF,N,K, &
                I,L,J,NI,LX,RN,RX,ZH,PO,PS,DPO,DPS,numResPairs,SGM,STM,DSGM,DSTM, &
                TMC,DTMC,SNC,XNSS,ITYPE,IHIST,LJX,LJN,JN2,JMX)

    implicit none
!             MUSS YIELDS THE MULTIPLE SCATTERING CORRECTION FOR A
!             SPHERICAL SHELL (MONTE CARLO)

    real(8), allocatable, dimension(:) :: AP,UM,A,SS,STT,abundance,BE,PE,eff_temp, &
                              AA,chan_rad,E,SG,SP,RN,RX,ZH,ZST,ZSG,ZT,ZF0, &
                              GNL,EE,STMC,SGMC,TM,PSM,POM
    integer, allocatable, dimension(:) :: NL,JN2,JMX
    real(8), allocatable, dimension(:,:) :: GG,R,EFF,DOP
    integer, allocatable, dimension(:,:) :: J2X,J2N
    real(8), allocatable, dimension(:,:,:) :: G,GNR,GNRIN,DJL,SINE,COSE,PLE,GGE
    integer, allocatable, dimension(:,:,:) :: LJX,LJN
    real(8), allocatable, dimension(:,:,:,:) :: DE,GNE,GNINE

    real(8) :: A1,A2,B11,B12,B23,BBOTH,C0,CG,CT,CTH, &
               DD,DEN,DPO,DPS,DSGM,DSSCF,DSTM,DTMC,DYN,EDD,EDGG,ETA, &
               FE,GNINS,GNS,GT,H,HNEG,HPOS,PO,POMG,PS,PSMG,   &
               RK,SAM,SC,SCC,SEM,SGM,SIM,SN,SNC,SOM,SQ,SSCF,ST,STH,STM,SUM,     &
               SUMF0,SUMSG,SUMST,SUMT,SYM,T,TMC,TX,UU,VL,VV,W,WG, &
               WI,WN,XI,XNSS,Z,Q
    integer :: I,IHIST,ITYPE,IZ,J,JL,JX,K,KZ,L,L1,L2,LP,LX,M,MCHD,MDIV,MM,MP,N, &
               NC,NH,NI,NP,NS,numResPairs

    allocate(SS(100),STT(100),ZST(101),ZSG(101),ZT(101),ZF0(101),GNL(8),        &
             COSE(10,4,10),SINE(10,4,10),EE(10),STMC(1000),SGMC(1000),TM(1000), &
             PSM(1000),POM(1000),PLE(10,4,10),GGE(10,4,10),DOP(10,10),          &
             DE(10,8,4,10),GNE(10,8,4,10),GNINE(10,8,4,10))
    
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
    do NC=1,10
        IF(NC.EQ.1)EE(NC)=E(K)
        IF(NC.GT.1)EE(NC)=EE(NC-1)*FE
        do I=1,NI
            DOP(NC,I)=SQRT(0.72531*A(I)/eff_temp(I)/EE(NC))
            LX=NL(I)
            do L=1,LX
        !
                CALL ENDEP(EE(NC),A(I),BE(I),PE(I),AP(I),UM(I),EDGG,EDD)
        !
                GGE(NC,L,I)=GG(L,I)*EDGG
                RK=chan_rad(I)*SQRT(EE(NC)/AA(I))/143.92
        !
                CALL PEPS(RK,L,PLE(NC,L,I),VL)
        !
                PLE(NC,L,I)=PLE(NC,L,I)*R(L,I)/chan_rad(I)
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
    NH=INT(ZH(K))
    NP=numResPairs
    do M=1,NH
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
                    C0=4.617E3/EE(NC)*AA(I)*abundance(I)*G(JL,L1,I)
                    MP=1
!                    FIND CENTRAL LEVEL INTERVAL AND ENERGY IN IT
!
                    CALL SPACE(DE(NC,J,L1,I),DD)
!
                    H=DD*random()
17                  GNS=0.
                    GNINS=0.
                    do L=L1,L2,2
!                         FIND REDUCED NEUTRON WIDTH
!
                        CALL PORTER(GNE(NC,J,L,I),GNL(L))
!
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
!                     CAPTURE CROSS SECTION
                    SCC=CG*UU
                    SN=SN+SCC
!                     EFFECTIVE CAPTURE CROSS SECTION
                    SC=SC+SCC*EFF(L1,I)
!                     TOTAL CROSS SECTION
                    do L=L1,L2,2
                        ST=ST+CT*GNL(L)*(UU*COSE(NC,L,I)+VV*SINE(NC,L,I))
                    end do
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
        IZ=int( 20.*ST/STT(K)+1. )
        IF(IZ.LT.1)  IZ=1
        IF(IZ.GT.101)IZ=101
        ZST(IZ)=ZST(IZ)+1.
        IZ=int( 20.*SC/SG(K)+1. )
        IF(IZ.LT.1)  IZ=1
        IF(IZ.GT.101)IZ=101
        ZSG(IZ)=ZSG(IZ)+1.
        IZ=int( 100.*T+1. )
        IF(IZ.LT.1)  IZ=1
        IF(IZ.GT.101)IZ=101
        ZT(IZ)=ZT(IZ)+1.
        IZ=int( 20.*(WG/ST)/(EXP(-XNSS*STT(K))*XNSS*SG(K))+1. )
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
7       IF(random().GT.WN)GO TO 19
!             CALCULATE RADIAL COORDINATE (Z) OF NC-TH COLLISION
        T=-LOG(1.-WI*random())/ST
        IF(T.GT.(Z*CTH-Q))  T=T+2.*Q
        Z=SQRT(Z**2+  T**2-2.*Z*  T*CTH)
!             ANGLE WITH RESPECT TO RADIUS AFTER COLLISION
        CTH=2.*random()-1.
        STH=SQRT(1.-CTH**2)
!             MAXIMUM PATH LENGTH IN MATTER
        Q=0.
        IF(CTH.GT.SQRT(1.-(RN(N)/Z)**2))Q=SQRT(RN(N)**2-(Z*STH)**2)
        TX=Z*CTH+SQRT(RX(N)**2-(Z*STH)**2)-2.*Q
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
        SGM   =SGM   +SGMC(M)
        STM   =STM   +STMC(M)
    end do
    SGM   =SGM   /FLOAT(NS)
    TMC=TMC/FLOAT(NS)
    STM   =STM   /FLOAT(NS)
    PO=SOM/ZH(K)
    PS=SUM/ZH(K)
    SNC=SNC/ZH(K)
!             STATISTICAL UNCERTAINTIES
    do M=1,NS
        DTMC=DTMC+(TM(M)-TMC)**2
        SIM=SIM+(SGMC(M)-SGM   )**2
        SYM=SYM+(STMC(M)-STM   )**2
    end do
    B11=0.
    B12=0.
    B23=0.
    BBOTH=0.
    do M=1,MM
        BBOTH=BBOTH+(POM(M)+PSM(M)-PO-PS)*(SGMC(M)-SGM)
                      B11=B11+(POM(M)-PO)*(PSM(M)-PS)
        SAM=SAM+(POM(M)-PO     )**2
        SEM=SEM+(PSM(M)-PS     )**2
    end do
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
subroutine moct(AP,UM,A,abundance,BE,PE,eff_temp,NL,AA,chan_rad,GG,R,EFF,J2X,J2N,XN,E,SG,SP,G,   &
                GNR,GNRIN,DJL,N,K,I,L,J,NI,LX,ZH,numResPairs,SGM,STM,DSGM,DSTM, &
                TMC,DTMC,TAU,DTAU,IHIST,LJX,LJN,JN2,JMX)

    implicit none
!             MOCT  YIELDS THE SELF-SHIELDING CORRECTION FACTOR AND
!             TRANSMISSION FOR A CYLINDRICAL SAMPLE

    real(8), allocatable, dimension(:) :: AP,UM,A,STT,abundance,BE,PE,eff_temp,AA,chan_rad,   &
                                          XN,E,SG,SP,ZH,ZST,ZSG,ZT,GNL,TMCG,TN, &
                                          SGMC,STMC,DOP
    integer, allocatable, dimension(:) :: NL,JN2,JMX
    real(8), allocatable, dimension(:,:) :: GG,R,EFF,PL,SINE,COSE,GGE,VL,EDGG,EDD
    integer, allocatable, dimension(:,:) :: J2X,J2N
    real(8), allocatable, dimension(:,:,:) :: G,GNR,GNRIN,DJL,GN,GNIN,DE
    integer, allocatable, dimension(:,:,:) :: LJX,LJN

    real(8) :: AK,AVT,AVTS,C0,CG,CT,DD,DEN,DSGM,DSTM,DT,DTS,DTS2,DTAU,DTMC,    &
               DUMMY,ETA,GNINS,GNS,GT,H,HNEG,HPOS,RK,SC,SCC,SGM,SIM,SN,ST,STM, &
               SUMSG,SUMST,SUMT,SYM,T,TAU,TMC,UU,VV,XI,TNUM,TS,TSSM,TSUM
    integer :: I,IHIST,IT,IZ,J,JL,JX,K,KZ,L,L1,L2,LP,LX,M,MCHD,MDIV,MM,MP,N,   &
               NH,NHIST,NI,NP,numResPairs

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
    do IT=2,101
        TN(IT)=TN(IT-1)+1.
        ZST(IT)=0.
        ZSG(IT)=0.
        ZT(IT)=0.
    end do
    do I=1,NI
        DOP(I)=SQRT(0.72531*A(I)/eff_temp(I)/E(K))
!
        CALL ENDEP(E(K),A(I),BE(I),PE(I),AP(I),UM(I),EDGG(I,K),EDD(I,K))
!
        LX=NL(I)
        do L=1,LX
            GGE(L,I)=GG(L,I)*EDGG(I,K)
            AK=chan_rad(I)*SQRT(E(K)/AA(I))/143.92
!
            CALL PEPS(AK,L,DUMMY,VL(L,I))
!
            RK=AK*R(L,I)/chan_rad(I)
!
            CALL PEPS(RK,L,PL(L,I),DUMMY)
!
            SINE(L,I)=SIN(2.*PL(L,I))
            COSE(L,I)=COS(2.*PL(L,I))
            JX=(J2X(L,I)-J2N(L,I))/2+1
            do JL=1,JX
                J=(J2N(L,I)-JN2(I))/2+JL
                DE(J,L,I)=DJL(JL,L,I)*EDD(I,K)
                GN  (J,L,I)=GNR  (JL,L,I)*SQRT(E(K)*1000.)*VL(L,I)*EDD(I,K)
                GNIN(J,L,I)=GNRIN(JL,L,I)*SQRT(E(K)*1000.)*VL(L,I)*EDD(I,K)
            end do
        end do
    end do
!             BEGIN LOOP OF NEUTRON HISTORIES
    NH=INT(ZH(K))
    NP=numResPairs
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
    do I=1,NI
        JX=JMX(I)
        do J=1,JX
            DO 93 LP=1,2
                L1=LJN(LP,J,I)
                L2=LJX(LP,J,I)
                IF(L1.LT.0)GO TO 93
            !           4.617E3=2.605E3*1.772454
                JL=J-(J2N(L1,I)-JN2(I))/2
                C0=4.617E3/E(K)*AA(I)*abundance(I)*G(JL,L1,I)
                MP=1
            !            FIND CENTRAL LEVEL INTERVAL AND ENERGY IN IT
            !
                    CALL SPACE(DE(J,L1,I),DD)
            !
                H=DD*random()
17              GNS=0.
                GNINS=0.
                do L=L1,L2,2
            !             FIND REDUCED NEUTRON WIDTH
            !
                    CALL PORTER(GN(J,L,I),GNL(L))
            !
                    GNS=GNS+GNL(L)
                    GNINS=GNINS+GNIN(J,L,I)
                end do
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
                do L=L1,L2,2
                    ST=ST+CT*GNL(L)*(UU*COSE(L,I)+VV*SINE(L,I))
                end do
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
18              IF(H.LT.0.)GO TO 16
                HPOS=H
                HNEG=H-DD
                H=HNEG
                GO TO 17
            !             INCLUDE ANOTHER PAIR OF RESONANCES IF REQUIRED
16              IF(MP.EQ.NP)GO TO 93
                MP=MP+1
            !
                    CALL WIGNER(DE(J,L1,I),DD)
            !
                HPOS=HPOS+DD
                H=HPOS
                GO TO 17
93          continue
        end do
    end do
    SN=ST-SN
    T=EXP(-XN(N)*ST)
    TS=T*ST
    write(99,'(2E12.5)')ST, T
    IF(IHIST.NE.1)GO TO 52
    IZ=int( 20.*ST/STT(K)+1. )
    IF(IZ.LT.1)  IZ=1
    IF(IZ.GT.101)IZ=101
    ZST(IZ)=ZST(IZ)+1.
    IZ=int( 20.*SC/SG(K)+1. )
    IF(IZ.LT.1)  IZ=1
    IF(IZ.GT.101)IZ=101
    ZSG(IZ)=ZSG(IZ)+1.
52  IZ=int( 2.*T/DT+1. )
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
    do M=1,MM
        SGM=SGM+SGMC(M)
        STM=STM+STMC(M)
        TMC=TMC+TMCG(M )
    end do
    TMC=TMC/FLOAT(MM)
    SGM=SGM/FLOAT(MM)
    STM=STM/FLOAT(MM)
!           STATISTICAL UNCERTAINTY
    DTMC=0.
    DTS=0.
    do M=1,MM
        SIM=SIM+(SGMC(M)-SGM)**2
        SYM=SYM+(STMC(M)-STM)**2
        DTS=DTS+(STMC(M)-STM)*(TMCG(M)-TMC)
        DTMC=DTMC+(TMCG(M)-TMC)**2
    end do
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
    do IZ=1,101
        SUMST=SUMST+ZST(IZ)
        SUMSG=SUMSG+ZSG(IZ)
        SUMT =SUMT + ZT(IZ)
    end do
    do IZ=1,101
        ZST(IZ)=ZST(IZ)/SUMST
        ZSG(IZ)=ZSG(IZ)/SUMSG
        ZT(IZ)=ZT(IZ)/SUMT
    end do
    WRITE(8,199)(IZ,ZST(IZ),ZSG(IZ),ZT(IZ),IZ=1,101)
199 FORMAT(1H1,///,' NO.    TALLY FOR HISTOGRAM OF                ST     SG     T           ',(I4,0PF8.4,2F7.4))
!
!C     CALL PLOT(TN,ZT,100,3,3,1,1,1,1,100.,0.,0.125,100.,0.,0.125,'DISTR
!C   1IBUTION OF TRANSMISSION VALUES',1)
53  CONTINUE
    RETURN
end subroutine moct
!
end module mc_routines
