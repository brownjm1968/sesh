module physics
    contains

!
subroutine peps(RK,L,PL,VL)
!             PEPS CALCULATES CENTRIFUGAL-BARRIER PENETRABILITIES, VL,
!             AND HARD-SPHERE PHASE SHIFTS, PL.
    
    implicit none

    real(8), intent(in) :: RK ! wave-number * radius = rho
    integer, intent(in) :: L
    real(8), intent(out) :: PL
    real(8) :: VL ! sammy penet. divided by rho

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

    real(8), intent(in) :: A,B,P,AP,UM
    real(8), intent(inout) :: E

    real(8) :: Y1,Y2,A1,E1,E2,ER,F1,F2,GR,SUM1,SUM2,U,U1X,U2X,UB,U1,U2,EDGG,EDD
    integer :: N

    DIMENSION Y1(21),Y2(21)
!    PRINT*, E, A, B, P, AP, UM, EDGG, EDD
    E=0.001*E ! keV to MeV
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

    real(8), intent(in) :: U,A,PE,AP,UM

    real(8) :: rho,VARJ2
    real(8) :: AI,AU4,T,UEFF

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

    real(8), intent(inout) :: SGMC, STMC
    real(8), intent(in) :: SC, ST
    integer, intent(inout) :: NS
    integer, intent(in) :: I

    real(8) :: corr,W
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
    RETURN
end subroutine aver
!
subroutine averp(PZ,PZG,PS,PSG,MM,I)

    implicit none

    real(8), intent(inout) :: PZ,PS
    real(8), intent(in) :: PZG,PSG
    integer, intent(inout) :: MM
    integer, intent(in) :: I

    real(8) :: W,corr
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
    RETURN
end subroutine averp
!
subroutine avert(TM,T,SGM,SC,STM,ST,MM,I)

    implicit none

    real(8), intent(in) :: T,SC,ST
    integer, intent(in) :: I
    real(8), intent(inout) :: TM
    integer, intent(inout) :: MM

    real(8) :: corr,W,SGM,STM
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
    RETURN
end subroutine avert
!
! ----------- Fadeeva function ---------------------------------------
! - If Y is less than 0.001 call the more efficient PFCN 
subroutine pfcn(X,Y,U,V,K)

    implicit none

    real(8), intent(in) :: X,Y
    real(8), intent(inout) :: U,V
    integer, intent(inout) :: K

    real(8), dimension(10) :: EMN
    real(8) :: AM,CXP2,CXY2,D,DM,DP,EYP2,EYX,F,P,PI,PI2,SU,SV,SXP2,SXY2, &
               XNM,XNP,XP2,XY,XY2,Y2
    integer :: I,IS,N

    PI = 3.14159d0
    PI2 = 6.28318d0
    N = 5
    IS = 0

    IF(Y.GT.0.01)GO TO 40
    IF(IS.NE.0)GO TO 5
    IS=1
    do I=2,N
        EMN(I)=EXP(-FLOAT(I-1)**2)/PI
    end do
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

    real(8), intent(in) :: X,Y
    real(8), intent(inout) :: U,V
    integer, intent(inout) :: L

    real(8) :: c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18, &
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
    real(8), intent(in) :: DAV
    real(8), intent(out) :: D
    real(8)::www

1   WWW=random()
    IF(WWW.EQ.0.) GO TO 1
! ---- 2/sqrt(pi) = 1.12837
    D=1.12837  *DAV*SQRT(-LOG(WWW))
    RETURN
end subroutine wigner
!
subroutine porter(GNAV,GN)
    implicit none
!             PORTER YIELDS NEUTRON WIDTHS WITH THE PORTER-THOMAS
!             DISTRIBUTION. GNAV IS THE AVERAGE WIDTH.
!             V. NEUMANN SCHEME
! ---------------------------------------------------------------------
! ---- Sample a simpler function that envelopes the porter distribution.
! ---- If the sample is within the Porter, keep. If not re-sample
! ---------------------------------------------------------------------
    real(8), intent(in) :: GnAv
    real(8), intent(out) :: Gn
    real(8) :: r,rr,x,norm

    norm = 0.560 ! normalize max PDF value to one

1   R=random()*0.929
    RR=1./(1.-R)**2
    X=R**2*RR
    IF(X.GT.170.) GO TO 1
    IF(random().GT.norm*EXP(-X)*RR)GO TO 1
    GN=2.*GNAV*X
    RETURN
end subroutine porter
!
subroutine space(DAV,D)
    implicit none
    real(8), intent(in) :: dav
    real(8), intent(out) :: d
    real(8)::r,x

!             SPACE YIELDS LEVEL SPACINGS WITH A D-WEIGHTED WIGNER
!             DISTRIBUTION. DAV IS THE AVERAGE OF THE WIGNER DISTRIBU-
!             TION. V. NEUMANN SCHEME
1   R=random()
    X=(R/(1.-R))**0.3333
    IF(X.GT.13.) GO TO 1
    IF(random().GT.0.494*EXP(-X**2)/(1.-R)**2)GO TO 1
! ---- 2/sqrt(pi) = 1.12837
    D=1.12837  *DAV*X
    RETURN
end subroutine space

function random()
    implicit none
    real(8) :: random

    CALL random_number(random) ! built-in Fortran
    RETURN
end function random
!
subroutine simp(Y,M,N,SUM)
    implicit none
!
!         SIMPSON'S RULE
!
    real(8), intent(in) :: y
    real(8), intent(out) :: sum
    integer(4), intent(in) :: m,n
    integer(4) :: k,k1,k2
    DIMENSION Y(21)
    SUM=0.
    K1=M+1
    K2=N-1
    do K=K1,K2,2
        SUM=SUM+Y(K-1)+4.*Y(K)+Y(K+1)
    end do
    SUM=SUM/3.
    RETURN
end subroutine simp
!
function cbrt(MANT)
    implicit none

    integer, intent(in) :: mant
    real(8) :: expo,cbrt

    EXPO=1./3.
    CBRT=MANT**EXPO
    RETURN
end function cbrt
!
!

end module physics
