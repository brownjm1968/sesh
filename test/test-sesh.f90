program test_sesh

    use physics

    implicit none

    logical :: pass=.false.

    pass = test_peps()
    pass = test_cbrt()
    if( pass ) then 
        print *, "passed"
    else
        print *, "failed"
    end if

contains

 function test_simp() result(pass)
    use physics
    implicit none
    logical :: pass
    real(8),dimension(21) :: function = 2.0d0
    real(8) :: sum_gold,sum

    function(1) = 3.0d0
    function(21) = 3.0d0
    ! simpsons rule but only summing heights, not area
    sum_gold = 7.d0/3.d0 + 19*2.d0 + 7.d0/3.d0

    call simp(function,1,21,sum)

    if( abs(sum-sum_gold)<1e-30 ) then
        pass = .true.
    else
        pass = .false.
        print *, "Failed test_simp."
        print *, "Gold, actual = ", sum_gold,sum
    end if
    
end function test_simp

function test_cbrt() result(pass)
    use physics
    implicit none
    logical :: pass
    real(8) :: cbrt_gold,radius,exponent
    integer :: A

    A = 181
    exponent = 1.0/3.0
    cbrt_gold = A**exponent

    radius = cbrt(A)
    
    if( abs(radius-cbrt_gold)<1e-30 ) then
        pass = .true.
    else
        pass = .false.
        print *, "Failed test_cbrt."
        print *, "Gold, actual = ", cbrt_gold,radius
    end if

    return
end function


function test_peps() result(pass)
    use physics
    implicit none
    logical :: pass,pl_good,vl_good
    integer :: L ! ang. momentum
    integer :: i
    real(8) :: RK,PL,VL
    real(8), dimension(4) :: pl_gold = (/2.0d0,0.d0,0.d0,0.d0/)
    real(8), dimension(4) :: vl_gold = (/1.0d0,0.d0,0.d0,0.d0/)

    RK = 2.0d0
    pass = .false.
    pl_good = .false.
    vl_good = .false.

    ! -- PL answers
    pl_gold(2) = 2.0d0 - atan(RK)
    pl_gold(3) = 2.0d0 - atan(3.d0*RK/(3-RK**2))
    pl_gold(4) = 2.0d0 - atan(RK*(15.d0-RK**2)/(15.d0-6.d0*RK**2))
    ! -- VL answers
    vl_gold(2) = RK**2/(1.d0+RK**2)
    vl_gold(3) = RK**4/(9.d0+3.d0*RK**2+RK**4)
    vl_gold(4) = RK**6/(225.d0+45.d0*RK**2+6.d0*RK**4+RK**6)


    do i=1,4
        call peps(RK,i,PL,VL)
        pl_good = abs(PL - pl_gold(i)) < 1e-30
        vl_good = abs(VL - vl_gold(i)) < 1e-30
        print *, "-------"
        print *, "l,PL,VL= ", i, PL, VL
        if( vl_good .and. pl_good ) then
            pass = .true.
        else
            pass = .false.
            print *, "Gold   = ", pl_gold(i), " & ", vl_gold(i)
            print *, "Actual = ", PL, " & ", VL
            print *, "Failed test_peps."
            return
        end if
    end do 
    
    return
end function test_peps

end program