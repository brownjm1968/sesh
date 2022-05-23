program test_sesh

    use physics

    implicit none

    logical :: pass=.false.

    pass = test_peps()
    if( pass ) then 
        print *, "passed"
    else
        print *, "failed"
    end if

contains

function test_peps() result(pass)
    use physics
    logical :: pass,pl_good,vl_good
    integer :: L ! ang. momentum
    integer :: i
    real(8) :: RK,PL,VL
    real(8), dimension(4) :: pl_answers = (/2.0d0,0.d0,0.d0,0.d0/)
    real(8), dimension(4) :: vl_answers = (/1.0d0,0.d0,0.d0,0.d0/)

    RK = 2.0d0
    pass = .false.
    pl_good = .false.
    vl_good = .false.

    ! -- PL answers
    pl_answers(2) = 2.0d0 - atan(RK)
    pl_answers(3) = 2.0d0 - atan(3.d0*RK/(3-RK**2))
    pl_answers(4) = 2.0d0 - atan(RK*(15.d0-RK**2)/(15.d0-6.d0*RK**2))
    ! -- VL answers
    vl_answers(2) = RK**2/(1.d0+RK**2)
    vl_answers(3) = RK**4/(9.d0+3.d0*RK**2+RK**4)
    vl_answers(4) = RK**6/(225.d0+45.d0*RK**2+6.d0*RK**4+RK**6)


    do i=1,4
        call peps(RK,i,PL,VL)
        pl_good = (PL - pl_answers(i)) < 1e-30
        vl_good = (VL - pl_answers(i)) < 1e-30
        print *, "-------"
        print *, "l,PL,VL= ", i, PL, VL
        if( vl_good .and. pl_good ) then
            pass = .true.
        else
            pass = .false.
            return
        end if
    end do 
    
    return
end function test_peps

end program