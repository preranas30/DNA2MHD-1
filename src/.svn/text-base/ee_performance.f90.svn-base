!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 29/12/2012                                                                !!
!!                              ee_performance.f90                           !!
!!                                                                           !!
!!  performance                                                              !!
!!  -- performance_compute                                                   !!
!!  -- performance_linear                                                    !!
!!  -- performance_nl                                                        !!
!!  -- performance_rhs                                                       !!
!!  -- performance_par                                                       !!
!!  -- start_clock                                                           !!
!!  -- end_clock                                                             !!
!!                                                                     1.000 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                               performance                                 !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Tools for determining the numerical performance of different subroutines.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE performance
  USE par_mod
  USE mpi 
  USE linear_rhs
  USE nonlinearity
  USE time_advance
  USE field_solver, only: get_phi
  IMPLICIT NONE
  
  PUBLIC :: start_clock, end_clock, performance_compute
  PRIVATE

  INTEGER :: perf_handle


  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                             performance_compute                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE performance_compute
    IMPLICIT NONE

    IF(perf_test_lin) THEN
      IF(verbose) WRITE(*,*) "Calling performance_linear.",mype
      CALL performance_linear
      IF(verbose) WRITE(*,*) "Done with performance_linear.",mype
    END IF
    IF(perf_test_nl) THEN
      IF(mype==0.and.verbose) WRITE(*,*) "Calling performance_nl.",mype
      CALL performance_nl
      IF(mype==0.and.verbose) WRITE(*,*) "Done with performance_nl.",mype
    END IF
    IF(perf_test_rhs) THEN
      IF(mype==0.and.verbose) WRITE(*,*) "Calling performance_rhs.",mype
      CALL performance_rhs
      IF(mype==0.and.verbose) WRITE(*,*) "Done with performance_rhs.",mype
    END IF
    IF(perf_test_par) THEN
      IF(mype==0.and.verbose) WRITE(*,*) "Calling performance_par.",mype
      CALL performance_par
      IF(mype==0.and.verbose) WRITE(*,*) "Done with performance_par.",mype
    END IF

  END SUBROUTINE performance_compute
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                             performance_linear                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE performance_linear
    IMPLICIT NONE

    INTEGER :: handle,i
    COMPLEX :: rhs_out(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
    COMPLEX :: rhs_out2(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
    REAL :: lin_time
    CHARACTER(len=40) :: ic
    INTEGER :: num
    INTEGER :: rhs_lin_version_bak
    REAL :: norm,diff

    IF(np_hank.gt.1) STOP "gather_g1 not yet implemented for np_hank.gt.1"
    IF(np_spec.gt.1) STOP "gather_g1 not yet implemented for np_spec.gt.1"
    IF(np_kz.gt.1) STOP "gather_g1 not yet implemented for np_kz.gt.1"

    num=10
    !IF(verbose) WRITE(*,*) "Start of performance linear.",mype

    ic='DEFAULT'
    
    CALL get_io_number
    handle=io_number
    IF(mype==0) OPEN(unit=handle,file=trim(diagdir)//'/performance_linear.dat',status='unknown') 

    !IF(verbose) WRITE(*,*) 87,mype
    CALL initial_condition(ic)
   
    !IF(verbose) WRITE(*,*) 90,mype
    CALL get_phi(g_1)

    !IF(mype==0) WRITE(*,*) "Testing accuracy of rhs_lin_version selection vs. version 1."
    !CALL get_rhs_lin(g_1,phi,rhs_out,0)
    !rhs_lin_version_bak=rhs_lin_version
    !rhs_lin_version=1
    !CALL get_rhs_lin(g_1,phi,rhs_out2,0)
    !norm=REAL( sum(sum(sum(sum(abs(rhs_out2),1),1),1),1) )
    !diff=REAL( sum(sum(sum(sum(abs(rhs_out1-rhs_out2),1),1),1),1) )
    !IF(norm.gt.1.0e-15) THEN
!
!    ELSE
!
!    END IF

    !IF(verbose) WRITE(*,*) 95,mype
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "Starting loop for linear performance test."

    CALL start_clock
    DO i=1,num
      CALL get_rhs_lin(g_1,phi,rhs_out,0)
    END DO
    CALL end_clock(lin_time)

    IF(mype==0) WRITE(*,*) "Done."
    IF(mype==0) WRITE(handle,*) "#",num, " calls to get_rhs_lin:"
    IF(mype==0) WRITE(handle,*) lin_time
    IF(mype==0) WRITE(*,*) num, " calls:",lin_time
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

    IF(mype==0) CLOSE(handle)

  END SUBROUTINE performance_linear


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                performance_nl                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE performance_nl
    IMPLICIT NONE

    INTEGER :: handle,i
    COMPLEX :: rhs_out(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
    REAL :: nl_time

    CHARACTER(len=40) :: ic
    INTEGER :: num

    IF(np_hank.gt.1) STOP "performance_nl not yet implemented for np_hank.gt.1"
    IF(np_spec.gt.1) STOP "performance_nl not yet implemented for np_spec.gt.1"
    IF(np_kz.gt.1) STOP "performance_nl not yet implemented for np_kz.gt.1"
    ic='DEFAULT'
    num=10
    
    CALL get_io_number
    handle=io_number
    IF(mype==0) OPEN(unit=handle,file=trim(diagdir)//'/performance_nonlinear.dat',status='unknown') 

    CALL initial_condition(ic)
    CALL initialize_fourier
   
    CALL get_phi(g_1)

    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "Starting loop for nonlinear performance test."

    CALL start_clock
    DO i=1,num
      CALL get_rhs_nl(g_1,phi,rhs_out)
    END DO
    CALL end_clock(nl_time)

    IF(mype==0) WRITE(*,*) "Done."
    IF(mype==0) WRITE(handle,*) "#",num, " calls to get_rhs_nl:"
    IF(mype==0) WRITE(handle,*) nl_time
    IF(mype==0) WRITE(*,*) num," calls:",nl_time
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

    CLOSE(handle)

  END SUBROUTINE performance_nl


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                               performance_rhs                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE performance_rhs
    IMPLICIT NONE

    INTEGER :: handle,i
    COMPLEX :: rhs_out(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
    REAL :: nl_time

    CHARACTER(len=40) :: ic
    INTEGER :: num

    IF(np_hank.gt.1) STOP "performance_rhs not yet implemented for np_hank.gt.1"
    IF(np_spec.gt.1) STOP "performance_rhs not yet implemented for np_spec.gt.1"
    IF(np_kz.gt.1) STOP "performance_rhs not yet implemented for np_kz.gt.1"

    ic='DEFAULT'
    num=10
    
    CALL get_io_number
    handle=io_number
    IF(mype==0) OPEN(unit=handle,file=trim(diagdir)//'/performance_rhs.dat',status='unknown') 

    CALL initial_condition(ic)
    CALL initialize_fourier

    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "Starting loop for rhs performance test."

    itime=0
    time=0.0
    max_itime=num
    max_time=1000.0
    istep_ffm=0
    istep_energy3d=0
    istep_energy=0
    istep_fmom3d=0
    istep_gamma=0
    istep_nltest=0
    istep_schpt=0
    adapt_dt_nl=.false.

    CALL start_clock
    CALL iv_solver
    CALL end_clock(nl_time)

    IF(mype==0) WRITE(*,*) "Done."
    IF(mype==0) WRITE(handle,*) "#",num, " iterations for iv solver:"
    IF(mype==0) WRITE(handle,*) nl_time
    IF(mype==0) WRITE(*,*) num," iterations for iv_solver:",nl_time
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

    CLOSE(handle)
    STOP

  END SUBROUTINE performance_rhs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                               performance_par                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE performance_par
    IMPLICIT NONE

    INTEGER :: handle,i
    COMPLEX :: g_bounds(0:nkx0-1,0:nky0-1,lkz1:lkz2,2)
    REAL :: nl_time

    CHARACTER(len=40) :: ic
    INTEGER :: num,ierr

    IF(np_hank.gt.1) STOP "performance_par not yet implemented for np_hank.gt.1"
    IF(np_spec.gt.1) STOP "performance_par not yet implemented for np_spec.gt.1"
    IF(np_kz.gt.1) STOP "performance_par not yet implemented for np_kz.gt.1"

    ic='DEFAULT'
    num=20000
    
    CALL get_io_number
    handle=io_number
    IF(mype==0) OPEN(unit=handle,file=trim(diagdir)//'/performance_par.dat',status='unknown') 

    CALL initial_condition(ic)

    !CALL initialize_fourier

    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "Starting loop for parallel performance test."

    CALL start_clock
    DO i=1,num
      CALL get_v_boundaries2(g_1,g_bounds)
    END DO
    CALL end_clock(nl_time)

    IF(mype==0) WRITE(*,*) "Done."
    IF(mype==0) WRITE(handle,*) "#number of processors:", np_herm
    IF(mype==0) WRITE(handle,*) "#nv0:", nv0
    IF(mype==0) WRITE(handle,*) "#",num, " iterations for iv solver:"
    IF(mype==0) WRITE(handle,*) nl_time
    IF(mype==0) WRITE(*,*) num," iterations for iv_solver:",nl_time
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

    CLOSE(handle)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    STOP

  END SUBROUTINE performance_par

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                               start_clock                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE start_clock
    IMPLICIT NONE

    INTEGER :: MCLOCK
    
    perf_monitor(1)=MCLOCK()

  END SUBROUTINE start_clock


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                end_clock                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE end_clock(seconds)
    IMPLICIT NONE

    REAL, INTENT(out) :: seconds
    INTEGER :: max_time,ierr
    INTEGER :: MCLOCK
    INTEGER :: diff

    perf_monitor(2)=MCLOCK()  
    diff=perf_monitor(2)-perf_monitor(1)
    CALL MPI_ALLREDUCE(diff,max_time,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    seconds=REAL(max_time)/1000.0

  END SUBROUTINE end_clock
  

END MODULE performance

