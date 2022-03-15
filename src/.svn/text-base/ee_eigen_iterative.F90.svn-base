!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 29/12/2012                                                                !!
!!                              cc_eigen_iterative.f90                       !!
!!                                                                           !!
!!  eigen_iterative                                                          !!
!!  -- ev_iterative                                                          !!
!!  -- analyze_eigenpairs                                                    !!
!!  -- matmul_k                                                              !!
!!  -- calc_dt_ev                                                            !!
!!  -- eigenpairs_dt                                                         !!
!!  -- gather_g1                                                             !!
!!  -- sort_evs                                                              !!
!!  -- round_to_prec                                                         !!
!!                                                                     1.000 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                               eigen_iterative                             !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  MODULE taken directly from GENE (with modifications)
!!  #include "redef.h"
!!  Iterative eigenvalue solver based on SLEPc
!!
!!  All input variables like the number of eigenvalues n_ev, specification of 
!!  the method and spectral transform which_ev etc. can be set in the parameters
!!  file and/or can be overwritten at runtime with the appropriate SLEPc runtime
!!  options (see SLEPc manual) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE eigen_iterative
  USE mpi
  USE par_mod
  USE petsc_aux
  USE slepc_aux
  USE field_solver, only: get_phi
  USE linear_rhs
  IMPLICIT NONE

  COMPLEX, ALLOCATABLE, DIMENSION(:), PUBLIC :: current_evs

  PUBLIC:: ev_iterative, n_ev, evfile, calc_dt_ev 
  
  PRIVATE

  INTEGER:: evfile

#include "petscversion.h"


  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                              ev_iterative                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Computes the desired eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ev_iterative
    
    !double precision :: time_1, time_2
    INTEGER :: num_ev!, i!,ierr
    !LOGICAL:: eof
    !REAL :: time_ev
!    CHARACTER(len=10)::pc_string,ksp_string,eps_string

    IF(mype.eq.0) THEN
       WRITE(*,"(a)") '******************************************************'
       WRITE(*,"(a)") 'Starting eigenvalue computation with PETSC/SLEPc'
       WRITE(*,"(a)") '******************************************************'
    END IF
    
    ! always USE linear version of CalFullRhs for eigenvalue computation
    !IF(verbose.and.mype==0) WRITE(*,*) 45
    CALL my_SlepcInitialize(mpi_comm_world,.false.)

    !initialize_petsc_mat can't be used because derived rhs-computations
    !are allowed in this routine (e.g. scalapack interface)
    !IF(verbose.and.mype==0) WRITE(*,*) 50
    CALL MatCreateShell(mpi_comm_world,ev_size_loc,ev_size_loc,ev_size,&
         ev_size,PETSC_NULL_INTEGER,shellmat,globerr)

    !CALL initialize_calc_k
    !IF(verbose.and.mype==0) WRITE(*,*) 55
    CALL MatShellSetOperation(shellmat,MATOP_MULT,matmul_k,globerr)
    CALL VecCreateMPI(mpi_comm_world,ev_size_loc,ev_size,glob_g,globerr)

    !IF(mype==0) WRITE(*,*) "ev_size_loc,ev_size",ev_size_loc,ev_size
    CALL VecCreateMPI(mpi_comm_world,ev_size_loc,ev_size,v0_,globerr)
!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!1
    !IF(mype==0) WRITE(*,*) "globerr",globerr
    !CALL VecGetArrayF90(v0_,ptest,globerr)   
    !IF(mype==0) WRITE(*,*) "globerr",globerr
    !IF(mype==0) WRITE(*,*) "ptest",ptest
!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!1

    !IF(verbose.and.mype==0) WRITE(*,*) 72
    CALL initial_condition(init_cond)
    !IF(mype==0) WRITE(*,*) "globerr",globerr
    !IF(mype==0) WRITE(*,*) "CALL fo2pe, in ev_iterative."
    CALL fo2pe(g_1,v0_)

    !IF(verbose.and.mype==0) WRITE(*,*) 78
    CALL initialize_slepc_eps!(.true.)

    CALL VecDestroy(v0_,globerr)

    !time measurement only for solver (avoiding file I/O)
    !CALL get_systime(time_1)

    !IF(mype==0) WRITE(*,*) "Calling EPSSolve"
    !IF(verbose.and.mype==0) WRITE(*,*) 87
    CALL EPSSolve(eps_obj,globerr)

    !CALL get_systime(time_2)
    !time_ev=time_2-time_1

    CALL eps_info(num_ev)
    
    !IF(mype==0) WRITE(*,*) "Calling analyze_eigenpairs."
    !IF(verbose.and.mype==0) WRITE(*,*) 96
    IF (num_ev.gt.0) CALL analyze_eigenpairs(num_ev)
    !CALL finalize_calc_k

    CALL finalize_slepc_eps
    CALL MatDestroy(shellmat,globerr)
    CALL VecDestroy(glob_g,globerr)
    CALL my_SlepcFinalize(.false.)
    !IF(verbose.and.mype==0) WRITE(*,*) 104

    !IF(mype.eq.0) THEN
    !   WRITE(*,"(a,f10.3,a)") 'time for eigenvalue computation:',time_ev,' sec'
    !   WRITE(*,"(a)") '******************************************************'
    !END IF
  END SUBROUTINE ev_iterative


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                            analyze_eigenpairs                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Analyzes the eigenpairs stored in eps_obj, i.e. writes the eigenvalues.dat 
!!  file and applies the chosen diagnostics to the eigenvectors.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE analyze_eigenpairs(num_ev)
    INTEGER, INTENT(in):: num_ev !<number of eigenpairs in eps_obj
    COMPLEX,DIMENSION(:),ALLOCATABLE:: eigenvalue
    !INTEGER :: ev_handle
    !REAL ::norm_energy
    INTEGER:: i,j
    !LOGICAL:: dummy
    COMPLEX, DIMENSION(:,:), ALLOCATABLE :: eigenvectors
    INTEGER :: evecfile,ffile
    REAL :: norm

    ALLOCATE(eigenvalue(0:num_ev-1))
    ALLOCATE(eigenvectors(ev_size,0:num_ev-1))
    
    !CALL initialize_all_diags

    !IF(ev_out) CALL ev_out_open(ev_handle,1)
    DO i=0,num_ev-1
       CALL EPSGetEigenpair(eps_obj,i,eigenvalue(i),PETSC_NULL_SCALAR,glob_g,&
            &PETSC_NULL_OBJECT,globerr)
       IF(which_ev.eq.'shift_invert_s') THEN
          !backtransform by hand for shift_invert_s
          eigenvalue(i)=conjg(1./eigenvalue(i)+ev_shift)
       END IF
       CALL pe2fo(glob_g,g_1) 
       CALL gather_g1(g_1,eigenvectors(:,i))
       
       time=i+1
       
       !normalize eigenvectors so energy is 1.0
!       IF (xy_local) THEN
!          CALL get_energy_total(g_1,norm_energy)
!          g_1=g_1/sqrt(norm_energy)
!       END IF

!       IF(ev_out) CALL ev_out_write(g_1,i+1,num_ev,ev_handle)
!       CALL exec_all_diags(0,time,dummy,dummy,dummy)
    enddo

    CALL sort_evs(eigenvalue,eigenvectors,num_ev)
    
    time=0.0
     
    IF(mype==0) THEN
     CALL get_io_number
     ffile=io_number
     CALL get_io_number
     evecfile=io_number
     CALL get_io_number
     evfile=io_number
     OPEN(evfile, file=Trim(diagdir)//'/eigenvalues',&
                  &form='formatted', status='replace')
     IF(left_ev) THEN
       OPEN(evecfile, file=Trim(diagdir)//'/left_eigenvectors',&
                  &form='formatted', status='replace')
       OPEN(ffile, file=Trim(diagdir)//'/left_eigenvectors_uf',&
                  &form='unformatted', status='replace',access='stream')
     ELSE
       OPEN(evecfile, file=Trim(diagdir)//'/eigenvectors',&
                  &form='formatted', status='replace')
       OPEN(ffile, file=Trim(diagdir)//'/eigenvectors_uf',&
                  &form='unformatted', status='replace',access='stream')
     END IF

     WRITE(*,"(a)") 'computed eigenvalues: '
     WRITE(evfile, "(a,3F7.3)") '#eigenvalues (re/im part) at kx,ky,kz', kxmin,kymin,kzmin

     DO i=0,num_ev-1
        WRITE(*,"(2f15.8)") eigenvalue(i)
        !WRITE(evfile, "(8ES20.8)") eigenvalue(i)
        WRITE(evfile, "(8ES20.8)") aimag(eigenvalue(i)),REAL(eigenvalue(i))

        norm=REAL(sqrt(sum(conjg(eigenvectors(:,i))*eigenvectors(:,i))))
        DO j=1,nv0*nkx0*nky0*nkz0
          WRITE(evecfile,*) REAL(j-1),REAL(eigenvectors(j,i))/norm,aimag(eigenvectors(j,i))/norm
        END DO
        WRITE(evecfile,*) ""
        WRITE(evecfile,*) ""
        WRITE(ffile) i 
        WRITE(ffile) eigenvectors(:,i) 
     enddo

     CLOSE(evfile);CLOSE(evecfile)
    END IF

    DEALLOCATE(eigenvalue)
    
  END SUBROUTINE analyze_eigenpairs
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 matmul_k                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE matmul_k(locmat,petsc_in,petsc_out,ierr)
#IF !((PETSC_VERSION_MAJOR>2) && (PETSC_VERSION_MINOR>0))
#include "finclude/petsc.h"
#endif
#include "finclude/petscvec.h"   
#include "finclude/petscmat.h"
    Mat locmat
    Vec petsc_in,petsc_out
    PetscErrorCode ierr
    COMPLEX, DIMENSION(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2):: fort_g,fort_rhs

    CALL pe2fo(petsc_in,fort_g)
    CALL get_phi(fort_g)
    CALL get_rhs_lin(fort_g,phi,fort_rhs,0)
    !CALL calc_k(fort_g,fort_rhs,0)
    !IF(mype==0) WRITE(*,*) "CALL fo2pe, in matmul_k."
    CALL fo2pe(fort_rhs,petsc_out)

  END SUBROUTINE matmul_k


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 calc_dt_ev                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE calc_dt_ev
    
    INTEGER :: num_ev!, i!,ierr
    !LOGICAL:: eof

    !IF(mype.eq.0) THEN
    !   WRITE(*,"(a)") '******************************************************'
    !   WRITE(*,"(a)") 'Starting dt computation with PETSC/SLEPc'
    !   WRITE(*,"(a)") '******************************************************'
    !END IF

    which_ev='largest_magnitude'
    
    !WRITE(*,*) "mype,218",mype
    ! always USE linear version of CalFullRhs for eigenvalue computation
    CALL my_SlepcInitialize(mpi_comm_world,.false.)

    !initialize_petsc_mat can't be used because derived rhs-computations
    !are allowed in this routine (e.g. scalapack interface)
    !WRITE(*,*) "mype,224",mype
    CALL MatCreateShell(mpi_comm_world,ev_size_loc,ev_size_loc,ev_size,&
         ev_size,PETSC_NULL_INTEGER,shellmat,globerr)

    !CALL initialize_calc_k
    !WRITE(*,*) "mype,229",mype
    CALL MatShellSetOperation(shellmat,MATOP_MULT,matmul_k,globerr)
    CALL VecCreateMPI(mpi_comm_world,ev_size_loc,ev_size,glob_g,globerr)

    !WRITE(*,*) "mype,233",mype
    !IF(mype==0) WRITE(*,*) "ev_size_loc,ev_size",ev_size_loc,ev_size
    CALL VecCreateMPI(mpi_comm_world,ev_size_loc,ev_size,v0_,globerr)

!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!1
    !IF(mype==0) WRITE(*,*) "globerr",globerr
    !CALL VecGetArrayF90(v0_,ptest,globerr)   
    !IF(mype==0) WRITE(*,*) "globerr",globerr
    !IF(mype==0) WRITE(*,*) "ptest",ptest
!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!1

    !WRITE(*,*) "mype,248",mype
    CALL initial_condition(init_cond)
    !IF(mype==0) WRITE(*,*) "globerr",globerr
    !IF(mype==0) WRITE(*,*) "CALL fo2pe, in ev_iterative."
    CALL fo2pe(g_1,v0_)

    !WRITE(*,*) "mype,254",mype
    CALL initialize_slepc_eps!(.true.)

    CALL VecDestroy(v0_,globerr)

    !time measurement only for solver (avoiding file I/O)
    !CALL get_systime(time_1)

    !IF(mype==0) WRITE(*,*) "Calling EPSSolve"
    !WRITE(*,*) "mype,263",mype
    CALL EPSSolve(eps_obj,globerr)

    !CALL get_systime(time_2)
    !time_ev=time_2-time_1

    !WRITE(*,*) "mype,269",mype
    CALL eps_info(num_ev)
    
    !IF(mype==0) WRITE(*,*) "Calling analyze_eigenpairs."
    !WRITE(*,*) "mype,273",mype
    IF (num_ev.gt.0) CALL eigenpairs_dt(num_ev)
    !CALL finalize_calc_k

    !WRITE(*,*) "mype,277",mype
    CALL finalize_slepc_eps
    CALL MatDestroy(shellmat,globerr)
    CALL VecDestroy(glob_g,globerr)
    CALL my_SlepcFinalize(.false.)

    !IF(mype.eq.0) THEN
    !   WRITE(*,"(a,f10.3,a)") 'time for eigenvalue computation:',time_ev,' sec'
    !   WRITE(*,"(a)") '******************************************************'
    !END IF
  END SUBROUTINE calc_dt_ev


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                               eigenpairs_dt                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE eigenpairs_dt(num_ev)
    
    INTEGER, INTENT(in):: num_ev !<number of eigenpairs in eps_obj
    COMPLEX,DIMENSION(:),ALLOCATABLE:: eigenvalue
    REAL,DIMENSION(:),ALLOCATABLE:: max_eig
    !INTEGER :: ev_handle
    !REAL ::norm_energy
    INTEGER:: i!,mloc(1)
    !LOGICAL:: dummy
    !REAL :: dt_temp
    REAL :: real_part,imag_part

    ALLOCATE(eigenvalue(num_ev))
    ALLOCATE(max_eig(num_ev))
    
    !CALL initialize_all_diags

    !IF(ev_out) CALL ev_out_open(ev_handle,1)
    current_evs=cmplx(0.0,0.0)
    DO i=0,num_ev-1
       !WRITE(*,*) i,"318",mype
       CALL EPSGetEigenpair(eps_obj,i,eigenvalue(i+1),PETSC_NULL_SCALAR,glob_g,&
            &PETSC_NULL_OBJECT,globerr)
       !WRITE(*,*) "eigenvalue(i+1)",eigenvalue(i+1)
       !WRITE(*,*) "321",mype
       CALL pe2fo(glob_g,g_1) 
       !WRITE(*,*) "323",mype
       !CALL round_to_prec(REAL(eigenvalue(i+1)),real_part)
       !CALL round_to_prec(aimag(eigenvalue(i+1)),imag_part)
       !WRITE(*,*) "325",mype
       current_evs(i+1)=eigenvalue(i+1)
       !current_ev=eigenvalue(i+1) 
       !IF(mype==0) WRITE(*,*) "Rounded ev",eigenvalue(i+1)
       !dt_temp=100.0
       !IF(real_part.ne.0.0.or.imag_part.ne.0.0) CALL compute_stability_criterion(eigenvalue(i+1),dt_temp)
       !IF(dt_temp.lt.dt_max) THEN
       !  dt_max=dt_temp
       !  IF(mype==0) WRITE(*,*) "New dt_max:",dt_max
       !  IF(mype==0) WRITE(*,*) "Omega",eigenvalue(i+1)
       !  IF(mype==0) WRITE(*,*) "kx,ky,kz",kxmin,kymin,kzmin
       !END IF
    enddo

    !DEALLOCATE(eigenvalue)
    
  END SUBROUTINE eigenpairs_dt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 gather_g1                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE gather_g1(g_in,g_out)
   USE par_mod
   USE mpi
   IMPLICIT NONE

   COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
   COMPLEX, INTENT(out) :: g_out(0:nkx0-1,0:nky0-1,0:nkz0-1,0:nv0-1)
   COMPLEX :: g_temp(0:nkx0-1,0:nky0-1,0:nkz0-1,0:nv0-1)
   INTEGER :: ierr

   IF(np_hank.gt.1) STOP "gather_g1 not yet implemented for np_hank.gt.1"
   IF(np_spec.gt.1) STOP "gather_g1 not yet implemented for np_spec.gt.1"
   IF(np_kz.gt.1) STOP "gather_g1 not yet implemented for np_kz.gt.1"

   g_temp=cmplx(0.0,0.0)
   g_temp(:,:,:,lv1:lv2)=g_in(:,:,:,lv1:lv2,0,0)
   !IF(mype==0) write(*,*) "gather_g1"
   CALL MPI_ALLREDUCE(g_temp,g_out,nkx0*nky0*nkz0*nv0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  END SUBROUTINE gather_g1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 sort_evs                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE sort_evs(eigenvalue,eigenvectors,num_ev)
    USE par_mod
    IMPLICIT NONE

    INTEGER, INTENT(in) :: num_ev
    COMPLEX, INTENT(inout) :: eigenvalue(num_ev)    
    COMPLEX, INTENT(inout) :: eigenvectors(0:nkx0-1,0:nky0-1,0:nkz0-1,0:nv0-1,num_ev)    
    COMPLEX :: ev_bak(num_ev)    
    COMPLEX :: evec_bak(0:nkx0-1,0:nky0-1,0:nkz0-1,0:nv0-1,num_ev)    
    INTEGER :: mloc(1)
    REAL :: growth_rate(num_ev)
    INTEGER :: i

    ev_bak=eigenvalue
    evec_bak=eigenvectors
    growth_rate=REAL(eigenvalue)
    DO i=1,num_ev
      mloc=maxloc(growth_rate)
      !IF(mype==0) WRITE(*,*) mloc
      eigenvalue(i)=ev_bak(mloc(1))
      eigenvectors(:,:,:,:,i)=evec_bak(:,:,:,:,mloc(1))
      growth_rate(mloc(1))=-1.0e15
    END DO

  END SUBROUTINE sort_evs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                              round_to_prec                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE round_to_prec(r_in,r_out)
    USE par_mod
    IMPLICIT NONE

    REAL, INTENT(in) :: r_in
    REAL, INTENT(out) :: r_out
    INTEGER :: integ
    REAL :: prec_factor

    prec_factor=0.1*ev_prec
    !WRITE(*,*) "mype,r_in",mype,r_in
    !WRITE(*,*) "ev_prec",ev_prec
    integ=nint(r_in/prec_factor)
    !WRITE(*,*) "mype,integ",mype,integ
    r_out=REAL(integ)*prec_factor
    !WRITE(*,*) "mype,r_out",mype,r_out

  END SUBROUTINE round_to_prec
  
END MODULE eigen_iterative

