!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 29/12/2012                                                                !!
!!                               ee_slepc_aux.f90                            !!
!!                                                                           !!
!!  slepc_aux                                                                !!
!!                                                                     1.000 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                                  slepc_aux                                !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Module taken largely from GENE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#include "redef.h"
#include "petscversion.h"

!>Routines for handling the EPS and ST objects used in SLEPc
!!
!!It provides a MODULE variable eps_obj which should be used by all routines based on SLEPc.
MODULE slepc_aux
  USE par_mod
  USE petsc_aux
  USE petscvec  
IMPLICIT NONE

#include "finclude/slepc.h"

#IF (PETSC_VERSION_MINOR==0)
#include "finclude/slepcst.h"
#include "finclude/slepceps.h"
#endif
  PUBLIC:: initialize_slepc_eps, finalize_slepc_eps, eps_info, eps_obj, &
       & my_SlepcInitialize, my_SlepcFinalize, n_ch !, slepc_restartable
  PUBLIC :: slepc_initialized_from_outside
  PRIVATE

  INTEGER:: n_ch
  LOGICAL:: write_mem_req=.true.,slepc_initialized_from_outside=.false.

  EPS eps_obj  !<SLEPc object to describe the eigenvalue problem
  ST st_obj  !<SLEPc object to describe the spectral transform (IF any)
  PC pc_obj

CONTAINS

!>Initializes the eps_obj with the appropriate problem type, dimensions etc. 
!!It also sets the corresponding spectral transform (st_obj) IF needed
  SUBROUTINE initialize_slepc_eps

    !STType st_type
    !MatStructure flag,flag1
    LOGICAL :: ev_run

    ev_run=.false.

    CALL EPSCreate(PETSC_COMM_WORLD,eps_obj,globerr)

    !problem definitions
    CALL EPSSetOperators(eps_obj,shellmat,PETSC_NULL_OBJECT,globerr)
    CALL EPSSetProblemType(eps_obj,EPS_NHEP,globerr)
  
    CALL EPSSetDimensions(eps_obj,n_ev,ev_n_test,PETSC_DECIDE,globerr)  
    CALL EPSSetTolerances(eps_obj,ev_prec,ev_max_it,globerr)

    !SELECT eigenvalues and method
    IF (.not.ev_run) THEN
       !slepc only used to determine timestep
       CALL EPSSetWhichEigenPairs(eps_obj,EPS_LARGEST_REAL,globerr)
    ELSE
#IF (PETSC_VERSION_MINOR>0)
       CALL EPSSetInitialSpace(eps_obj,n_ch,v0_,globerr)
#ELSE
       CALL EPSSetInitialVector(eps_obj,v0_,globerr)
#endif
       IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
       IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
       IF(mype==0) WRITE(*,*) "Using which_ev=",which_ev
       IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
       IF(mype==0) WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
       SELECT CASE (which_ev)
       CASE('largest_real')
          CALL EPSSetWhichEigenPairs(eps_obj,EPS_LARGEST_REAL,globerr)
          CALL EPSSetType(eps_obj,EPSKRYLOVSCHUR,globerr)
       CASE('smallest_real')
          CALL EPSSetWhichEigenPairs(eps_obj,EPS_SMALLEST_REAL,globerr)
          CALL EPSSetType(eps_obj,EPSKRYLOVSCHUR,globerr)
       CASE('largest_magnitude')
          CALL EPSSetWhichEigenPairs(eps_obj,EPS_LARGEST_MAGNITUDE,globerr)
          CALL EPSSetType(eps_obj,EPSKRYLOVSCHUR,globerr)
       CASE('harmonic')
#IF (PETSC_VERSION_MINOR>0)
          CALL EPSSetWhichEigenPairs(eps_obj,EPS_TARGET_MAGNITUDE,globerr)
#endif
          CALL EPSSetTarget(eps_obj,ev_shift,globerr)
          CALL EPSSetExtraction(eps_obj,EPS_HARMONIC,globerr)
          CALL EPSSetType(eps_obj,EPSKRYLOVSCHUR,globerr)
#IF (PETSC_VERSION_MINOR>0)
       CASE('jd')
          !CALL EPSSetWhichEigenPairs(eps_obj,EPS_TARGET_MAGNITUDE,globerr)
          CALL EPSSetWhichEigenPairs(eps_obj,EPS_LARGEST_IMAGINARY,globerr)
          !CALL EPSSetTarget(eps_obj,ev_shift,globerr)
          CALL EPSSetType(eps_obj,EPSJD,globerr)
!          CALL EPSSetExtraction(eps_obj,EPS_HARMONIC_RIGHT,globerr)
          CALL EPSSetExtraction(eps_obj,EPS_HARMONIC,globerr)
          CALL EPSJDSetRestart(eps_obj,8,0,globerr)
          CALL EPSGetST(eps_obj,st_obj,globerr)
          CALL STSetShift(st_obj,ev_shift,globerr)
!          CALL STSetUp(st_obj,globerr)
          CALL STGetKSP(st_obj,ksp_obj,globerr)
          CALL KSPSetType(ksp_obj,KSPBCGS,globerr)
          IF(ksp_max_it.eq.0) THEN
             IF (pc_type.eq.'none') THEN
                ksp_max_it=80
             ELSE
                ksp_max_it=5
             END IF
          END IF
          CALL KSPSetTolerances(ksp_obj,1e-16,1e-10,1e4,ksp_max_it,globerr) 
       CASE('gd')
          CALL EPSSetWhichEigenPairs(eps_obj,EPS_TARGET_MAGNITUDE,globerr)
          CALL EPSSetTarget(eps_obj,ev_shift,globerr)
          CALL EPSSetExtraction(eps_obj,EPS_HARMONIC,globerr)
          CALL EPSSetType(eps_obj,EPSGD,globerr)
          CALL EPSGetST(eps_obj,st_obj,globerr)
          CALL STSetShift(st_obj,ev_shift,globerr)
#endif
       CASE('shift_invert')
#IF (PETSC_VERSION_MINOR>0)
          CALL EPSSetWhichEigenPairs(eps_obj,EPS_TARGET_MAGNITUDE,globerr)
          CALL EPSSetTarget(eps_obj,ev_shift,globerr)
          CALL EPSSetType(eps_obj,EPSKRYLOVSCHUR,globerr);
#ELSE
          CALL EPSSetWhichEigenPairs(eps_obj,EPS_LARGEST_MAGNITUDE,globerr)
#endif
          CALL EPSGetST(eps_obj,st_obj,globerr)
#IF (PETSC_VERSION_MINOR>0)
          !CALL EPSSetTrueResidual(eps_obj,PETSC_TRUE,globerr)
          CALL EPSSetTrueResidual(eps_obj,PETSC_FALSE,globerr)
          CALL STSetType(st_obj,STSINVERT,globerr)
          CALL STSetMatMode(st_obj,ST_MATMODE_SHELL,globerr)
#ELSE
          CALL STSetType(st_obj,STSINV,globerr)
          CALL STSetMatMode(st_obj,STMATMODE_SHELL,globerr)
#endif
          CALL STSetFromOptions(st_obj,globerr)
          CALL STSetShift(st_obj,ev_shift,globerr)
          CALL STGetKSP(st_obj,ksp_obj,globerr)
          !CALL KSPSetFromOptions(ksp_obj,globerr)
          CALL KSPSetType(ksp_obj,KSPGMRES,globerr)
          !CALL STSetFromOptions(st_obj,globerr)
          CALL STSetUp(st_obj,globerr)
       CASE ('runtime')  
          !DO nothing; all options taken from command line
       CASE('all_lapack')
          CALL EPSSetType(eps_obj,EPSLAPACK,globerr)
#ifdef WITHSCAL
       CASE('shift_invert_s')
          !spectral transform is done explicitly
          CALL EPSSetWhichEigenPairs(eps_obj,EPS_LARGEST_MAGNITUDE,globerr)
          CALL EPSSetType(eps_obj,EPSKRYLOVSCHUR,globerr);
#endif
       CASE DEFAULT
          IF (mype.eq.0) PRINT*, 'no valid which_ev selected'
          STOP
       END SELECT

       !the eps options can be overwritten at runtime!
       CALL EPSSetFromOptions(eps_obj,globerr)
   
       !setup the preconditioner
       CALL EPSGetST(eps_obj,st_obj,globerr)
       CALL STGetKSP(st_obj,ksp_obj,globerr)
       CALL KSPGetPC(ksp_obj,pc_obj,globerr)
       CALL PCSetType(pc_obj,PCNONE,globerr)
       !CALL KSPSetType(ksp_obj,KSPGMRES,globerr)
       !CALL KSPGMRESSetRestart(ksp_obj,12,globerr)
       !CALL KSPSetTolerances(ksp_obj,1e-5,1e-35,100000,globerr)
!       IF (pc_type .ne. PCNONE) THEN
!          CALL STGetType(st_obj,st_type,globerr)
!          !IF (mype.eq.0) WRITE(*,*) 'der ST_TYPE ist',st_type
!          SELECT CASE(st_type)
!#IF (PETSC_VERSION_MINOR>0)
!          CASE(STPRECOND)
!             CALL initialize_L_g
!!            CALL STGetShift(st_obj,ev_shift,globerr)
!!            CALL MatShift(L_g_mat,ev_shift,globerr)
!             CALL PCSetType(pc_obj,PCJACOBI,globerr)
!             CALL STPrecondSetMatForPC(st_obj,L_g_mat,globerr)
!             CALL EPSSetUp(eps_obj,globerr)
!             !CALL finalize_L_g
!!            CALL PCGetOperators(pc_obj,b_mat,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,globerr)
!             !CALL KSPGetOperators(ksp_obj,a_mat,b_mat,flag1,globerr)
!             !CALL KSPSetOperators(ksp_obj,a_mat,L_g_mat,flag,globerr)
!             !CALL KSPSetUp(ksp_obj,globerr)
!             CALL set_pc_obj
!             CALL finalize_L_g
!             !CALL PCSetUp(pc_obj,globerr)
!             !CALL KSPSetUp(ksp_obj,globerr)
!             !CALL EPSSetUp(eps_obj,globerr)
!          CASE(STSINVERT)
!#ELSE
!          CASE(STSINV)
!#endif
!            ! CALL EPSSetUp(eps_obj,globerr)
!             CALL initialize_L_g
!             flag = DIFFERENT_NONZERO_PATTERN
!             CALL KSPGetOperators(ksp_obj,a_mat,b_mat,flag1,globerr)
!             !CALL PetscObjectReference(a_mat,globerr)
!             CALL KSPSetOperators(ksp_obj,a_mat,L_g_mat,flag,globerr)
!            ! CALL KSPGetPC(ksp_obj,pc_obj,globerr)
!            ! CALL PCSetOperators(pc_obj,a_mat,L_g_mat,flag,globerr)
!            !CALL MatDestroy(a_mat,globerr)
!            ! CALL finalize_L_g
!             CALL KSPSetFromOptions(ksp_obj,globerr)
!             CALL KSPSetUp(ksp_obj,globerr)
!             CALL set_pc_obj
!             CALL finalize_L_g
!             CALL EPSSetUp(eps_obj,globerr)
!          END SELECT
!       END IF
    END IF
    !IF (mype.eq.0) CALL EPSView(eps_obj,PETSC_VIEWER_STDOUT_SELF,globerr)
  END SUBROUTINE initialize_slepc_eps

  !>Deletes the eps_obj
  SUBROUTINE finalize_slepc_eps
!    IF (L_g_initialized) CALL finalize_L_g
    CALL EPSDestroy(eps_obj,globerr)
  END SUBROUTINE finalize_slepc_eps

  !>Routine to get and output the number of computed eigenvalues and SLEPc iterations needed
  SUBROUTINE eps_info(ev_number)
    INTEGER, INTENT(out):: ev_number !< number of eigenvalues contained in eps_obj
    
    CALL EPSGetIterationNumber(eps_obj,it_ev,globerr)
    !IF (mype.eq.0) WRITE(*,"(a,i6)") 'number of iterations:',it_ev
    
    CALL EPSGetConverged(eps_obj,ev_number,globerr)
    !IF (mype.eq.0) WRITE(*,"(a,i6,a)") 'calculated ',ev_number,' eigenvectors'
       
    IF (ev_number.eq.0) THEN
       IF (mype.eq.0) THEN
          WRITE(*,"(a)") '***NO EIGENVALUES CALCULATED***'
       END IF
    END IF

  END SUBROUTINE eps_info
  
  SUBROUTINE my_SlepcInitialize(slepc_communicator,call_from_outside)
    INTEGER,INTENT(in):: slepc_communicator
    LOGICAL, INTENT(in):: call_from_outside

    !IF(mype==0) WRITE(*,*) "Initializing slepc."

    IF (.not.slepc_initialized_from_outside) THEN
       PETSC_COMM_WORLD=slepc_communicator
       CALL SlepcInitialize(PETSC_NULL_CHARACTER,globerr)
    END IF
    IF(call_from_outside) slepc_initialized_from_outside=.true.

  END SUBROUTINE my_SlepcInitialize

  SUBROUTINE my_SlepcFinalize(call_from_outside) 
    LOGICAL, INTENT(in):: call_from_outside

    IF(call_from_outside.or.(.not.slepc_initialized_from_outside)) THEN
       CALL SlepcFinalize(globerr)      
    END IF
  END SUBROUTINE my_SlepcFinalize


END MODULE slepc_aux
