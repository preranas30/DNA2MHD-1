!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 29/12/2012                                                                !!
!!                              ee_petsc_aux.f90                             !!
!!                                                                           !!
!!  petsc_aux                                                                !!
!!                                                                     1.000 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                                  petsc_aux                                !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Module modified from GENE version.
!!  Note taken from GENE with modifications
!!  Routines for handling the KSP and PC objects used in PETSc (and SLEPc) 
!!  and routines for the conversion between fortran and PETSc vector formats
!!
!!  It provides a MODULE variable ksp_obj which should be used by all routines 
!!  based on PETSc and a variable shellmat which defines the operator for the 
!!  matrix-free PETSc methods used in GENE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
#include "petscversion.h"

MODULE petsc_aux
  USE par_mod
  !USE par_other, only: dt, prec, print_ini_msg
  !USE par_in
  !USE eigen_parameters
  !USE discretization
  !USE calc_rhs
  USE petsc
  USE petscvec
  USE petscmat
  USE petscksp
  USE petscpc
  USE field_solver, only: get_phi
  USE linear_rhs
  IMPLICIT NONE

  PUBLIC:: pe2fo, fo2pe, initialize_petsc_mat, finalize_petsc_mat, initialize_petsc_vec, finalize_petsc_vec!,&
       !initialize_petsc_flags
  PUBLIC:: glob_g, glob_rhs, v0_, ksp_obj, shellmat, globerr,a_mat,b_mat!, impl_shift
  PUBLIC:: matop_mult, different_nonzero_pattern, same_preconditioner, same_nonzero_pattern, petsc_null_object,&
       petsc_null_character, petsc_null_scalar, petsc_decide, petsc_null, petsc_null_integer, petsc_comm_world,&
       kspreason,ptest!,petsc_t,petsc_f

  !INTEGER, PUBLIC :: ev_size_loc,ev_size,n_ev,ev_n_test
  INTEGER, PUBLIC :: ev_n_test
  INTEGER, PUBLIC :: ksp_max_it=0
  CHARACTER(len=100), PUBLIC :: pc_type='DEFAULT'
  INTEGER, PUBLIC :: it_ev=-1

  PRIVATE
#IF (PETSC_VERSION_MINOR>0)
#include "finclude/petscsysdef.h"
#ELSE
#include "finclude/petscdef.h"
#endif
#include "finclude/petscvecdef.h"   
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscpcdef.h"
 
  Mat shellmat, L_g_mat,a_mat,b_mat
  Vec glob_g, glob_rhs, v0_
  KSP ksp_obj
  KSPConvergedReason kspreason
  PetscErrorCode globerr
  PetscScalar,POINTER:: ptest(:)
  !PetscBool petsc_t, petsc_f

!petsc<3.2 does not know PetscBool but PetscTruth
#ifndef PetscBool
#define PetscBool PetscTruth
#endif


  !LOGICAL:: impl_shift=.false.

CONTAINS
  !>Converts an array of type g_1 from PETSc to Fortran format
  SUBROUTINE pe2fo(petg,locg)
    Vec petg !<petsc representation of a g_1-type array
    !>fortran array of type g_1
    COMPLEX, DIMENSION(ev_size_loc),INTENT(out):: locg 
    PetscScalar,pointer:: arr_point(:)
    INTEGER:: ierr


    !IF(mype==0) WRITE(*,*) "Start of pe2fo"
    CALL VecGetArrayF90(petg,arr_point,ierr)   

    IF (evenyz) THEN
      CALL insert_zero_kx(arr_point,locg)
    ELSE
       locg=arr_point
    END IF

    CALL VecRestoreArrayF90(petg,arr_point,ierr)
    
  END SUBROUTINE pe2fo

  !>Converts a vector of type g_1 from Fortran to PETSc format
  SUBROUTINE fo2pe(locg,petg) !convert g from fortran to petsc format
    !>fortran array of type g_1
    COMPLEX, DIMENSION(ev_size_loc),INTENT(in):: locg 
    Vec petg !<petsc representation of a g_1-type array
    PetscScalar,POINTER:: arr_point(:)
    INTEGER:: ierr

    !IF(mype==0) WRITE(*,*) "Start of fo2pe"
    CALL VecGetArrayF90(petg,arr_point,ierr)   
    !IF(mype==0) WRITE(*,*) "ierr",ierr
    !IF(mype==0) WRITE(*,*) "arr_point",arr_point

    IF (evenyz) THEN
      CALL remove_zero_kx(locg,arr_point)
    ELSE
       !WRITE(*,*) "sum of locg",sum(locg)
       arr_point=locg
    END IF

    CALL VecRestoreArrayF90(petg,arr_point,ierr)

  END SUBROUTINE fo2pe

  SUBROUTINE remove_zero_kx(locg,locsmallg)
    COMPLEX, DIMENSION(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2), INTENT(in):: locg 
    COMPLEX, DIMENSION(0:nkx0-1,0:nky0-2,lkz1:lkz2,lv1:lv2), INTENT(out):: locsmallg 
    INTEGER :: i,j,k
    INTEGER :: j_shift,k_shift

    !DO ind=1,lklmn0
    !   locsmallg(lg1:hkx,ind)=locg(lg1:hkx,ind)
    !   locsmallg(lkx-1:lg2-1,ind)=locg(lkx:lg2,ind)
    !enddo
    !remove extra k modes in y and z coordinates
    !STOP "Must check locsmallg g first--check nky0-2,nkz0-1"

   IF(np_hank.gt.1) STOP "remove_zero_kx not yet implemented for np_hank.gt.1"
   IF(np_spec.gt.1) STOP "remove_zero_kx not yet implemented for np_spec.gt.1"
   IF(np_kz.gt.1) STOP "remove_zero_kx not yet implemented for np_kz.gt.1"
    
    DO i=0,nkx0-1
      DO j=0,nky0-1

        IF(j.gt.hky_ind) THEN
          j_shift=-1    
        ELSE
          j_shift=0 
        END IF 

        DO k=0,nkz0-1

          IF(k.gt.hkz_ind) THEN
            k_shift=-1    
          ELSE
            k_shift=0 
          END IF 
          IF(j.ne.hky_ind+1.and.k.ne.hkz_ind+1) locsmallg(i,j+j_shift,k+k_shift,:)=locg(i,j,k,:)

        END DO
      END DO
    END DO

  END SUBROUTINE remove_zero_kx

  SUBROUTINE insert_zero_kx(locsmallg, locg)
    COMPLEX, DIMENSION(0:nkx0-1,0:nky0-2,0:nkz0-2,lv1:lv2), INTENT(in):: locsmallg 
    COMPLEX, DIMENSION(0:nkx0-1,0:nky0-1,0:nkz0-1,lv1:lv2), INTENT(out):: locg 
    !COMPLEX, DIMENSION(lg1:lg2-1,lklmn0), INTENT(in):: locsmallg 
    !COMPLEX, DIMENSION(lg1:lg2,lklmn0), INTENT(out):: locg 
    INTEGER:: i,j,k,j_shift,k_shift

    !STOP "Must check locsmallg g first--check nky0-2,nkz0-1"
    !remove extra k modes in y and z coordinates

   IF(np_hank.gt.1) STOP "remove_zero_kx not yet implemented for np_hank.gt.1"
   IF(np_spec.gt.1) STOP "remove_zero_kx not yet implemented for np_spec.gt.1"
   IF(np_kz.gt.1) STOP "remove_zero_kx not yet implemented for np_kz.gt.1"

    DO i=0,nkx0-1
      DO j=0,nky0-1

        IF(j.gt.hky_ind) THEN
          j_shift=-1    
        ELSE
          j_shift=0 
        END IF 

        DO k=0,nkz0-1

          IF(k.gt.hkz_ind) THEN
            k_shift=-1    
          ELSE
            k_shift=0 
          END IF 

          IF(j.ne.hky_ind+1.and.k.ne.hkz_ind+1) THEN
              locg(i,j,k,:)=locsmallg(i,j+j_shift,k+k_shift,:)
          ELSE
              locg(i,j,k,:)=0.0
          END IF
        END DO
      END DO
    END DO

    !DO ind=1,lklmn0
    !   locg(lg1:hkx,ind)=locsmallg(lg1:hkx,ind)
    !   locg(hkx+1,ind)=0
    !   locg(lkx:lg2,ind)=locsmallg(lkx-1:lg2-1,ind)
    !enddo

  END SUBROUTINE insert_zero_kx

  !>Initializes flags for calling petsc routines (only logicals right now..)
  !SUBROUTINE initialize_petsc_flags
  ! petsc_t = PETSC_TRUE
  !  petsc_f = PETSC_FALSE
  !END SUBROUTINE initialize_petsc_flags

  !>Initializes shellmat and creates the necessary PETSc vectors
  SUBROUTINE initialize_petsc_mat
    !LOGICAL,INTENT(in):: mat_invert 
    INTEGER :: petsc_scalar=0, petsc_complex=0


   IF(np_hank.gt.1) STOP "initialize_petsc_mat not yet implemented for np_hank.gt.1"
   IF(np_spec.gt.1) STOP "initialize_petsc_mat not yet implemented for np_spec.gt.1"
   IF(np_kz.gt.1) STOP "initialize_petsc_mat not yet implemented for np_kz.gt.1"
    !IF(mype==0) WRITE(*,*) "Start of initialize_petsc_mat."
#ifdef PETSC_SCALAR
    petsc_scalar = PETSC_SCALAR
    petsc_complex = PETSC_COMPLEX
#endif

    IF ((mype==0)) &
         & WRITE(*,'(A,I1,A,I1,A,I1,A,I2.2)') "Using PETSc version ",PETSC_VERSION_MAJOR,&
         &'.',PETSC_VERSION_MINOR,'.',PETSC_VERSION_SUBMINOR,'-p',PETSC_VERSION_PATCH
    IF (petsc_scalar.ne.petsc_complex) &
       & STOP 'Sorry, you need to configure PETSc using -with-scalar-type=COMPLEX'

!    CALL initialize_CalFullRhs(.true.)
    CALL MatCreateShell(mpi_comm_world,ev_size_loc,ev_size_loc,ev_size,&
         ev_size,PETSC_NULL_INTEGER,shellmat,globerr)
    CALL MatShellSetOperation(shellmat,MATOP_MULT,imp_matmul,globerr)

    CALL VecCreateMPI(MPI_COMM_WORLD,ev_size_loc,ev_size,glob_g,globerr)

    !impl_shift=mat_invert

  END SUBROUTINE initialize_petsc_mat

  !>Deletes shellmat and the PETSc vectors
  SUBROUTINE finalize_petsc_mat

    !IF(mype==0) WRITE(*,*) "Start of finalize_petsc_mat."
    CALL MatDestroy(shellmat,globerr)
    CALL VecDestroy(glob_g,globerr)
    !CALL finalize_CalFullRhs

  END SUBROUTINE finalize_petsc_mat

  SUBROUTINE initialize_petsc_vec
    !IF(mype==0) WRITE(*,*) "Start of initialize_petsc_vec."
    CALL VecCreateMPI(MPI_COMM_WORLD,ev_size_loc,ev_size,glob_rhs,globerr)
  END SUBROUTINE initialize_petsc_vec

  SUBROUTINE finalize_petsc_vec
    !IF(mype==0) WRITE(*,*) "Start of finalize_petsc_vec."
    CALL VecDestroy(glob_rhs,globerr)
  END SUBROUTINE finalize_petsc_vec

  !>Defines the matrix-free matrix vector multiplication by directly using CalFullRhs
  SUBROUTINE imp_matmul(locmat,vec_in,vec_out,ierr)
    
    Mat locmat !<PETSc matrix object which is defined in this routine
    Vec vec_in !<PETSc input vector
    Vec vec_out !<PETSc vector containing the result of locmat*vec_in
    PetscErrorCode ierr !<PETSC error code
    COMPLEX, DIMENSION(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2):: loc_g,loc_res

    !IF(mype==0) WRITE(*,*) "Start of imp_matmul."
    CALL pe2fo(vec_in,loc_g)

    CALL get_phi(loc_g)
    CALL get_rhs_lin(loc_g,phi,loc_res,0)

    !!!CALL CalFullRhs(loc_g,loc_res,0)
    !IF (impl_shift) THEN
    !   !compute (1-dt*L)g
    !   loc_res=loc_g-dt*loc_res
    !END IF

    !IF(mype==0) WRITE(*,*) "Calling fo2pe in imp_matmul."
    CALL fo2pe(loc_res,vec_out)
    ierr=0

  END SUBROUTINE imp_matmul

END MODULE petsc_aux


