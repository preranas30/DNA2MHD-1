MODULE eigen_direct
   USE par_mod
   USE linear_rhs
   USE field_solver, only: get_phi


   PUBLIC :: calc_dt_ev_lapack

   PRIVATE
   COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: mat
   COMPLEX, ALLOCATABLE, DIMENSION(:) :: current_evs_lp

   CONTAINS

SUBROUTINE calc_dt_ev_lapack(evs_out)

   COMPLEX, INTENT(OUT) :: evs_out(nv0*nh0*nspec)
   INTEGER :: i

   ALLOCATE(mat(nv0*nh0*nspec,nv0*nh0*nspec))
   ALLOCATE(current_evs_lp(nv0*nh0*nspec))

   CALL set_full_mat_from_rhs
   !DO i=0,nv0*nh0*nspec
   !  IF(mype==0) WRITE(*,*) mat(:,i)
   !END DO
   CALL serial_eigensolve   
   evs_out=current_evs_lp

   DEALLOCATE(mat)
   DEALLOCATE(current_evs_lp)

END SUBROUTINE calc_dt_ev_lapack

SUBROUTINE serial_eigensolve

  IMPLICIT NONE

  CHARACTER(LEN=1)  :: JOBVL,JOBVR
  INTEGER :: info
  REAL, ALLOCATABLE, DIMENSION(:) :: rwork
  COMPLEX, ALLOCATABLE, DIMENSION(:) :: work
  INTEGER :: lwork
  COMPLEX :: dummy

  !WRITE(*,*) "Beginning serial eigensolve"
  JOBVL='N'
  JOBVR='N'
  !IF(left_vec) JOBVL='V'
  !IF(right_vec) JOBVR='V'

  ALLOCATE(rwork(2*nv0*nh0*nspec))
  ALLOCATE(work(1))
  lwork=-1
  info=0
  CALL ZGEEV(JOBVL,JOBVR,nv0*nh0*nspec,mat,nv0*nh0*nspec,current_evs_lp,dummy,nv0*nh0*nspec,&
      dummy,nv0*nh0*nspec, work,lwork,rwork,info)

  lwork=work(1)
  !WRITE(*,*) "lwork",lwork
  DEALLOCATE(work)
  ALLOCATE(work(lwork))

  CALL ZGEEV(JOBVL,JOBVR,nv0*nh0*nspec,mat,nv0*nh0*nspec,current_evs_lp,dummy,nv0*nh0*nspec,&
      dummy,nv0*nh0*nspec,work,lwork,rwork,info)
  !write(*,*) "Finished serial eigensolve"

END SUBROUTINE serial_eigensolve

SUBROUTINE set_full_mat_from_rhs

  USE par_mod

  COMPLEX :: rhs_(lv1:lv2,lh1:lh2,ls1:ls2)
  INTEGER :: i,rhs_lin_version_bak,mat_ind

  rhs_lin_version_bak=rhs_lin_version
  rhs_lin_version=1
  DO i=0,nv0-1
    DO j=0,nh0-1
      DO k=0,nspec-1
        mat_ind=k*nh0*nv0+j*nv0+i+1
        g_1=cmplx(0.0,0.0)
        g_1(0,0,0,i,j,k)=cmplx(1.0,0.0)
        CALL get_phi(g_1)
        CALL get_rhs_lin(g_1,phi,rhs_,0) 
        !mat(:,mat_ind)=rhs_(0:nv0-1)
        call put_rhs_in_mat(rhs_,mat_ind)

      END DO
    END DO
  END DO
  rhs_lin_version=rhs_lin_version_bak

END SUBROUTINE set_full_mat_from_rhs

SUBROUTINE put_rhs_in_mat(rhs_in,mat_ind)

  implicit none

  COMPLEX, INTENT(IN) :: rhs_in(nv0*nh0*nspec)
  INTEGER, INTENT(IN) :: mat_ind

  mat(:,mat_ind)=rhs_in

END SUBROUTINE put_rhs_in_mat

END MODULE eigen_direct
