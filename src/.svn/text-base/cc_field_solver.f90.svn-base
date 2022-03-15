!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 29/12/2012                                                                !!
!!                             cc_field_solver.f90                           !!
!!                                                                           !!
!!  field_solver                                                             !!
!!  -- get_phi                                                               !!
!!  -- initialize_phi                                                        !!
!!                                                                     1.000 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                               field_solver                                !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE field_solver
  USE par_mod
  USE aux_func 
  USE mpi
  USE flr_effects
  USE hk_effects
  PUBLIC :: initialize_phi,get_phi


  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                  get_phi                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_phi(g_in)

  IMPLICIT NONE

  COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)

  IF (nspec==1) THEN
     IF (mu_integrated) THEN
      CALL get_phi_ae_mu0(g_in)
     ELSE
        IF (hankel) THEN
          CALL get_phi_ae_hk(g_in)
        ELSE
          CALL get_phi_ae_mu(g_in)
        END IF
     END IF 
  ELSE
!not yet ready for other cases
      STOP "field solver only available for with adiabatic e's."
  END IF

END SUBROUTINE get_phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                  get_phi_ae_mu0                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!field solver for adiabatic electrons (ae) and mu-integrated (mu0) version
  SUBROUTINE get_phi_ae_mu0(g_in)
    IMPLICIT NONE

    COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
    INTEGER :: i,j,k,ierr
    COMPLEX :: phi_fsa
  
  IF(mype==0) THEN

    !DO i=0,nkx0-1
    !  IF((kxgrid(i).ne.0.0).and.(lkz1==0)) THEN
    !     phi_fsa=pi**0.25*J0a(i,0)*g_in(i,0,0,0,0,0)/(1.0-gamma0(kxgrid(i)**2))
    !  END IF
    !END DO
    phi_fsa=0.0

    DO i=0,nkx0-1
      IF((kxgrid(i).ne.0.0).and.(lkz1==0)) THEN
        phi_fsa=pi**0.25*J0a(i,0)*g_in(i,0,0,0,0,0)/(1.0-gamma0(kxgrid(i)**2))
      END IF
       DO j=0,nky0-1
         DO k=lkz1,lkz2
          phi(i,j,k)=pi**0.25*J0a(i,j)*g_in(i,j,k,0,0,0)/phi_denom(i,j)       
          IF(spatial2d.and.kygrid(j)==0) THEN
            phi(i,j,k)=phi(i,j,k)+(Ti0Te*phi_fsa*etg_factor)/phi_denom(i,j)
          ELSEIF(kygrid(j)==0.and.kzgrid(k)==0) THEN
            phi(i,j,k)=phi(i,j,k)+(Ti0Te*phi_fsa*etg_factor)/phi_denom(i,j)
          END IF
          !WRITE(100,*) i,j,k,abs(phi(i,j,k)),phi_denom(i,j) 
        END DO
      END DO
    END DO 

  END IF !mype==0

  !Broadcast phi to all procs
  CALL MPI_BCAST(phi,nkx0*nky0*nkz0,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr) 

 END SUBROUTINE get_phi_ae_mu0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                  get_phi_ae_hk                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!field solver for adiabatic electrons (ae) and hk version
  SUBROUTINE get_phi_ae_hk(g_in)
    IMPLICIT NONE

    COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
    COMPLEX             :: g_glob(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,0:nh0-1,ls1:ls2)
    REAL, DIMENSION(0:nkx0-1,0:nky0-1) :: kperp
    INTEGER :: i,j,k,ierr
    INTEGER, DIMENSION(1) :: p_kperp
    COMPLEX :: phi_fsa

    kperp = sqrt(2*kperp2)  
    CALL MPI_GATHER(g_in,size(g_in),MPI_DOUBLE_COMPLEX,g_glob,size(g_in),MPI_DOUBLE_COMPLEX,0,MPI_COMM_HANK,ierr)

  IF(mype==0) THEN

    phi_fsa=0.0

    DO i=0,nkx0-1
      IF((kxgrid(i).ne.0.0).and.(lkz1==0)) THEN
        call nearest_value(kperp(i,0),p_kperp)
        !write(*,*) "value", kperp(i,0),hkgrid(p_kperp(1)), p_kperp(1)
        phi_fsa=pi**(0.25)*g_glob(i,0,0,0,p_kperp(1),0)/(1.0-gamma0(kxgrid(i)**2))
      END IF
       DO j=0,nky0-1
         DO k=lkz1,lkz2
          call nearest_value(kperp(i,j),p_kperp)
          !write(*,*) "value", kperp(i,j),hkgrid(p_kperp(1))
          phi(i,j,k)=2*pi**(0.25)*g_glob(i,j,k,0,p_kperp(1),0)/phi_denom(i,j)       
          IF(kygrid(j)==0.and.spatial2d) THEN
            phi(i,j,k)=phi(i,j,k)+(2*Ti0Te*phi_fsa*etg_factor)/phi_denom(i,j)
          ELSEIF(kygrid(j)==0.and.kzgrid(k)==0)  THEN
            phi(i,j,k)=phi(i,j,k)+(2*Ti0Te*phi_fsa*etg_factor)/phi_denom(i,j)
          END IF
        END DO
      END DO
    END DO! 

  END IF !mype==0

  !Broadcast phi to all procs
  CALL MPI_BCAST(phi,nkx0*nky0*nkz0,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr) 

 END SUBROUTINE get_phi_ae_hk


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                  get_phi_ae_mu                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!field solver for adiabatic electrons (ae)
  SUBROUTINE get_phi_ae_mu(g_in)
    IMPLICIT NONE

    COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
    COMPLEX  :: g_mu0(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,ls1:ls2)
    REAL  :: pref(0:nkx0-1,0:nky0-1,lh1:lh2)
    INTEGER :: i,j,k,h,ierr
    COMPLEX :: phi_fsa
  
    pref(:,:,lh1:lh2) = J0_fac(:,:,lh1:lh2)
    g_mu0  = CMPLX(0.0,0.0)
    phi_fsa= CMPLX(0.0,0.0)
    !First I integrated over vperp and then I do the flux surface
    CALL integral_v(g_in,pi*pref,g_mu0)
    
  IF(lv1.le.0.and.lv2.ge.0) THEN
     DO i=0,nkx0-1
       IF((kxgrid(i).ne.0.0).and.(lkz1==0)) THEN
         phi_fsa=pi**(0.25)*g_mu0(i,0,0,0,0)/(1.0-gamma0(kxgrid(i)**2))
       END IF
       DO j=0,nky0-1
         DO k=lkz1,lkz2
          phi(i,j,k)=pi**(0.25)*g_mu0(i,j,k,0,0)/phi_denom(i,j)       
          IF(kygrid(j)==0.and.spatial2d) THEN
            phi(i,j,k)=phi(i,j,k)+(Ti0Te*phi_fsa*etg_factor)/phi_denom(i,j)
          ELSEIF(kygrid(j)==0.and.kzgrid(k)==0)  THEN
            phi(i,j,k)=phi(i,j,k)+(Ti0Te*phi_fsa*etg_factor)/phi_denom(i,j)
          END IF
       END DO
      END DO
     END DO 
 END IF !lv1 condition 
  !Broadcast phi to all procs
  CALL MPI_BCAST(phi,nkx0*nky0*nkz0,MPI_DOUBLE_COMPLEX,0,MPI_COMM_HERM,ierr) 
  
!  IF (mype==0) WRITE(*,*) "phi",Real(phi(10,10,10))

 END SUBROUTINE get_phi_ae_mu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                           initialize_phi                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE initialize_phi
    !Define phi_denom here
    IMPLICIT NONE

    INTEGER :: i,j

    IF(.not.flr_on) THEN
      phi_denom=Ti0Te
    ELSE
      DO i=0,nkx0-1
        DO j=0,nky0-1
         phi_denom(i,j)=Ti0Te+(1.0-gamma0((kxgrid(i)**2+kygrid(j)**2)))
        END DO
      END DO
    END IF

  END SUBROUTINE initialize_phi

END MODULE field_solver


