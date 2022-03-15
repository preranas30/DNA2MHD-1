 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 22/03/2013                                                                !!
!!                         cc_hk_effects.f90                                 !!
!!                                                                           !!
!!  hankel terms                                                             !!
!!  -- get_hk                                                                !!
!!  -- finalize_hk                                                           !!
!!                                                                     1.000 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                               hk_effects                                  !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE hk_effects
  USE mpi
  USE par_mod
  USE aux_func 
  USE communications 

  PUBLIC :: get_hk, Gethankelgrid, hankel_transform, finalize_hk, J0_fac, I0a, I1a, J0zeros,&
            integral_hk, nearest_value, integral_v, get_phiavg, get_cfgamma, make_global_hk,&
            get_phiavg_hk, get_cfgamma_hk

  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: J0_fac
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: I0a
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: I1a

  PRIVATE


  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                  get_hk                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_hk
    IMPLICIT NONE

    INTEGER :: i,j,m,h
    INTEGER :: J0_fac_handle
    INTEGER :: kperp_handle
    INTEGER :: I0a_handle
    INTEGER :: I1a_handle
    INTEGER :: v_handle

    IF(.not.allocated(J0_fac)) ALLOCATE(J0_fac(0:nkx0-1,0:nky0-1,0:nh0-1))
    IF(.not.allocated(I0a)) ALLOCATE(I0a(0:nkx0-1,0:nky0-1,0:nh0-1))
    IF(.not.allocated(I1a)) ALLOCATE(I1a(0:nkx0-1,0:nky0-1,0:nh0-1))
 

    CALL get_io_number
    v_handle=io_number
    IF (mype == 0) THEN
        OPEN(unit=v_handle,file=trim(diagdir)//'/vgrid.dat',status='unknown')
        DO h=0,nh0-1
          WRITE (v_handle, "(3ES12.4)") delta_v(h),vgrid(h)
        ENDDO
        CLOSE(v_handle)
    ENDIF
   
    CALL get_io_number
    J0_fac_handle=io_number
    IF(mype==0) OPEN(unit=J0_fac_handle,file=trim(diagdir)//'/J0fac.dat',status='unknown')

    CALL get_io_number
    I0a_handle=io_number
    IF(mype==0) OPEN(unit=I0a_handle,file=trim(diagdir)//'/I0a.dat',status='unknown')

    CALL get_io_number
    I1a_handle=io_number
    IF(mype==0) OPEN(unit=I1a_handle,file=trim(diagdir)//'/I1a.dat',status='unknown')

    CALL get_io_number
    kperp_handle=io_number
    IF(mype==0) OPEN(unit=kperp_handle,file=trim(diagdir)//'/kperpgrid.dat',status='unknown')

    DO i=0,nkx0-1
        DO j=0,nky0-1
            IF(mype==0) WRITE(kperp_handle,*) kxgrid(i),kygrid(j),sqrt(kperp2(i,j))
        ENDDO
    ENDDO

  IF (.not.mu_integrated) THEN 
    IF(hankel) THEN
      DO i=0,nkx0-1
       DO j=0,nky0-1
         Do h=0,nh0-1
            J0_fac(i,j,h)=bessel0(sqrt(2*kperp2(i,j))*vgrid(h))
            I0a(i,j,h)= 0.5*pi**(-1)*e**(-0.5*kperp2(i,j)-0.25*hkgrid(h)**2)*&
                                        bessel_mod0(sqrt(kperp2(i,j))*hkgrid(h)/sqrt(2.))
            I1a(i,j,h)=-0.5*pi**(-1)*e**(-0.5*kperp2(i,j)-0.25*hkgrid(h)**2)*((-0.5*kperp2(i,j)**2 -0.25*hkgrid(h)**2)*&
                bessel_mod0(sqrt(kperp2(i,j))*hkgrid(h)/sqrt(2.)) +& 
                sqrt(kperp2(i,j))*hkgrid(h)/sqrt(2.)*bessel_mod1(sqrt(kperp2(i,j))*hkgrid(h)/sqrt(2.)))
            IF(mype==0) WRITE(J0_fac_handle,"(3I3,ES12.4)") i,j,h,J0_fac(i,j,h)
            IF(mype==0) WRITE(I0a_handle,"(3I3,ES12.4)") i,j,h,I0a(i,j,h)
            IF(mype==0) WRITE(I1a_handle,"(3I3,ES12.4)") i,j,h,I1a(i,j,h)
         END DO
        END DO
      END DO
    ELSE
      DO i=0,nkx0-1
       DO j=0,nky0-1
         DO h=0,nh0-1
           J0_fac(i,j,h)=bessel0(sqrt(2*kperp2(i,j)*vgrid(h)))
           I0a(i,j,h)=pi**(-1)*e**(-vgrid(h))*J0_fac(i,j,h)
           I1a(i,j,h)=vgrid(h)*I0a(i,j,h)
           IF(mype==0) WRITE(J0_fac_handle,"(3I3,ES12.4)") i,j,h,J0_fac(i,j,h)
           IF(mype==0) WRITE(I0a_handle,"(3I3,ES12.4)") i,j,h,I0a(i,j,h)
           IF(mype==0) WRITE(I1a_handle,"(3I3,ES12.4)") i,j,h,I1a(i,j,h)
         END DO
        END DO
      END DO
    ENDIF
  ELSE
     J0_fac(:,:,:)=1.0
     I0a(:,:,:)=1.0
     I1a(:,:,:)=1.0
     DO h=0,nh0-1
       IF(mype==0) WRITE(J0_fac_handle,"(ES12.4)") J0_fac(0,0,h)
       IF(mype==0) WRITE(I0a_handle,"(ES12.4)") I0a(0,0,h)
       IF(mype==0) WRITE(I1a_handle,"(ES12.4)") I1a(0,0,h)
     END DO
  END IF

  IF(mype==0) CLOSE(J0_fac_handle)
  IF(mype==0) CLOSE(kperp_handle)
  IF(mype==0) CLOSE(I0a_handle)
  IF(mype==0) CLOSE(I1a_handle)

  END SUBROUTINE get_hk
  
  SUBROUTINE finalize_hk
    IMPLICIT NONE

    DEALLOCATE(J0_fac)
    DEALLOCATE(I0a)
    DEALLOCATE(I1a)

  END SUBROUTINE finalize_hk


  SUBROUTINE Gethankelgrid(hk_max,v_max,delta_hk,delta_v,hkgrid,vgrid,m1,m2,Tmn)
    IMPLICIT NONE

    REAL, INTENT(IN) :: hk_max 
    REAL, DIMENSION (0:nh0) :: J0zeros
    REAL, DIMENSION (0:nh0-1), INTENT(OUT) :: m1, m2 
    REAL, DIMENSION (0:nh0-1), INTENT(OUT) :: delta_hk, delta_v, hkgrid, vgrid
    !REAL, DIMENSION (lh1:lh2) :: test_loc
    !REAL, DIMENSION (0:nh0-1) :: test_glob
    REAL, DIMENSION (0:nh0-1,0:nh0-1), INTENT(OUT) :: Tmn
    REAL, INTENT(OUT) :: v_max 
    REAL, DIMENSION (0:nh0-1,0:nh0-1) :: Jn, Jm
    INTEGER :: h,n,m,ierr
    INTEGER :: J0_handle, hk_handle

    CALL get_io_number
    J0_handle = io_number
    IF (mype==0) THEN
      OPEN(J0_handle, file='J0zeros.txt', status = 'unknown', action='read' )
      DO h=0,nh0
        READ (J0_handle, *) J0zeros(h)
      ENDDO
      CLOSE(J0_handle)

      v_max = J0zeros(nh0)/(hk_max)
      hkgrid(0:nh0-1) = J0zeros(0:nh0-1)*hk_max/J0zeros(nh0)
      vgrid(0:nh0-1) = J0zeros(0:nh0-1)/(hk_max)

    ENDIF

    !Broadcast hkgrid, vgrid, v_max,and zeros of Bessel function  to all processors
    CALL MPI_BCAST(hkgrid,nh0,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
    CALL MPI_BCAST(vgrid,nh0,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
    CALL MPI_BCAST(v_max,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
    CALL MPI_BCAST(J0zeros,nh0+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 

    !So far a very simple discretization
    DO h=0,nh0-1
         IF (h==0) THEN
             delta_hk(h) = hkgrid(h)
             delta_v(h) = vgrid(h)
         ELSE 
             delta_hk(h) = hkgrid(h) - hkgrid(h-1)
             delta_v(h) = vgrid(h) - vgrid(h-1)
         END IF
    END DO
 
    !Routines for the Inverse Hankel transform, so far only  mype_hank == 0 should do it

 IF (mype_hank == 0) THEN
     
     DO h=0,nh0-1
         Jn(:,h) = J0zeros(h)
         Jm(h,:) = J0zeros(h)
     END DO     

     DO n=0,nh0-1
       DO  m=0,nh0-1
         Tmn(n,m) = (2/J0zeros(nh0))*bessel0(Jn(n,m)*Jm(n,m)/J0zeros(nh0))/(abs(bessel1(Jn(n,m)))*abs(bessel1(Jm(n,m))))
       ENDDO
     ENDDO 

     DO h=0,nh0-1
        m1(h) = (abs(bessel1(J0zeros(h)))/v_max) 
        m2(h) = m1(h)*v_max/hk_max
     ENDDO
 END IF

 END SUBROUTINE Gethankelgrid

 SUBROUTINE hankel_transform(f_in,hk_trans)
    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: hk_trans
    COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2), INTENT(INOUT) :: f_in
    COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,0:nh0-1,ls1:ls2) :: f_glob
    COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,0:nh0-1,ls1:ls2) :: f_glob_int
    COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2) :: f_out
    COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,0:nh0-1,ls1:ls2) :: F, F2 
    INTEGER :: i,j,k,l,s,ierr
    INTEGER :: hk_handle
  
   f_glob = CMPLX(0.0,0.0) 
   f_glob_int = CMPLX(0.0,0.0) 
   f_out  = CMPLX(0.0,0.0) 
   CALL MPI_GATHER(f_in,size(f_in),MPI_DOUBLE_COMPLEX,f_glob,size(f_in),MPI_DOUBLE_COMPLEX,0,MPI_COMM_HANK,ierr)

   DO i=0,nkx0-1
    DO j=0,nky0-1
     DO k=lkz1,lkz2
      DO l=lv1,lv2
       DO s=ls1,ls2
        IF (mype_hank==0) THEN
           IF (hk_trans) THEN
           !Obtain inverse hankel transform
              F2(i,j,k,l,:,s) = f_glob(i,j,k,l,:,s)/m2_hk
              F(i,j,k,l,:,s) = matmul(T_hkv,F2(i,j,k,l,:,s))
              f_glob_int(i,j,k,l,:,s) = F(i,j,k,l,:,s)*m1_v
           ELSE
              F(i,j,k,l,:,s) = f_glob(i,j,k,l,:,s)/m1_v
              F2(i,j,k,l,:,s) = matmul(T_hkv,F(i,j,k,l,:,s))
              f_glob_int(i,j,k,l,:,s) = F2(i,j,k,l,:,s)*m2_hk
           ENDIF
        ENDIF
      END DO
     END DO
   END DO
  END DO
 END DO

 CALL MPI_SCATTER(f_glob_int,size(f_out),MPI_DOUBLE_COMPLEX,f_out,size(f_out),MPI_DOUBLE_COMPLEX,0,MPI_COMM_HANK,ierr)
 f_in = f_out

!     CALL get_io_number
!    hk_handle = io_number
!    IF (hk_trans) THEN
!    IF (mype == 0) THEN
!        OPEN(unit=hk_handle,file=trim(diagdir)//'/v.dat',status='unknown')
!        DO h=lh1,lh2
!            WRITE (hk_handle, "(2ES12.4)") f_in(h)
!        ENDDO
!    ENDIF
!    ELSE
!    IF (mype == 0) THEN
!        OPEN(unit=hk_handle,file=trim(diagdir)//'/h.dat',status='unknown')
!        DO h=lh1,lh2
!            WRITE (hk_handle, "(2ES12.4)") f_in(h)
!        ENDDO
!    ENDIF
!    ENDIF 

  END SUBROUTINE hankel_transform

 SUBROUTINE make_global_hk(f_in,f_glob)
    IMPLICIT NONE

    COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2), INTENT(IN) :: f_in
    COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,0:nh0-1,ls1:ls2), INTENT(OUT) :: f_glob
    INTEGER :: i,j,k,l,s,ierr
    INTEGER :: hk_handle
  
    f_glob = CMPLX(0.0,0.0) 
    DO i=0,nkx0-1
     DO j=0,nky0-1
      DO k=lkz1,lkz2
       DO l=lv1,lv2
        DO s=ls1,ls2
        CALL MPI_GATHER(f_in(i,j,k,l,0,s),lh0,MPI_DOUBLE_PRECISION,f_glob(i,j,k,l,lh1,s), &
                        lh0,MPI_DOUBLE_PRECISION,0,MPI_COMM_HANK,ierr)
        END DO
       END DO
      END DO
     END DO
    END DO
 
  END SUBROUTINE make_global_hk

  SUBROUTINE nearest_value(value,idx)
   
   IMPLICIT NONE
   REAL, INTENT(IN) :: value
   INTEGER, DIMENSION(1),INTENT(OUT) :: idx
   REAL,DIMENSION(0:nh0-1)  :: diff
    
   IF (mype==0) THEN
      diff = abs(hkgrid - value)
      idx = minloc(diff)
      idx = idx -1
   ENDIF

  END SUBROUTINE nearest_value


  SUBROUTINE integral_hk(f_in,prefac,f_out)
    !Not the best efficient routine 
    IMPLICIT NONE

    COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2), INTENT(IN) :: f_in
    REAL, DIMENSION (0:nkx0-1,0:nky0-1,lh1:lh2), INTENT(IN) :: prefac
    COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,ls1:ls2), INTENT(OUT) :: f_out
    INTEGER :: i,j,h

    f_out = cmplx(0.0,0.0)
     DO i=0,nkx0-1
      DO j =0,nky0-1
       DO h=lh1,lh2
          f_out(i,j,:,:,:) = prefac(i,j,h)*hkgrid(h)*delta_hk(h)*f_in(i,j,:,:,h,:) + f_out(i,j,:,:,:)
       ENDDO
      ENDDO
     ENDDO

    Call my_complex_sum_hank(f_out,size(f_out))
 
  END SUBROUTINE integral_hk

  SUBROUTINE get_phiavg_hk(phi_in,phi_avg)

   IMPLICIT NONE
   COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2), INTENT(IN) :: phi_in
   COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2), INTENT(OUT) :: phi_avg
   INTEGER :: i,j,k,l,h,s

   phi_avg = CMPLX(0.0,0.0)
    DO  k=lkz1,lkz2
       DO l=lv1,lv2
         DO h =lh1,lh2
            DO s=ls1,ls2
               phi_avg(:,:,k,l,h,s) = I0a(:,:,h)*phi_in(:,:,k)
            ENDDO
         ENDDO   
       ENDDO
    ENDDO

  END SUBROUTINE get_phiavg_hk

  SUBROUTINE get_cfgamma_hk(f_in,phi_in,cfgamma1,cfgamma2)

  IMPLICIT NONE
  COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2), INTENT(IN) :: f_in
  COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2), INTENT(IN) :: phi_in
  COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2), INTENT(OUT) :: cfgamma1,cfgamma2
  COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,0:nh0-1,ls1:ls2) :: phi_avg_glob
  COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2) :: phi_avg
  REAL, DIMENSION(0:nkx0-1,0:nky0-1) :: kperp
  INTEGER :: i,j,k,l,h, ierr
  INTEGER, DIMENSION(1) :: p_kperp
  
  kperp = sqrt(2*kperp2)  
  phi_avg = CMPLX(0.0,0.0)
  phi_avg_glob = CMPLX(0.0,0.0)
  cfgamma1 = CMPLX(0.0,0.0)
  cfgamma2 = CMPLX(0.0,0.0)
    DO i=0,nkx0-1
      DO j = 0,nky0-1
        DO k=0,nkz0 -1
          CALL nearest_value(kperp(i,j),p_kperp)
          IF (mype == 0) THEN  
             phi_avg_glob(i,j,k,0,p_kperp(1),:) = phi(i,j,k)
          END IF
        END DO
      END DO
    END DO 

    CALL MPI_SCATTER(phi_avg_glob,size(phi_avg),MPI_DOUBLE_COMPLEX,phi_avg,size(phi_avg),MPI_DOUBLE_COMPLEX,0,MPI_COMM_HANK,ierr)
    IF(lv1.le.0.and.lv2.ge.0) cfgamma2(:,:,:,0,:,:) = pi**(5.0/4.0)*conjg(phi_avg(:,:,:,0,:,:))
    cfgamma1 = pi**(5.0/2.0)*conjg(f_in) 

  END SUBROUTINE get_cfgamma_hk

!Routines for the mu dimension (althoug I called it vperp)
  SUBROUTINE integral_v(f_in,prefac,f_out)
   
    IMPLICIT NONE

    COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2), INTENT(IN) :: f_in
    REAL, DIMENSION (0:nkx0-1,0:nky0-1,lh1:lh2), INTENT(IN) :: prefac
    COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,ls1:ls2), INTENT(OUT) :: f_out
    INTEGER :: i,j,h

    f_out = cmplx(0.0,0.0)

     DO i=0,nkx0-1
      DO j =0,nky0-1
       DO h=lh1,lh2
          f_out(i,j,:,:,:) = prefac(i,j,h)*delta_v(h)*f_in(i,j,:,:,h,:) + f_out(i,j,:,:,:)
       ENDDO
      ENDDO
     ENDDO
   
   Call my_complex_sum_hank(f_out,size(f_out))

  END SUBROUTINE integral_v

  SUBROUTINE get_phiavg(phi_in,phi_avg)

   IMPLICIT NONE
   COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2), INTENT(IN) :: phi_in
   COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2), INTENT(OUT) :: phi_avg
   INTEGER :: i,j,k,l,h,s


   phi_avg = CMPLX(0.0,0.0)
    DO  k=lkz1,lkz2
       DO l=lv1,lv2
         DO h =lh1,lh2
            DO s=ls1,ls2
               phi_avg(:,:,k,l,h,s) = J0_fac(:,:,h)*phi_in(:,:,k)
            ENDDO
         ENDDO   
       ENDDO
    ENDDO

  END SUBROUTINE get_phiavg


  SUBROUTINE get_cfgamma(f_in,phi_in,cfgamma1,cfgamma2)

  IMPLICIT NONE
  COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2), INTENT(IN) :: f_in
  COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2), INTENT(IN) :: phi_in
  COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2), INTENT(OUT) :: cfgamma1,cfgamma2
  COMPLEX, DIMENSION (0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2) :: phi_avg
  INTEGER :: l,h

  phi_avg = CMPLX(0.0,0.0)
  cfgamma1 = CMPLX(0.0,0.0)
  cfgamma2 = CMPLX(0.0,0.0)
  call get_phiavg(phi_in, phi_avg)

    DO l=lv1,lv2
      DO h=lh1,lh2
         IF(lv1.le.0.and.lv2.ge.0) cfgamma2(:,:,:,l,h,:) = pi**(5.0/4.0)*conjg(phi_avg(:,:,:,0,h,:))
         cfgamma1(:,:,:,l,h,:)  = pi**(5.0/2.0)*conjg(f_in(:,:,:,l,h,:))*e**(vgrid(h)) 
      ENDDO 
    ENDDO

  END SUBROUTINE get_cfgamma



! SUBROUTINE hankel_grid(hk_max,v_max,delta_hk,delta_v,hk_grid,v_grid,m1,m2,Tmn)
!    IMPLICIT NONE
!
!    REAL, INTENT(IN) :: hk_max 
!    REAL, DIMENSION (lh1:lh2+1) :: J0zeros
!    REAL, DIMENSION (lh1:lh2), INTENT(OUT) :: m1, m2 
!    REAL, DIMENSION (lh1:lh2), INTENT(OUT) :: delta_hk, delta_v, hk_grid, v_grid
!    REAL, DIMENSION (lh1:lh2,lh1:lh2), INTENT(OUT) :: Tmn
!    REAL, INTENT(OUT) :: v_max 
!    REAL, DIMENSION (lh1:lh2,lh1:lh2) :: Jn, Jm
!    INTEGER :: h,n,m
!    INTEGER :: J0_handle, hk_handle
!
!    CALL get_io_number
!    J0_handle = io_number
!    OPEN(J0_handle, file='J0zeros.txt', status = 'unknown', action='read' )
!    DO h=lh1,lh2+1
!        READ (J0_handle, *) J0zeros(h)
!    ENDDO
!    CLOSE(J0_handle)
!
!    v_max = J0zeros(nh0)/(hk_max)
!    hk_grid(lh1:lh2) = J0zeros(lh1:lh2)*hk_max/J0zeros(nh0)
!    v_grid(lh1:lh2) = J0zeros(lh1:lh2)/(hk_max)
!
!
!    DO h=lh1,lh2
!         IF (h==lh1) THEN
!             delta_hk(h) = hk_grid(h)
!             delta_v(h) = v_grid(h)
!         ELSE 
!             delta_hk(h) = hk_grid(h) - hk_grid(h-1)
!             delta_v(h) = v_grid(h) - v_grid(h-1)
!         END IF
!    END DO
! 
!    DO h=lh1,lh2
!        Jn(:,h) = J0zeros(h)
!        Jm(h,:) = J0zeros(h)
!    END DO     
!
!
!    DO n=lh1,lh2
!      DO  m=lh1,lh2
!        Tmn(n,m) = (2/J0zeros(nh0))*bessel0(Jn(n,m)*Jm(n,m)/J0zeros(nh0))/(abs(bessel1(Jn(n,m)))*abs(bessel1(Jm(n,m))))
!      ENDDO
!    ENDDO 
!
!   DO h=lh1,lh2
!      m1(h) = (abs(bessel1(J0zeros(h)))/v_max) 
!      m2(h) = m1(h)*v_max/hk_max
!   ENDDO
!
!
!
!  END SUBROUTINE hankel_grid


! SUBROUTINE hankel_transform(f_in,hk_trans)
!    IMPLICIT NONE
!
!    LOGICAL, INTENT(IN) :: hk_trans
!    REAL, DIMENSION (lh1:lh2), INTENT(INOUT) :: f_in
!    REAL, DIMENSION (lh1:lh2) :: F, F2
!    INTEGER :: h
!    INTEGER :: hk_handle
!
! 
!   !Obtain inverse hankel transform
!   IF (hk_trans) THEN
!       F2 = f_in/m2_hk
!       F = matmul(T_hkv,F2)
!       f_in = F*m1_v
!   ELSE
!   !Obtain Hankel transform
!      F = f_in/m1_v
!     F2 = matmul(T_hkv,F)
!      f_in = F2*m2_hk
!   ENDIF
!
!    CALL get_io_number
!    hk_handle = io_number
!    IF (hk_trans) THEN
!    IF (mype == 0) THEN
!        OPEN(unit=hk_handle,file=trim(diagdir)//'/vgrid.dat',status='unknown')
!       WRITE (hk_handle, "(1ES12.4)") vmax
!        DO h=lh1,lh2
!            WRITE (hk_handle, "(3ES12.4)") vgrid(h),delta_v(h),f_in(h)
!        ENDDO
!    ENDIF
!    ELSE
!    IF (mype == 0) THEN
!        OPEN(unit=hk_handle,file=trim(diagdir)//'/hkgrid.dat',status='unknown')
!        WRITE (hk_handle, "(1ES12.4)") hkmax
!        DO h=lh1,lh2
!            WRITE (hk_handle, "(3ES12.4)") hkgrid(h),delta_hk(h),f_in(h)
!        ENDDO
!    ENDIF
!    ENDIF 
!
!  END SUBROUTINE hankel_transform
!
!
! SUBROUTINE hankel_transform1(f_in,f_hk,f_inv,Tmn,m1,m2,hkgrid,vgrid)
!    IMPLICIT NONE
!
!    REAL, DIMENSION (lh1:lh2), INTENT(IN) :: f_in
!    REAL, DIMENSION (lh1:lh2), INTENT(IN) :: m1, m2 
!    REAL, DIMENSION (lh1:lh2,lh1:lh2), INTENT(IN) :: Tmn
!    REAL, DIMENSION (lh1:lh2),INTENT(IN) :: hkgrid, vgrid
!    REAL, DIMENSION (lh1:lh2), INTENT(OUT) :: f_hk, f_inv
!    REAL, DIMENSION (lh1:lh2) :: F, F2, Fretrieved
!    INTEGER :: h
!    INTEGER :: hk_handle
!
! 
!   !Obtain hankel transform
!   F2 = f_in/m2
!   F = matmul(Tmn,F2)
!   f_inv = F*m1
!
!    !Obtain inverse Hankel transform
!
!   Fretrieved  = matmul(Tmn,F)
!   f_hk = Fretrieved*m2
!
!    CALL get_io_number
!    hk_handle = io_number
!    IF (mype == 0) THEN
!        OPEN(unit=hk_handle,file=trim(diagdir)//'/hkgrid.dat',status='unknown')
!        DO h=lh1,lh2
!            WRITE (hk_handle, "(5ES12.4)") hkgrid(h),vgrid(h),f_in(h),f_hk(h),f_inv(h)
!        ENDDO
!    ENDIF
!
!  END SUBROUTINE hankel_transform1

!  SUBROUTINE J0_zeros
!    IMPLICIT NONE
!
!    REAL :: vmax, hkmax 
!    REAL, DIMENSION (lh1:lh2+1) :: J0zeros
!    REAL, DIMENSION (lh1:lh2) :: m1, m2 
!    REAL, DIMENSION (lh1:lh2) :: hkgrid, vgrid
!    REAL, DIMENSION (lh1:lh2) :: f_in, F, F2, Fretrieved, f_inv, f_hk
!    REAL, DIMENSION (lh1:lh2,lh1:lh2) :: Jn, Jm, Tmn
!    INTEGER :: h,n,m
!    INTEGER :: J0_handle, hk_handle
!
!    hkmax = 3.0
!    CALL get_io_number
!    J0_handle = io_number
!   IF (mype == 0) THEN
!        OPEN(J0_handle, file='J0zeros.txt', status = 'unknown', action='read' )
!        DO h=lh1,lh2+1
!            READ (J0_handle, *) J0zeros(h)
!        ENDDO
!
!    CLOSE(J0_handle)
!    ENDIF
!
!    vmax = J0zeros(nh0)/(2*pi*hkmax)
!    hkgrid(lh1:lh2) = J0zeros(lh1:lh2)*hkmax/J0zeros(nh0)
!    vgrid(lh1:lh2) = J0zeros(lh1:lh2)/(2*pi*hkmax)
! 
!    !Test function
!    f_in(lh1:lh2) = 2*pi*e**(-2.*pi**2*hkgrid(lh1:lh2)**2)
!
!   
!    DO h=lh1,lh2
!        Jn(:,h) = J0zeros(h)
!        Jm(h,:) = J0zeros(h)
!    END DO     
!
!
!   DO n=lh1,lh2
!      DO  m=lh1,lh2
!        Tmn(n,m) = (2/J0zeros(nh0))*bessel0(Jn(n,m)*Jm(n,m)/J0zeros(nh0))/(abs(bessel1(Jn(n,m)))*abs(bessel1(Jm(n,m))))
!      ENDDO
!    ENDDO 
!
!   DO h=lh1,lh2
!      m1(h) = (abs(bessel1(J0zeros(h)))/hkmax) 
!      m2(h) = m1(h)*hkmax/vmax
!   ENDDO
!
!   !Obtain hankel transform
!   F = f_in/m1
!   F2 = matmul(Tmn,F)
!   f_inv = F2*m2
!
!    !Obtain inverse Hankel transform
!
!   Fretrieved  = matmul(Tmn,F2)
!   f_hk = Fretrieved*m1
!
!    CALL get_io_number
!    hk_handle = io_number
!    IF (mype == 0) THEN
!        OPEN(unit=hk_handle,file=trim(diagdir)//'/hkgrid.dat',status='unknown')
!        DO h=lh1,lh2
!            WRITE (hk_handle, "(6ES12.4)") J0zeros(h),hkgrid(h),vgrid(h),f_in(h),f_hk(h),f_inv(h)
!        ENDDO
!    ENDIF
!
!  END SUBROUTINE J0_zeros

END MODULE hk_effects


