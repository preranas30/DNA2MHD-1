!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 29/12/2012                                                                !!
!!                                cc_flr.f90                                 !!
!!                                                                           !!
!!  flr_effects                                                              !!
!!  -- get_J0                                                                !!
!!  -- finalize_flr                                                          !!
!!  -- GetMuWeightsAndKnots                                                  !!
!!  -- gauleg                                                                !!
!!  -- gaulag                                                                !!
!!                                                                     1.000 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                               flr_effects                                 !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE flr_effects
  USE par_mod
  USE aux_func 

  PUBLIC :: get_J0,finalize_flr,J0a

  REAL, ALLOCATABLE, DIMENSION(:,:) :: J0a

  PRIVATE


  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                  get_J0                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_J0
    IMPLICIT NONE

    REAL :: mumax,J0avg0
    INTEGER :: i,j,m
    INTEGER :: flr_handle
    INTEGER :: nmu0

    IF(.not.allocated(J0a)) ALLOCATE(J0a(0:nkx0-1,0:nky0-1))
    CALL get_io_number
    flr_handle=io_number
    IF(mype==0) OPEN(unit=flr_handle,file=trim(diagdir)//'/J0avg.dat',status='unknown')

    IF(flr_on) THEN

     IF(flr_version==1) THEN  !e(-b/2)
  
      DO i=0,nkx0-1
       DO j=0,nky0-1
        J0a(i,j)=e**(-0.5*kperp2(i,j))
        IF(mype==0.and.i==0) WRITE(flr_handle,*) kygrid(j),J0a(i,j)
       END DO
      END DO
    
     ELSE IF(flr_version==2) THEN !sqrt(b)

      DO i=0,nkx0-1
       DO j=0,nky0-1
        J0a(i,j)=sqrt(gamma0(kperp2(i,j)))
        IF(mype==0.and.i==0) WRITE(flr_handle,*) kygrid(j),J0a(i,j)
       END DO
      END DO

     ELSE
       STOP  "Error in get_J0!"
     END IF

    ELSE
      J0a(:,:)=1.0
      DO j=0,nky0-1
        IF(mype==0) WRITE(flr_handle,*) kygrid(j),J0a(0,j)
      END DO
    END IF

    IF(mype==0) CLOSE(flr_handle)

  END SUBROUTINE get_j0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                               finalize_flr                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE finalize_flr
    IMPLICIT NONE

    DEALLOCATE(J0a)

  END SUBROUTINE finalize_flr



END MODULE flr_effects
