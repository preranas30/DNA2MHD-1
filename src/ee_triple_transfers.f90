!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 04/02/2013                                                                !!
!!                           ee_triple_transfers.f90                         !!
!!                                                                           !!
!!  triple_transfers                                                         !!
!!  -- get_triple_transfers                                                  !!
!!  -- make_shell_transfers                                                  !!
!!  -- shell_filter                                                          !!   
!!                                                                     1.000 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                               triple_transfers                            !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE triple_transfers
  USE par_mod
  USE mpi
  IMPLICIT NONE

  PUBLIC :: get_triple_transfers

  PRIVATE

  INTEGER   :: shell_type
  INTEGER   :: kinit
  INTEGER   :: Nshell
  REAL, DIMENSION(0:100)    :: kgiven 
  NAMELIST /triple_transfers_para/ shell_type, kinit, Nshell, kgiven
  
  INTEGER, ALLOCATABLE, DIMENSION(:)  :: sss_modes !! Number of modes in a shell
  REAL, ALLOCATABLE, DIMENSION(:)     :: kshell    !! Shell lower wavenumber
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: sss_pgg   !! Triple shell transfers 

  
  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                             get_triple_transfers                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_triple_transfers
    USE field_solver  
    IMPLICIT NONE
    
    INTEGER :: s1, s2, s3
    INTEGER :: ierr
    REAL :: gshell, kfirst
    

    IF(np_hank.gt.1) STOP "get_triple_transfers not yet implemented for np_hank.gt.1"
    IF(np_spec.gt.1) STOP "get_triple_transfers not yet implemented for np_spec.gt.1"
    IF(np_kz.gt.1) STOP "get_triple_transfers not yet implemented for np_kz.gt.1"
    
    !! READ PARAMETERS
    !!!!!!!!!!!!!!!!!!!!
    shell_type=4
    kinit=1
    Nshell=5
    kgiven=0.0
    CALL get_io_number
    OPEN(io_number, file='parameters', status = 'unknown')
    READ(io_number, nml=triple_transfers_para, iostat = ierr)
    CLOSE(io_number)


    !! DEFINE SHELL BORDERS
    !!!!!!!!!!!!!!!!!!!!!!!!!
    ALLOCATE(kshell(1:Nshell+1))  
    SELECT CASE(kinit)
    CASE(1)
      gshell=2.0**(0.2)
	  kfirst=0.1
      kshell(1)=0.0
	  DO s1 = 2,Nshell
        kshell(s1)=kfirst*gshell**(s1-2)
      END DO
      kshell(Nshell+1) = 10000000.0
    CASE(2)
      gshell=2.0**(0.2)
	  kfirst=kgiven(1)
      kshell(1)=0.0
	  DO s1 = 2,Nshell
        kshell(s1)=kfirst*gshell**(s1-2)
      END DO
      kshell(Nshell+1) = 10000000.0
    END SELECT
    IF(mype==0)  WRITE(*,*) kshell 
    
    
    !! COMPUTE TRANSFERS
    !!!!!!!!!!!!!!!!!!!!!! 

    ALLOCATE(sss_pgg(Nshell, Nshell, Nshell))  
    CALL get_phi(g_1)
    CALL make_shell_transfers(phi, g_1, sss_pgg)    
        
        
    !! WRITE OUTPUT
    !!!!!!!!!!!!!!!!!    
    CALL get_io_number
    IF (mype==0) THEN
      OPEN(unit=io_number,file=trim(diagdir)//'/sss.txt',status='unknown')
      WRITE (io_number,*) "kfirst=",kfirst
      WRITE (io_number,*) "kmax=",kshell(Nshell)
      DO s1=1,Nshell
        DO s2=1,Nshell
          WRITE (io_number,*) "Transfer spectra T(L3) for shells L1, L2"
          WRITE (io_number,'(2I5)') s1, s2
          DO s3=1,Nshell
            WRITE (io_number,'(I5,E16.8)') s3, sss_pgg(s1,s2,s3)
          END DO
        END DO  
      END DO
      CLOSE(io_number)
    END IF
  
    DEALLOCATE(kshell,sss_pgg)
  
  END SUBROUTINE get_triple_transfers
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                              make_shell_transfers                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE make_shell_transfers(phi_in, g_in, sss)
    USE nonlinearity
    IMPLICIT NONE
    
    COMPLEX, INTENT(IN) :: phi_in(0:nkx0-1,0:nky0-1,lkz1:lkz2)
    COMPLEX, INTENT(IN) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
    REAL, DIMENSION(Nshell,Nshell,Nshell), INTENT(OUT) :: sss
  
    COMPLEX, ALLOCATABLE, DIMENSION(:,:,:)       :: Zf !! helper filtered
    COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: Yf !! giver filtered
    COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: Xf !! reciver filtered
    COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: NL !! Nonlinear term  
    INTEGER :: s1, s2, s3, n, i, j, k, ierr
    REAL :: local

    IF(np_hank.gt.1) STOP "make_shell_transfers not yet implemented for np_hank.gt.1"
    IF(np_spec.gt.1) STOP "make_shell_transfers not yet implemented for np_spec.gt.1"
    IF(np_kz.gt.1) STOP "make_shell_transfers not yet implemented for np_kz.gt.1"
     
    ALLOCATE(Zf(0:nkx0-1,0:nky0-1,0:nkz0-1)) 
    ALLOCATE(Yf(0:nkx0-1,0:nky0-1,0:nkz0-1,lv1:lv2,lh1:lh2,ls1:ls2))  
    ALLOCATE(Xf(0:nkx0-1,0:nky0-1,0:nkz0-1,lv1:lv2,lh1:lh2,ls1:ls2)) 
    ALLOCATE(NL(0:nkx0-1,0:nky0-1,0:nkz0-1,lv1:lv2,lh1:lh2,ls1:ls2)) 
        
    sss=0.0
    
    DO s1=1,Nshell
      IF(mype==0) WRITE(*,*) 'loop 1: radial index is ', s1
      CALL shell_filter(phi,kshell(s1),kshell(s1+1),Zf)
      DO s2=1,Nshell
        IF(mype==0) WRITE(*,*) 'loop 2: radial index is ', s2
        DO n=lv1,lv2
          CALL shell_filter(g_in(:,:,:,n,0,0),kshell(s2),kshell(s2+1),Yf(:,:,:,n,0,0))
        END DO  
        NL=cmplx(0.0,0.0)
        CALL get_rhs_nl(Yf,Zf,NL)
        IF(mype==0) WRITE(*,*) 'Nonlinear term done'
        DO s3=1,Nshell
          IF(mype==0) WRITE(*,*) 'loop 3: radial index is ', s3
          DO n=lv1,lv2
            CALL shell_filter(g_in(:,:,:,n,0,0),kshell(s3),kshell(s3+1),Xf(:,:,:,n,0,0))
          END DO  
          !!
          local=0.0
          DO k=0,nkz0-1
            DO j=0,nky0-1
              local=local+0.5*REAL(pi**0.5*SUM(NL(0,j,k,:,0,0)*CONJG(Xf(0,j,k,:,0,0))))
              DO i=1,nkx0-1
                local=local+REAL(pi**0.5*SUM(NL(i,j,k,:,0,0)*CONJG(Xf(i,j,k,:,0,0))))
              END DO 
            END DO
          END DO  
          CALL MPI_ALLREDUCE(local, sss(s1,s2,s3) ,1 , &
                             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
        END DO
      END DO  
    END DO 
  
    DEALLOCATE(Zf) 
    DEALLOCATE(Yf)  
    DEALLOCATE(Xf) 
    DEALLOCATE(NL) 
      
  END SUBROUTINE make_shell_transfers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                shell_filter                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE shell_filter(phi, kinf, ksup, phi_f)
    IMPLICIT NONE

    COMPLEX, DIMENSION(0:nkx0-1,0:nky0-1,0:nkz0-1), INTENT(IN)  :: phi   
    COMPLEX, DIMENSION(0:nkx0-1,0:nky0-1,0:nkz0-1), INTENT(OUT) :: phi_f 
    REAL, INTENT(IN) :: kinf, ksup

    INTEGER :: i, j, k
    REAL :: knorm
    
    phi_f=0.0 

    DO i=0,nkx0-1
      DO j=0,nky0-1
        DO k=0,nkz0-1
          !!
          SELECT CASE(shell_type)
          CASE(1)
            knorm=SQRT(kxgrid(i)**2)
          CASE(2)
            knorm=SQRT(kygrid(j)**2)
          CASE(3)
            knorm=SQRT(kzgrid(k)**2)
           CASE(4)
            knorm=SQRT(kxgrid(i)**2+kygrid(j)**2)
           CASE(5)
            knorm=SQRT(kxgrid(i)**2+kygrid(j)**2+kzgrid(k)**2)       
          END SELECT
          !!
          IF (kinf<knorm .AND. knorm <= ksup) phi_f(i,j,k)=phi(i,j,k)
        END DO
      END DO
    END DO
 
    IF(mype==0)  WRITE (*,*) 'Filter done!'

  END SUBROUTINE shell_filter


END MODULE triple_transfers







