! module containing all the FGS group simulations routines
! dynamic procedure routine
! Free energy spectral analysis (kx,ky,kperp ...)
! Free energy consevation analysis vs time
! Transfer function
Module Gyro_LES

  USE par_mod
  USE mpi
  USE nonlinearity
  USE field_solver, only: get_phi
  !USE diagnostics, only: sum3d_real
  USE flr_effects, only: J0a

    Implicit None
    Public :: initialize_GyroLES, exec_GyroLES, exec_GyroLES_herm, exec_GyroLES_z, Correlation, finalize_GyroLES,&
                initialize_corr, finalize_corr

    PRIVATE
    
    Character(Len=8) :: filestat='replace', filepos='rewind'
    Complex, Dimension(:,:,:,:,:,:), Allocatable  :: filter_g_1
    Complex, Dimension(:,:,:), Allocatable  :: filter_phi
    Complex, Dimension(:,:,:,:,:,:), Allocatable :: nlt
    Complex, Dimension(:,:,:,:,:,:), Allocatable :: filter_nlt    
    Complex, Dimension(:,:,:,:,:,:), Allocatable :: nlt_filter_g_1 
    Complex, Dimension(:,:,:,:,:,:), Allocatable:: filter_nlt_filter_g_1  
    Complex, Dimension(:,:,:,:,:,:), Allocatable  :: t_fgs
    complex, dimension(:,:,:,:,:,:), allocatable :: m_x, m_y, m_z
    real, dimension(:,:,:,:,:,:), allocatable :: tmp_submod
    complex, dimension(:,:,:,:,:,:), allocatable :: m_x_mod, m_y_mod
    real, dimension(:,:,:), allocatable :: tmp_3D_dp
    real, dimension(:,:,:), allocatable :: energy3d
    
    Integer  :: GyroLES_FILE, CORRELATION_FILE


Contains

subroutine initialize_corr
  
    implicit none

    IF (mype==0) then
        CALL get_io_number
        CORRELATION_FILE=io_number
        open(CORRELATION_FILE, file=trim(diagdir)//'/Correlation.dat', form='formatted', &
        &status=filestat, position=filepos)
          write(CORRELATION_FILE,*) "#  1. time"
          write(CORRELATION_FILE,*) "#  2. Subgrid dissipation"
          write(CORRELATION_FILE,*) "#  3. Correlation x model" 
          write(CORRELATION_FILE,*) "#  4. Correlation y model"
          write(CORRELATION_FILE,*) "#  5. Correlation z model"
    ENDIF
  
end subroutine initialize_corr

subroutine initialize_GyroLES
  
    implicit none

    If (mype==0) then
        CALL get_io_number
        GyroLES_FILE=io_number
        open(GyroLES_FILE, file=trim(diagdir)//'/GyroLES.dat', form='formatted', &
        &status=filestat, position=filepos)
        IF (Gyroherm) THEN
          write(GyroLES_FILE,*) "#  1. time"
          write(GyroLES_FILE,*) "#  2. Hermite"
          write(GyroLES_FILE,*) "#  3. hyp_x amplitude" 
          write(GyroLES_FILE,*) "#  4. hyp_y amplitude"
        ELSE
          write(GyroLES_FILE,*) "#  1. time"
          write(GyroLES_FILE,*) "#  2. hyp_x amplitude" 
          write(GyroLES_FILE,*) "#  3. hyp_y amplitude"
        ENDIF
    endif
  
end subroutine initialize_GyroLES

subroutine exec_GyroLES_herm(g_in)

    complex,dimension(0:nkx0-1,0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2),intent(in) :: g_in
    real,dimension(lv1:lv2) :: mymy, mxmx, mxmy, Tmx, Tmy, c_x, c_y
    real,dimension(0:nv0-1) :: hyp_x_herm_glob, hyp_y_herm_glob
    real :: eps_alpha_x, eps_alpha_y, tgs, mmod
    Integer :: i, j, l, k,ierr

    Allocate(filter_g_1(0:nkx0-1,0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))
    Allocate(filter_phi(0:nkx0-1,0:nky0-1, lkz1:lkz2))
    Allocate(nlt(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2)) 
    Allocate(filter_nlt(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))     
    Allocate(nlt_filter_g_1(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))  
    Allocate(filter_nlt_filter_g_1(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))   
    Allocate(t_fgs(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))
    Allocate(m_x(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2)) 
    Allocate(m_x_mod(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2)) 
    Allocate(m_y(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2)) 
    Allocate(m_y_mod(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2)) 
    allocate(tmp_submod(0:nkx0-1, 0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
    allocate(tmp_3D_dp(0:nkx0-1, 0:nky0-1,lkz1:lkz2))
    allocate(energy3d(0:nkx0-1, 0:nky0-1,lkz1:lkz2))

    !Apply the filter to g_1, compute the filtered fields
    Call get_phi(g_in)
    Call filter(g_in, filter_g_1)

    ! compute the T sub grid term
    nlt = cmplx(0.0,0.0)
    Call get_rhs_nl(g_in, phi, nlt)
    Call filter(nlt, filter_nlt)
    
    nlt_filter_g_1 = cmplx(0.0,0.0)
    Call filter(g_in, filter_g_1)
    Call filter_in_phi(phi, filter_phi)
    Call get_rhs_nl(filter_g_1, filter_phi, nlt_filter_g_1)
    Call filter(nlt_filter_g_1, filter_nlt_filter_g_1)     
    t_fgs = filter_nlt - filter_nlt_filter_g_1

!Some test to calculate the sub-grid contribution
    tmp_3d_dp = cmplx(0.0,0.0)
    energy3d =  0.0
    tgs = 0.0
  
    tmp_3d_dp(:,:,:)=real(sqrt(pi)*sum(conjg(g_in(:,:,:,:,0,0))*t_fgs(:,:,:,:,0,0),4))
    if(mype==0) then
      do k=0,nkz0-1
        tmp_3d_dp(:,:,k)=tmp_3d_dp(:,:,k)+real(pi**(0.25)*j0a(:,:)*conjg(phi(:,:,k))*t_fgs(:,:,k,0,0,0))
      end do
    end if
    call mpi_allreduce(tmp_3d_dp,energy3d,nkx0*nky0*nkz0 &
      ,mpi_double_precision,mpi_sum,mpi_comm_herm,ierr)
    call sum3d_real(energy3d,tgs)
!Finish test


    ! compute hyper-diff as models, with amplitude set to 1
    
    m_x = cmplx(0.0,0.0)
    m_y = cmplx(0.0,0.0)

    DO i =0,nkx0-1
        m_x(i,:,:,:,:,:) = filter_g_1(i,:,:,:,:,:)*(kxgrid(i)/kxmax0)**hypx_order
        !m_x(i,:,:,:,:,:) = g_in(i,:,:,:,:,:)*(kxgrid(i)/kxmax0)**hypx_order
    END DO

    DO j=0,nky0-1
        m_y(:,j,:,:,:,:) = filter_g_1(:,j,:,:,:,:)*(kygrid(j)/kymax0)**hypy_order
        !m_y(:,j,:,:,:,:) = g_in(:,j,:,:,:,:)*(kygrid(j)/kymax0)**hypy_order
    ENDDO

    
    ! terms needed for coefficients of the dynamic procedure
    tmp_submod = 0.0
    tmp_3D_dp = 0.0
    energy3d =  0.0
    mxmx = 0.0
    tmp_submod = real(m_x*conjg(m_x))
  
    DO l=lv1,lv2
       CALL sum3d_real(tmp_submod(:,:,:,l,0,0),mxmx(l))
    ENDDO

    
    tmp_submod = 0.0
    tmp_3D_dp = 0.0
    energy3d = 0.0
    mxmy = 0.0
    tmp_submod = real(m_x*conjg(m_y))
 
    DO l=lv1,lv2
       CALL sum3d_real(tmp_submod(:,:,:,l,0,0),mxmy(l))
    ENDDO

    tmp_submod = 0.0
    tmp_3D_dp = 0.0
    mymy = 0.0
    tmp_submod = real(m_y*conjg(m_y))

    DO l=lv1,lv2
       CALL sum3d_real(tmp_submod(:,:,:,l,0,0),mymy(l))
    ENDDO

    tmp_submod = 0.0
    tmp_3D_dp = 0.0
    energy3d  = 0.0
    Tmx  = 0.0
    tmp_submod = real(t_fgs*conjg(m_x))
 
    DO l=lv1,lv2
       CALL sum3d_real(tmp_submod(:,:,:,l,0,0),Tmx(l))
    ENDDO

    tmp_submod = 0.0
    tmp_3D_dp = 0.0
    energy3d = 0.0
    Tmy = 0.0
    tmp_submod = real(t_fgs*conjg(m_y))
 
    DO l=lv1,lv2
       CALL sum3d_real(tmp_submod(:,:,:,l,0,0),Tmy(l))
    ENDDO

    !eps_alpha = 1.0 - (sqrt(fracx*fracy)**(-(hypx_order + 1./3.)))
    eps_alpha_x = 1.0 - fracx**(-(hypx_order + 1./3.))
    eps_alpha_y = 1.0 - fracy**(-(hypy_order + 1./3.))
    !IF (mype==0) WRITE(*,*) "arrays", mymy,Tmx,mxmy,Tmy
    ! Calculation of the amplitude 
    c_x = 0.0
    c_y = 0.0

    DO l=lv1,lv2
            IF (mype==0) WRITE(*,*) "dyn",mxmy(l),mxmx(l),mymy(l), mxmy(l)*mxmy(l) - mymy(l)*mxmx(l)
            IF (abs(mxmy(l)*mxmy(l) - mymy(l)*mxmx(l)).gt.0.00001) THEN
                c_x(l) = (mymy(l)*Tmx(l) - mxmy(l)*Tmy(l)) /(mxmy(l)*mxmy(l) - mymy(l)*mxmx(l))/eps_alpha_x
                c_y(l) = (mxmx(l)*Tmy(l) - mxmy(l)*Tmx(l)) /(mxmy(l)*mxmy(l) - mymy(l)*mxmx(l))/eps_alpha_y
            ENDIF
    enddo

   !Update the values of the amplitude
    Do l=lv1,lv2
        if (c_x(l).gt.0.0) then
            hyp_x_herm1(l) = c_x(l)
        else
            hyp_x_herm1(l) = 0.0001
        endif
    
        if (c_y(l).gt.0.0) then
            hyp_y_herm1(l) = c_y(l)
        else
            hyp_y_herm1(l) = 0.0001
        endif
    enddo

!Some test
    tmp_3D_dp = 0.0
    energy3d =  0.0
    mmod = 0.0
    m_x_mod = CMPLX(0.0,0.0)
    m_y_mod = CMPLX(0.0,0.0)

    DO i =0,nkx0-1
       Do l=lv1,lv2
        m_x_mod(i,:,:,l,:,:) = -hyp_x_herm1(l)*g_in(i,:,:,l,:,:)*(kxgrid(i)/kxmax0)**hypx_order
       END DO
    END DO

    DO j=0,nky0-1
       DO l=lv1,lv2
        m_y_mod(:,j,:,l,:,:) = -hyp_y_herm1(l)*g_in(:,j,:,l,:,:)*(kygrid(j)/kymax0)**hypy_order
       END DO
    ENDDO

    tmp_3D_dp(:,:,:)=REAL(sqrt(pi)*sum(conjg(g_in(:,:,:,:,0,0))*(m_x_mod(:,:,:,:,0,0) + m_y_mod(:,:,:,:,0,0)),4))
    IF(mype==0) THEN
      DO k=0,nkz0-1
        tmp_3D_dp(:,:,k)=tmp_3D_dp(:,:,k)+REAL(pi**(0.25)*J0a(:,:)*conjg(phi(:,:,k))*(m_x_mod(:,:,k,0,0,0) + m_y_mod(:,:,k,0,0,0)))
      END DO
    END IF
 

   CALL MPI_ALLREDUCE(tmp_3D_dp,energy3d,nkx0*nky0*nkz0 &
      ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
    CALL sum3d_real(energy3d,mmod)

!Finish test

    CALL mpi_barrier(mpi_comm_world,ierr)

    !Transfer all coefficientes to mype_herm ==0
    CALL MPI_GATHER(hyp_x_herm1(lv1), lv0, MPI_DOUBLE_PRECISION, hyp_x_herm_glob(0), lv0, MPI_DOUBLE_PRECISION,0,MPI_COMM_HERM,ierr)
    CALL MPI_GATHER(hyp_y_herm1(lv1), lv0, MPI_DOUBLE_PRECISION, hyp_y_herm_glob(0), lv0, MPI_DOUBLE_PRECISION,0,MPI_COMM_HERM,ierr)
  

    !Write the amplitude to the file
    IF (mype==0) then
     !   DO l=0,nv0-1
       !     WRITE(GyroLES_FILE,"(ES12.4,I3,2E12.4)") time, l, hyp_x_herm_glob(l), hyp_y_herm_glob(l), tgs
            WRITE(GyroLES_FILE,"(5E12.4)") time, hyp_x_herm_glob(0), hyp_y_herm_glob(0), tgs, mmod
       !ENDDO
    ENDIF       
 
    Deallocate(filter_g_1)
    Deallocate(filter_phi)
    Deallocate(nlt)
    Deallocate(filter_nlt)
    Deallocate(nlt_filter_g_1)
    Deallocate(filter_nlt_filter_g_1)
    Deallocate(t_fgs)
    Deallocate(m_x)
    Deallocate(m_x_mod)
    Deallocate(m_y)
    Deallocate(m_y_mod)
    Deallocate(tmp_submod)
    Deallocate(tmp_3D_dp)
    Deallocate(energy3d)
    
End subroutine exec_GyroLES_herm

Subroutine exec_GyroLES

    real :: mymy, mxmx, mxmy, Tmx, Tmy, c_x, c_y
    real,dimension(0:nv0-1) :: hyp_x_herm_glob, hyp_y_herm_glob
    real :: eps_alpha_x, eps_alpha_y 
    Integer :: i, j, l, ierr

    Allocate(filter_g_1(0:nkx0-1,0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))
    Allocate(filter_phi(0:nkx0-1,0:nky0-1, lkz1:lkz2))
    Allocate(nlt(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2)) 
    Allocate(filter_nlt(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))     
    Allocate(nlt_filter_g_1(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))  
    Allocate(filter_nlt_filter_g_1(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))   
    Allocate(t_fgs(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))
    Allocate(m_x(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2)) 
    Allocate(m_y(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2)) 
    allocate(tmp_submod(0:nkx0-1, 0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
    allocate(tmp_3D_dp(0:nkx0-1, 0:nky0-1,lkz1:lkz2))
    allocate(energy3d(0:nkx0-1, 0:nky0-1,lkz1:lkz2))

    !Apply the filter to g_1, compute the filtered fields
    Call get_phi(g_1)
    Call filter(g_1, filter_g_1)

    ! compute the T sub grid term
    nlt = cmplx(0.0,0.0)
    Call get_rhs_nl(g_1, phi, nlt)
    Call filter(nlt, filter_nlt)
    
    nlt_filter_g_1 = cmplx(0.0,0.0)
    Call filter(g_1, filter_g_1)
    Call filter_in_phi(phi, filter_phi)
    Call get_rhs_nl(filter_g_1, filter_phi, nlt_filter_g_1)
    Call filter(nlt_filter_g_1, filter_nlt_filter_g_1)     
    t_fgs = filter_nlt - filter_nlt_filter_g_1

    ! compute hyper-diff as models, with amplitude set to 1
    
    m_x = cmplx(0.0,0.0)
    m_y = cmplx(0.0,0.0)

    DO i =0,nkx0-1
        m_x(i,:,:,:,:,:) = filter_g_1(i,:,:,:,:,:)*(kxgrid(i)/kxmax0)**hypx_order
    END DO

    DO j=0,nky0-1
        m_y(:,j,:,:,:,:) = filter_g_1(:,j,:,:,:,:)*(kygrid(j)/kymax0)**hypy_order
    ENDDO

    ! terms needed for coefficients of the dynamic procedure
    tmp_submod = 0.0
    tmp_3D_dp = 0.0
    energy3d =  0.0
    mxmx = 0.0
    tmp_submod = REAL(m_x*conjg(m_x))
  
    tmp_3D_dp(:,:,:) = sum(tmp_submod(:,:,:,:,0,0),4) 
    CALL MPI_ALLREDUCE(tmp_3d_dp,energy3d,nkx0*nky0*nkz0 &
        ,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_HERM, ierr)
 
    CALL sum3d_real(energy3d,mxmx)
    
    tmp_submod = 0.0
    tmp_3D_dp = 0.0
    energy3d = 0.0
    mxmy = 0.0
    tmp_submod = REAL(m_x*conjg(m_y))
 
    tmp_3d_dp(:,:,:) =sum(tmp_submod(:,:,:,:,0,0),4) 
    CALL MPI_ALLREDUCE(tmp_3d_dp,energy3d,nkx0*nky0*nkz0 &
        ,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_HERM, ierr)
 
    CALL sum3d_real(energy3d,mxmy)
    
    tmp_submod = 0.0
    tmp_3D_dp = 0.0
    mymy = 0.0
    tmp_submod = REAL(m_y*conjg(m_y))

    tmp_3D_dp(:,:,:) = sum(tmp_submod(:,:,:,:,0,0),4) 
    CALL MPI_ALLREDUCE(tmp_3D_dp,energy3d,nkx0*nky0*nkz0 &
        ,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_HERM, ierr)
 
    CALL sum3d_real(energy3d,mymy)
 
    tmp_submod = 0.0
    tmp_3D_dp = 0.0
    energy3d  = 0.0
    Tmx  = 0.0
    tmp_submod = REAL(t_fgs*conjg(m_x))
 
    tmp_3D_dp(:,:,:) = sum(tmp_submod(:,:,:,:,0,0),4) 
    CALL MPI_ALLREDUCE(tmp_3D_dp,energy3d,nkx0*nky0*nkz0 &
        ,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_HERM, ierr)
 
    CALL sum3d_real(energy3d,Tmx)
 
    tmp_submod = 0.0
    tmp_3D_dp = 0.0
    energy3d = 0.0
    Tmy = 0.0
    tmp_submod = REAL(t_fgs*conjg(m_y))
 
    tmp_3D_dp(:,:,:) = sum(tmp_submod(:,:,:,:,0,0),4) 
    CALL MPI_ALLREDUCE(tmp_3D_dp,energy3d,nkx0*nky0*nkz0 &
        ,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_HERM, ierr)
 
    CALL sum3d_real(energy3d,Tmy)
 
    eps_alpha_x = 1.0 - fracx**(-(hypx_order + 1./3.))
    eps_alpha_y = 1.0 - fracy**(-(hypy_order + 1./3.))
    
    c_x = 0.0
    c_y = 0.0

    c_x = (mymy*Tmx - mxmy*Tmy) /(mxmy*mxmy - mymy*mxmx)/eps_alpha_x
    c_y = (mxmx*Tmy - mxmy*Tmx) /(mxmy*mxmy - mymy*mxmx)/eps_alpha_y

    !Update the values of the amplitude
    if (c_x.gt.0.0) then
        hyp_x_herm = c_x
    else
        hyp_x_herm = 0.0001
    endif
    
    if (c_y.gt.0.0) then
       hyp_y_herm = c_y
    else
       hyp_y_herm = 0.0001
    endif

    CALL mpi_barrier(mpi_comm_world,ierr)

    IF (mype==0) then
       WRITE(GyroLES_FILE,"(3E12.4)") time, hyp_x_herm, hyp_y_herm
    ENDIF       
 
    Deallocate(filter_g_1)
    Deallocate(filter_phi)
    Deallocate(nlt)
    Deallocate(filter_nlt)
    Deallocate(nlt_filter_g_1)
    Deallocate(filter_nlt_filter_g_1)
    Deallocate(t_fgs)
    Deallocate(m_x)
    Deallocate(m_y)
    Deallocate(tmp_submod)
    Deallocate(tmp_3D_dp)
    Deallocate(energy3d)
    
End subroutine exec_GyroLES

Subroutine exec_GyroLES_z

    real :: mxmx, Tmx, c_x
    real :: eps_alpha_x 
    Integer :: i, j, k, l, ierr

    Allocate(filter_g_1(0:nkx0-1,0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))
    Allocate(filter_phi(0:nkx0-1,0:nky0-1, lkz1:lkz2))
    Allocate(nlt(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2)) 
    Allocate(filter_nlt(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))     
    Allocate(nlt_filter_g_1(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))  
    Allocate(filter_nlt_filter_g_1(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))   
    Allocate(t_fgs(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))
    Allocate(m_x(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2)) 
    allocate(tmp_submod(0:nkx0-1, 0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
    allocate(tmp_3D_dp(0:nkx0-1, 0:nky0-1,lkz1:lkz2))
    allocate(energy3d(0:nkx0-1, 0:nky0-1,lkz1:lkz2))

    !Apply the filter to g_1, compute the filtered fields
    Call get_phi(g_1)
    Call filter(g_1, filter_g_1)

    ! compute the T sub grid term
    nlt = cmplx(0.0,0.0)
    Call get_rhs_nl(g_1, phi, nlt)
    Call filter(nlt, filter_nlt)
    
    nlt_filter_g_1 = cmplx(0.0,0.0)
    Call filter(g_1, filter_g_1)
    Call filter_in_phi(phi, filter_phi)
    Call get_rhs_nl(filter_g_1, filter_phi, nlt_filter_g_1)
    Call filter(nlt_filter_g_1, filter_nlt_filter_g_1)     
    t_fgs = filter_nlt - filter_nlt_filter_g_1

    ! compute hyper-diff as models, with amplitude set to 1
    
    m_x = cmplx(0.0,0.0)

    DO k =0,nkz0-1
        m_x(:,:,k,:,:,:) = filter_g_1(:,:,k,:,:,:)*(kzgrid(k)/kzmax0)**hypz_order
    END DO

    ! terms needed for coefficients of the dynamic procedure
    tmp_submod = 0.0
    tmp_3D_dp = 0.0
    energy3d =  0.0
    mxmx = 0.0
    tmp_submod = REAL(m_x*conjg(m_x))
  
    tmp_3D_dp(:,:,:) = sum(tmp_submod(:,:,:,:,0,0),4) 
    CALL MPI_ALLREDUCE(tmp_3d_dp,energy3d,nkx0*nky0*nkz0 &
        ,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_HERM, ierr)
 
    CALL sum3d_real(energy3d,mxmx)
    
    tmp_submod = 0.0
    tmp_3D_dp = 0.0
    energy3d  = 0.0
    Tmx  = 0.0
    tmp_submod = REAL(t_fgs*conjg(m_x))
 
    tmp_3D_dp(:,:,:) = sum(tmp_submod(:,:,:,:,0,0),4) 
    CALL MPI_ALLREDUCE(tmp_3D_dp,energy3d,nkx0*nky0*nkz0 &
        ,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_HERM, ierr)
 
    CALL sum3d_real(energy3d,Tmx)
 
    eps_alpha_x = 1.0 - fracx**(-(hypz_order + 1./3.))
    
    c_x = 0.0

    c_x = Tmx /mxmx/eps_alpha_x

    !Update the values of the amplitude
    if (c_x.gt.0.0) then
        hyp_z = c_x
    else
        hyp_z = 0.0001
    endif
    
!    CALL mpi_barrier(mpi_comm_world,ierr)

    IF (mype==0) then
       WRITE(GyroLES_FILE,"(2E12.4)") time, hyp_z
    ENDIF       
 
    Deallocate(filter_g_1)
    Deallocate(filter_phi)
    Deallocate(nlt)
    Deallocate(filter_nlt)
    Deallocate(nlt_filter_g_1)
    Deallocate(filter_nlt_filter_g_1)
    Deallocate(t_fgs)
    Deallocate(m_x)
    Deallocate(tmp_submod)
    Deallocate(tmp_3D_dp)
    Deallocate(energy3d)
    
End subroutine exec_GyroLES_z

Subroutine Correlation 

    real :: mxmx, Tmx, TT, mymy, Tmy, mzmz, Tmz, c_x, c_y, c_z, subgrid
    Integer :: i, j, k, l, ierr

    Allocate(filter_g_1(0:nkx0-1,0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))
    Allocate(filter_phi(0:nkx0-1,0:nky0-1, lkz1:lkz2))
    Allocate(nlt(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2)) 
    Allocate(filter_nlt(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))     
    Allocate(nlt_filter_g_1(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))  
    Allocate(filter_nlt_filter_g_1(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))   
    Allocate(t_fgs(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2))
    Allocate(m_x(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2)) 
    Allocate(m_y(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2)) 
    Allocate(m_z(0:nkx0-1, 0:nky0-1, lkz1:lkz2, lv1:lv2, lh1:lh2, ls1:ls2)) 
    allocate(tmp_submod(0:nkx0-1, 0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
    allocate(tmp_3D_dp(0:nkx0-1, 0:nky0-1,lkz1:lkz2))
    allocate(energy3d(0:nkx0-1, 0:nky0-1,lkz1:lkz2))

    !Apply the filter to g_1, compute the filtered fields
    Call get_phi(g_1)
    Call filter(g_1, filter_g_1)

    ! compute the T sub grid term
    nlt = cmplx(0.0,0.0)
    Call get_rhs_nl(g_1, phi, nlt)
    Call filter(nlt, filter_nlt)
    
    nlt_filter_g_1 = cmplx(0.0,0.0)
    Call filter(g_1, filter_g_1)
    Call filter_in_phi(phi, filter_phi)
    Call get_rhs_nl(filter_g_1, filter_phi, nlt_filter_g_1)
    Call filter(nlt_filter_g_1, filter_nlt_filter_g_1)     
    t_fgs = filter_nlt - filter_nlt_filter_g_1

    ! compute hyper-diff as models, with amplitude set to 1
    
    m_x = cmplx(0.0,0.0)
    m_y = cmplx(0.0,0.0)
    m_z = cmplx(0.0,0.0)

    DO i =0,nkx0-1
        m_x(i,:,:,:,:,:) = filter_g_1(i,:,:,:,:,:)*(kxgrid(i)/kxmax0)**hypx_order
    END DO

    DO j =0,nky0-1
        m_y(:,j,:,:,:,:) = filter_g_1(:,j,:,:,:,:)*(kygrid(j)/kymax0)**hypy_order
    END DO

    DO k =0,nkz0-1
        m_z(:,:,k,:,:,:) = filter_g_1(:,:,k,:,:,:)*(kzgrid(k)/kzmax0)**hypz_order
    END DO


    ! terms needed for the correlations
    ! Mx correlation
    tmp_submod = 0.0
    tmp_3D_dp = 0.0
    energy3d =  0.0
    mxmx = 0.0
    tmp_submod = REAL(m_x*conjg(m_x))
  
    tmp_3D_dp(:,:,:) = sum(tmp_submod(:,:,:,:,0,0),4) 
    CALL MPI_ALLREDUCE(tmp_3d_dp,energy3d,nkx0*nky0*nkz0 &
        ,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_HERM, ierr)
 
    CALL sum3d_real(energy3d,mxmx)
    
    tmp_submod = 0.0
    tmp_3D_dp = 0.0
    energy3d  = 0.0
    Tmx  = 0.0
    tmp_submod = REAL(t_fgs*conjg(m_x))
 
    tmp_3D_dp(:,:,:) = sum(tmp_submod(:,:,:,:,0,0),4) 
    CALL MPI_ALLREDUCE(tmp_3D_dp,energy3d,nkx0*nky0*nkz0 &
        ,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_HERM, ierr)
 
    CALL sum3d_real(energy3d,Tmx)
 
    tmp_submod = 0.0
    tmp_3D_dp = 0.0
    energy3d  = 0.0
    TT  = 0.0
    tmp_submod = REAL(t_fgs*conjg(t_fgs))
 
    tmp_3D_dp(:,:,:) = sum(tmp_submod(:,:,:,:,0,0),4) 
    CALL MPI_ALLREDUCE(tmp_3D_dp,energy3d,nkx0*nky0*nkz0 &
        ,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_HERM, ierr)
 
    CALL sum3d_real(energy3d,TT)   
    
    !My correlation
    tmp_submod = 0.0
    tmp_3D_dp = 0.0
    energy3d =  0.0
    mymy = 0.0
    tmp_submod = REAL(m_y*conjg(m_y))
  
    tmp_3D_dp(:,:,:) = sum(tmp_submod(:,:,:,:,0,0),4) 
    CALL MPI_ALLREDUCE(tmp_3d_dp,energy3d,nkx0*nky0*nkz0 &
        ,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_HERM, ierr)
 
    CALL sum3d_real(energy3d,mymy)
    
    tmp_submod = 0.0
    tmp_3D_dp = 0.0
    energy3d  = 0.0
    Tmy  = 0.0
    tmp_submod = REAL(t_fgs*conjg(m_y))
 
    tmp_3D_dp(:,:,:) = sum(tmp_submod(:,:,:,:,0,0),4) 
    CALL MPI_ALLREDUCE(tmp_3D_dp,energy3d,nkx0*nky0*nkz0 &
        ,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_HERM, ierr)
 
    CALL sum3d_real(energy3d,Tmy)
 
    !Mz correlation
    tmp_submod = 0.0
    tmp_3D_dp = 0.0
    energy3d =  0.0
    mzmz = 0.0
    tmp_submod = REAL(m_z*conjg(m_z))
  
    tmp_3D_dp(:,:,:) = sum(tmp_submod(:,:,:,:,0,0),4) 
    CALL MPI_ALLREDUCE(tmp_3d_dp,energy3d,nkx0*nky0*nkz0 &
        ,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_HERM, ierr)
 
    CALL sum3d_real(energy3d,mzmz)
    
    tmp_submod = 0.0
    tmp_3D_dp = 0.0
    energy3d  = 0.0
    Tmz  = 0.0
    tmp_submod = REAL(t_fgs*conjg(m_z))
 
    tmp_3D_dp(:,:,:) = sum(tmp_submod(:,:,:,:,0,0),4) 
    CALL MPI_ALLREDUCE(tmp_3D_dp,energy3d,nkx0*nky0*nkz0 &
        ,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_HERM, ierr)
 
    CALL sum3d_real(energy3d,Tmz)
 
    !Subgrid contribution
    tmp_3D_dp = 0.0
    energy3d  = 0.0
    tmp_3d_dp(:,:,:)=real(sqrt(pi)*sum(conjg(g_1(:,:,:,:,0,0))*t_fgs(:,:,:,:,0,0),4))
    if(mype==0) then
      do k=0,nkz0-1
        tmp_3d_dp(:,:,k)=tmp_3d_dp(:,:,k)+real(pi**(0.25)*j0a(:,:)*conjg(phi(:,:,k))*t_fgs(:,:,k,0,0,0))
      end do
    end if
    call mpi_allreduce(tmp_3d_dp,energy3d,nkx0*nky0*nkz0 &
      ,mpi_double_precision,mpi_sum,mpi_comm_herm,ierr)
    call sum3d_real(energy3d,subgrid)

    
    c_x = Tmx /(sqrt(mxmx)*sqrt(TT))
    c_y = Tmy /(sqrt(mymy)*sqrt(TT))
    c_z = Tmz /(sqrt(mymy)*sqrt(TT))

    IF (mype==0) then
       WRITE(CORRELATION_FILE,"(5E12.4)") time, subgrid, c_x, c_y, c_z
    ENDIF       
 
    Deallocate(filter_g_1)
    Deallocate(filter_phi)
    Deallocate(nlt)
    Deallocate(filter_nlt)
    Deallocate(nlt_filter_g_1)
    Deallocate(filter_nlt_filter_g_1)
    Deallocate(t_fgs)
    Deallocate(m_x)
    Deallocate(m_y)
    Deallocate(m_z)
    Deallocate(tmp_submod)
    Deallocate(tmp_3D_dp)
    Deallocate(energy3d)
    
End subroutine Correlation

subroutine finalize_GyroLES

        IF (mype==0) then    
            Close (GyroLES_FILE)
        endif
        
end subroutine finalize_GyroLES

subroutine finalize_corr

        IF (mype==0) then    
            Close (CORRELATION_FILE)
        endif
        
end subroutine finalize_corr


!!!**************************************************************************!!!
!!!********************* Subroutine that the module needs ************************!!! 

!!!**************************************************************************!!!
!!!*********************** Calculate the filter  ************************
subroutine filter(f_in,f_out)

    implicit none
    
    complex, dimension(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2), intent(in) :: f_in
    complex, dimension(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2), intent(out) :: f_out
    Integer :: fkx, fky
    
    ! Filter in the kx direction, it gives you the indice of the maximum value of the ky mode you keep
        fky = int((fracy*nky0/2.))
        ! Filter in the kx direction, it gives you the indice of the maximum value of the kx mode you keep 
        fkx = int(fracx*(nkx0-1))

        !Copy the input function to the output function FGS
        f_out = f_in
        
        !Filter FGS in the kx direction for modes higher than fkx and lower than -fkx 

        f_out(fkx + 1 : nkx0 -1,:,:,:,:,:) = CMPLX(0.,0.)
        f_out( :,fky + 1 : nky0 -1 - fky,:,:,:,:) = CMPLX(0.,0.)
   
End Subroutine filter

subroutine filter_in_phi(f_in,f_out)

    implicit none
    
    complex, dimension(0:nkx0-1,0:nky0-1,lkz1:lkz2), intent(in) :: f_in
    complex, dimension(0:nkx0-1,0:nky0-1,lkz1:lkz2), intent(out) :: f_out
    real, dimension(0:nkx0-1,0:nky0) :: int_test
    real, dimension(0:nkx0-1,0:nky0) :: out_test
    Integer :: fkx, fky, i,j
   
    int_test = 1.0 
    out_test = 1.0 
    ! Filter in the kx direction, it gives you the indice of the maximum value of the ky mode you keep
        fky = int((fracy*nky0/2.))
        ! Filter in the kx direction, it gives you the indice of the maximum value of the kx mode you keep 
        fkx = int(fracx*(nkx0-1))

        !Copy the input function to the output function FGS
        f_out = f_in
        
        !Filter FGS in the kx direction for modes higher than fkx and lower than -fkx 
 
        f_out(fkx + 1 : nkx0 -1,:,:) = CMPLX(0.,0.)
        f_out( :,fky + 1 : nky0 -1 - fky,:) = CMPLX(0.,0.)
        
        !Some test to check if it doing the correct filtering
        !out_test(fkx + 1 : nkx0 -1,:) = 0.
        !out_test( :,fky + 1 : nky0 -1 - fky) = 0.
   
        !DO j = 0,nky0-1
        !   IF (mype==0) WRITE(*,*) j, fky, out_test(0,j), int_test(0,j)
        !ENDDO

End Subroutine filter_in_phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                               sum3d_real                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE sum3d_real(vec3d,scal)
  IMPLICIT NONE

  REAL, INTENT(in) :: vec3d(0:nkx0-1,0:nky0-1,0:nkz0-1)
  REAL, INTENT(out) :: scal
  INTEGER :: i

  scal=sum(sum(vec3d(0,:,:),1),1)
  DO i=1,nkx0-1
    scal=scal+2.0*sum(sum(vec3d(i,:,:),1),1)
  END DO

END SUBROUTINE sum3d_real




End Module Gyro_LES
