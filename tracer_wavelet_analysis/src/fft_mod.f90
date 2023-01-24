!#################################################
! 3D Fast Fourier Transform                      #
! Reference: Numerical Recipes, Press et' al'    #
!#################################################

MODULE fft_mod

USE constants_mod

IMPLICIT NONE


INTEGER :: wx, wy, wz, wmax



CONTAINS


SUBROUTINE fft_init(nwx_in, nwy_in, nwz_in)

    INTEGER, INTENT(IN) :: nwx_in, nwy_in, nwz_in

    ! set fwt grid size
    wx = nwx_in
    wy = nwy_in
    wz = nwz_in
    
    wmax = MAX(wx,wy,wz)
    
END SUBROUTINE fft_init


! this subroutines performs FFT on a 3d array of scalars
SUBROUTINE fft_3d(fx, fk_re, fk_im)

    REAL, INTENT(IN) :: fx(1:wx,1:wy,1:wz) 
    REAL, INTENT(INOUT) :: fk_re(1:wx,1:wy,1:wz), fk_im(1:wx,1:wy,1:wz) 
    REAL :: dft_buffer_in_re(1:wmax), dft_buffer_in_im(1:wmax), dft_buffer_out_re(1:wmax),dft_buffer_out_im(1:wmax)   
    INTEGER :: ix, iy, iz
  
    
    ! clear output array    
    fk_re = 0.d0
    fk_im = 0.d0
        
    !PRINT*,' Doing x pass..'    
        
    ! FFT in x-direction
    DO iz = 1, wz
        DO iy = 1, wy
                            
            ! copy strips-along x into 1d buffer    
            DO ix = 1, wx
                dft_buffer_in_re(ix) = fx(ix,iy,iz)
                dft_buffer_in_im(ix) = 0.d0
            END DO
    
            ! perform 1D inverse DFT 
            !CALL direct_dft_1d(dft_buffer_in_re,dft_buffer_in_im,dft_buffer_out_re,dft_buffer_out_im)
            CALL cfft_1d(dft_buffer_in_re,dft_buffer_in_im,dft_buffer_out_re,dft_buffer_out_im,wx)
    
            ! copy back into delv array 
            DO ix = 1, wx
                fk_re(ix,iy,iz) = dft_buffer_out_re(ix)  
                fk_im(ix,iy,iz) = dft_buffer_out_im(ix)  
            END DO

        END DO       
    END DO

    !PRINT*,' Doing y pass..'    

    ! DFT in y-direction
    DO iz = 1, wz
        DO ix = 1, wx
                            
            ! copy strips-along y into 1d buffer    
            DO iy = 1, wy
                dft_buffer_in_re(iy) = fk_re(ix,iy,iz)
                dft_buffer_in_im(iy) = fk_im(ix,iy,iz)
            END DO
    
            ! perform 1D inverse DFT 
            !CALL direct_dft_1d(dft_buffer_in_re,dft_buffer_in_im,dft_buffer_out_re,dft_buffer_out_im)
            CALL cfft_1d(dft_buffer_in_re,dft_buffer_in_im,dft_buffer_out_re,dft_buffer_out_im,wy)
                
            ! copy back into delvk array 
            DO iy = 1, wy
                fk_re(ix,iy,iz) = dft_buffer_out_re(iy)  
                fk_im(ix,iy,iz) = dft_buffer_out_im(iy)  
            END DO            
            
        END DO
    END DO

    !PRINT*,' Doing z pass..'    

    ! DFT in z-direction. (also apply the IDFT prefactor of 1/nx_tot*nx_tot*nx_tot)
    DO iy = 1, wy
        DO ix = 1, wx
                            
            ! copy strips-along z into 1d buffer    
            DO iz = 1,wz
                dft_buffer_in_re(iz) = fk_re(ix,iy,iz)
                dft_buffer_in_im(iz) = fk_im(ix,iy,iz)
            END DO
    
            ! perform 1D inverse DFT 
            !CALL direct_dft_1d(dft_buffer_in_re,dft_buffer_in_im,dft_buffer_out_re,dft_buffer_out_im)
            CALL cfft_1d(dft_buffer_in_re,dft_buffer_in_im,dft_buffer_out_re,dft_buffer_out_im,wz)
                
            DO iz = 1, wz
                fk_re(ix,iy,iz) = dft_buffer_out_re(iz)  
                fk_im(ix,iy,iz) = dft_buffer_out_im(iz)               
            END DO
                
        END DO        
    END DO


END SUBROUTINE fft_3d


! This subroutine performs a 1D discrete fourier transform via direct summation
SUBROUTINE direct_dft_1d(in_re, in_im, out_re, out_im)


    REAL(4), INTENT(IN) :: in_re(0:wx-1), in_im(0:wx-1)  
    REAL(4), INTENT(INOUT)  :: out_re(0:wx-1), out_im(0:wx-1)    
    INTEGER :: sgn = 1   ! DFT for sgn = 1 and Inverse DFT for sgn = -1
    INTEGER :: ix, kx
    REAL(4) :: theta, theta0
    
    theta0 =  TWOPI / DBLE(wx)
    
    ! clear output arrays
    out_re = 0.d0
    out_im = 0.d0
        
    DO kx = 0, wx-1
        DO ix = 0, wx-1 
            theta = theta0 * DBLE(ix * kx)
            out_re(kx) = out_re(kx) + (in_re(ix) * COS(theta) + DBLE(sgn) * in_im(ix) * SIN(theta))
            out_im(kx) = out_im(kx) + (in_im(ix) * COS(theta) - DBLE(sgn) * in_re(ix) * SIN(theta))
        END DO
    END DO   
  

END SUBROUTINE direct_dft_1d



! 1D Daniel-Lanczos FFT algorithm for complex input
SUBROUTINE cfft_1d(in_re, in_im, out_re, out_im, nx)

    REAL(4), INTENT(IN) :: in_re(:), in_im(:)  
    REAL(4), INTENT(INOUT) :: out_re(:), out_im(:)  
    INTEGER, INTENT(IN) :: nx
    INTEGER :: sgn = 1   ! DFT for sgn = 1 and Inverse DFT for sgn = -1
    INTEGER :: ix, kx
 

    REAL(4) :: buffer(1:2*nx)   ! input array gets replaced by output
    REAL(4) :: tempr, tempi, theta, wi, wr, wpi, wpr, wtemp
    INTEGER :: i, j, n, m, mmax, istep
    INTEGER :: i1, i2, i3, i4, n2p3
    REAL*4 :: c1, c2, h1i, h1r, h2i, h2r, wis, wrs
    
 
    ! clear output arrays
    out_re = 0.d0
    out_im = 0.d0
    
    
    ! load input array into work buffer
    ix = 1
    DO i = 1, nx
        buffer(ix)   = in_re(i)   
        buffer(ix+1) = in_im(i) 
        ix = ix+2        
    END DO
  
  
    !***************************************
    ! Sort input array in bit-reversed order  
    !***************************************    
    n = 2*nx
    j = 1
    
    DO i = 1, n, 2
    
        IF(j .GT. i) THEN  ! swap the two complex numbers
            tempr = buffer(j)
            tempi = buffer(j+1)
            buffer(j) = buffer(i)
            buffer(j+1) = buffer(i+1)
            buffer(i) = tempr
            buffer(i+1) = tempi                    
        END IF
    
        m = n/2
        DO WHILE ((m .GT. 2) .AND. (j .GT. m))
            j = j - m 
            m = m / 2
        END DO
        j = j + m
    END DO
      
    !********************************************************************************
    ! Using Danielson-Laczos lemma, compute the DFT by summing up the 1-pt base DFT's
    !********************************************************************************     
    mmax = 2
    
    DO WHILE(n .GT. mmax) 
    
        ! initialize for trigonometric recurrence
        istep = 2 * mmax
        theta = TWOPI / (sgn*mmax)  
        wpr = -2.d0 * SIN(0.5d0 * theta)**2 
        wpi =  SIN(theta)
        wr = 1.d0
        wi = 0.d0
        
        DO m = 1, mmax, 2 
            DO i = m, n, istep
       
                j = i + mmax
                
                ! apply Danielson-Lanczos lemma
                tempr = wr*buffer(j) - wi*buffer(j+1)
                tempi = wr*buffer(j+1) + wi*buffer(j)
                buffer(j) = buffer(i) - tempr 
                buffer(j+1) = buffer(i+1) - tempi 
                buffer(i) = buffer(i) + tempr 
                buffer(i+1) = buffer(i+1) + tempi 
                
            END DO
            
            wtemp = wr
            wr = wr*wpr - wi*wpi + wr
            wi = wi*wpr + wtemp*wpi + wi
            
        END DO
        mmax = istep
        
    END DO
      
    ix = 1
    DO kx = 1, nx
        out_re(kx) =  buffer(ix)
        out_im(kx) =  buffer(ix+1) 
        ix = ix + 2
    END DO


END SUBROUTINE cfft_1d




END MODULE fft_mod