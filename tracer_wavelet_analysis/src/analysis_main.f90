PROGRAM analysis

USE constants_mod
USE fft_mod
USE fwt_mod
USE readfile_mod
USE grid_arrays_mod
USE OMP_LIB
USE tracer_reconstruction_mod

IMPLICIT NONE

!##########################################################################################

                       
REAL*4 :: Lx, Ly, Lz, dx
REAL*4 :: x, y, z, rho, bmag, bhat(3), vdotb
REAL*4 :: total_ke, v_rms, b_rms
INTEGER :: i, j, k, nd
REAL(8) :: t0, t1, t2, t2a, t3, t4, t5, t6, t7, t8, t9, t10

!##########################################################################################
!##########################################################################################

t0 = OMP_GET_WTIME()


! initialize grid arrays and acquire tracer data from file
CALL create_patch_arrays()


! grid size
Lx = 1.d0
Ly = 1.d0
Lz = 1.d0
dx = Lx/DBLE(nranks_x*nx)


t1 = OMP_GET_WTIME()

!CALL readfile_wombat(dump_num, fx)

! perform tracer wavelet reconstruction of the velocity field

t2 =  OMP_GET_WTIME()

!CALL tracer_reconstruct()

t3 = OMP_GET_WTIME()

!PRINT*,'Performing MHD mode decomposition..'
!CALL mode_decomposition_wavelet(nd)

! FFT of velocity field

!PRINT*,'Computing FFT of v...'
!CALL fft_init(nxtot, nytot, nztot)

!buffer(:,:,:) = fx(:,:,:,2)
!CALL fft_3d(buffer, vxk_re, vxk_im)

!buffer(:,:,:) = fx(:,:,:,3)
!CALL fft_3d(buffer,vyk_re, vyk_im)

!buffer(:,:,:) = fx(:,:,:,4)
!CALL fft_3d(buffer, vzk_re, vzk_im)

!PRINT*,'Computing Pk_v...'

!CALL compute_pk(dump_num, 1, vxk_re, vxk_im, vyk_re, vyk_im, vzk_re, vzk_im)

!PRINT*,'Computing fourier mode decomposition...'
!CALL mode_decomposition_fourier(dump_num, .FALSE.)


!GO TO 111

!PRINT*,'Computing FFT of wavelet reconstructed v...'

!CALL fft_init(nxtot, nytot, nztot)

!buffer(:,:,:) = fx(:,:,:,8)
!CALL fft_3d(buffer, vxk_re, vxk_im)

!buffer(:,:,:) = fx(:,:,:,9)
!CALL fft_3d(buffer,vyk_re, vyk_im)

!buffer(:,:,:) = fx(:,:,:,10)
!CALL fft_3d(buffer, vzk_re, vzk_im)

!PRINT*,'Computing Pk_v_wavelet...'

!CALL compute_pk(dump_num, 3, vxk_re, vxk_im, vyk_re, vyk_im, vzk_re, vzk_im)

!PRINT*,'Computing fourier mode decomposition...'
!CALL mode_decomposition_fourier(dump_num, .TRUE.)

!111 CONTINUE


!PRINT*,'Computing tracer wavelet mode decomposition.'
!CALL mode_decomposition_tracer_wavelets(dump_num)


!PRINT*,'Computing FFT of B...'

!CALL fft_init(nxtot, nytot, nztot)

!buffer(:,:,:) = fx(:,:,:,5)
!CALL fft_3d(buffer, bxk_re, bxk_im)

!buffer(:,:,:) = fx(:,:,:,6)
!CALL fft_3d(buffer,byk_re, byk_im)

!buffer(:,:,:) = fx(:,:,:,7)
!CALL fft_3d(buffer, bzk_re, bzk_im)

!PRINT*,'Computing Pk_b...'

!CALL compute_pk(dump_num, 2, bxk_re, bxk_im, byk_re, byk_im, bzk_re, bzk_im)


!PRINT*,'Computing FFT of rho...'


!buffer(:,:,:) = fx(:,:,:,1)
!bxk_re = 0.d0
!bxk_im = 0.d0
!byk_re = 0.d0
!byk_im = 0.d0
!bzk_re = 0.d0
!bzk_im = 0.d0
!CALL fft_3d(buffer, bxk_re, bxk_im)

!PRINT*,'Computing Pk_rho..'

!CALL compute_pk(dump_num, 0, bxk_re, bxk_im, byk_re, byk_im, bzk_re, bzk_im)


t4 = OMP_GET_WTIME()

t5 = OMP_GET_WTIME()

!##############################################################################################################

t6 = OMP_GET_WTIME()

t7 = OMP_GET_WTIME()

!##############################################################################################################

t8 = OMP_GET_WTIME()



t9 = OMP_GET_WTIME()
!##############################################################################################################


!##############################################################################################################


    PRINT*,''
    PRINT*,'Workpool initialization time (sec) = ', t1-t0
    PRINT*,'File read time (sec)               = ', t2-t1
    PRINT*,'Tracer reconstruction time (sec)   = ', t3-t2
    PRINT*,'Mode decomposition time (sec)      = ', t4-t3
    WRITE(*,'(" Total time elapsed = ",i3," hour ",i3," min ",i3, " sec")') INT(t9-t0)/3600 , MOD(INT(t9-t0),3600)/60, MOD(INT(t9-t0),60) 
    PRINT*,''


CALL destroy_patch_arrays()

    
PRINT*,'Done!'
PRINT*,''


!##########################################################################################


CONTAINS


! fourier MHD wave mode decomposition (Ref: Cho & Lazarian, MNRAS, 345, 2003)
SUBROUTINE mode_decomposition_fourier(t_dump, wavelet_reconstruction)

    INTEGER, INTENT(IN) :: t_dump 
    LOGICAL, INTENT(IN) :: wavelet_reconstruction
    REAL*4 :: bmag, B0(3), B0_hat(3), alpha, D, zet_a(3), zet_f(3), zet_s(3), khat(3), &
              kdotB, kvec(3), kpar(3), kperp(3), kpar_hat(3), kperp_hat(3), tmp1, tmp2, tmp3, tmp4, rhok_f(2), rhok_s(2), bk_a(2), bk_f(2), bk_s(2), &
              delvk_f, delvk_s, ca, cf, cs, icf, ics, vsqr_a, vsqr_f, vsqr_s, rhosqr_f, rhosqr_s, bsqr_a, bsqr_f, bsqr_s, sound_speed  
    INTEGER :: kx, ky, kz, kbins, ik
    REAL*4 :: kmin, kmax, dk, k, kpar_sqr, kperp_sqr, kpark_sqr, kperpk_sqr, k0, dv_shell 
    CHARACTER(LEN=300) :: filename
    CHARACTER(LEN=6) :: uniti
    REAL*4, ALLOCATABLE :: Pk_v(:,:), Pk_b(:,:), Pk_rho(:,:) 


    IF(t_dump<10) THEN
        WRITE(uniti,'(I1.1)') t_dump
    ELSE IF(t_dump>=10 .and. t_dump<100) THEN
        WRITE(uniti,'(I2.2)') t_dump
    ELSE IF(t_dump>=100 .and. t_dump<1000) THEN
        WRITE (uniti,'(I3.3)') t_dump
    ELSE IF(t_dump>=1000 .and. t_dump<10000) THEN
        WRITE (uniti,'(I4.3)') t_dump
    ELSE IF(t_dump>=10000 .and. t_dump<100000) THEN
        WRITE (uniti,'(I5.3)') t_dump  
    END IF 
    
    ! set number of bins (i.e. k-shells)
    kbins = nmax
    
    ALLOCATE(Pk_v(1:kbins,3),Pk_b(1:kbins,3),Pk_rho(1:kbins,2))
    Pk_v = 0.d0
    Pk_b = 0.d0
    Pk_rho = 0.d0
    
    ! k band parameters
    k0 = TWOPI/Lx
    kmin = 0.d0
    kmax = SQRT(3.d0)*k0*(nmax/2)    
    dk = (kmax-kmin)/(kbins) ! uinform bin width

    B0(:) = 0.d0

    ! compute global mean magnetic field
    DO k = 1, nztot
        DO j = 1, nytot
            DO i = 1, nxtot

              B0(1) = B0(1) + fx(i,j,k,5)
              B0(2) = B0(2) + fx(i,j,k,6)
              B0(3) = B0(3) + fx(i,j,k,7)
            
            END DO
        END DO
    END DO
    B0(:) = B0(:) / DBLE(nxtot*nytot*nztot) 
    
    ! mean B-field hard coded for now...
    !B0 = (/ 1.0, 0.0, 0.0 /)
    
    bmag = SQRT(SUM(B0**2))
    B0_hat(:) = B0(:) / bmag
    
    alpha = 1.d0/SUM(B0**2)
    sound_speed = 1.D0

    PRINT*,''
    PRINT*,'<B0> = ',B0
    PRINT*,''

    ! loop over k space grid
    DO kz = 0, nmax/2-1
        DO ky = 0, nmax/2-1
            DO kx = 0, nmax/2-1

                IF((kx .EQ. 0) .AND. (ky .EQ. 0)  .AND. (kz .EQ. 0)) CYCLE

                k = k0 * SQRT(SNGL(kx**2 + ky**2 + kz**2))
                
                kvec(1) = k0*kx
                kvec(2) = k0*ky
                kvec(3) = k0*kz
                
                khat(:) = kvec(:)/k 


                ! paralell and perpendicular components of wave vector     
                kpar_hat(:)  = B0(:) / bmag
                kpar(:)      = (kvec(1)*kpar_hat(1)+kvec(2)*kpar_hat(2)+kvec(3)*kpar_hat(3)) * kpar_hat(:)
                kperp(:)     = kvec(:) - kpar(:)
                kperp_hat(:) = 0.d0
                IF(SUM(kperp**2) .GT. 1.d-20) kperp_hat(:) = kperp(:) / SQRT(SUM(kperp**2))
               
                kpar_sqr   = SUM(kpar**2) 
                kpark_sqr  = SUM(kpar**2)/k**2 
                kperp_sqr  = SUM(kperp**2)
                kperpk_sqr = SUM(kperp**2)/k**2              
              
                !*******************************************************
                ! check for parallel and perp propagation limiting cases
                !*******************************************************  
                IF(kperpk_sqr .LT. 1.d-6*kpark_sqr) THEN       ! parallel propagation
                    kperp(:) = 0.d0
                    kperp_hat(:) = 0.d0
                ELSE IF(kpark_sqr .LT. 1.d-6*kperpk_sqr) THEN  ! perp propagation
                    kpar(:) = 0.d0
                    kpar_hat(:) = 0.d0
                END IF
                
                
                ! compute displacement unit vectors for all three wave modes
                D = MAX(0.d0, (1.d0+alpha)**2 - 4.d0*alpha*SUM(kpar**2)/(k**2))              
                
                zet_a(1) = kperp_hat(2)*B0_hat(3) - kperp_hat(3)*B0_hat(2)
                zet_a(2) = kperp_hat(3)*B0_hat(1) - kperp_hat(1)*B0_hat(3)
                zet_a(3) = kperp_hat(1)*B0_hat(2) - kperp_hat(2)*B0_hat(1)
                                
                zet_f(:) = (-1.d0+alpha+SQRT(D))*kpar(:) + (1.d0+alpha+SQRT(D))*kperp(:) 
                IF(SUM(zet_f**2) .GT. 1.d-20) zet_f(:) = zet_f(:)/(SQRT(SUM(zet_f**2)))
                              
                zet_s(:) = (-1.d0+alpha-SQRT(D))*kpar(:) + (1.d0+alpha-SQRT(D))*kperp(:)  
                IF(SUM(zet_s**2) .GT. 1.d-20) zet_s(:) = zet_s(:)/(SQRT(SUM(zet_s**2)))

                !********************************************************************************                
                ! project velocity fourier components onto displacement unit vectors of each mode
                !******************************************************************************** 

                ! alfven mode
                va_k_re(kx,ky,kz) = vxk_re(kx,ky,kz)*zet_a(1) + vyk_re(kx,ky,kz)*zet_a(2) + vzk_re(kx,ky,kz)*zet_a(3)
                va_k_im(kx,ky,kz) = vxk_im(kx,ky,kz)*zet_a(1) + vyk_im(kx,ky,kz)*zet_a(2) + vzk_im(kx,ky,kz)*zet_a(3)      
                
                ! fast mode
                vf_k_re(kx,ky,kz) = vxk_re(kx,ky,kz)*zet_f(1) + vyk_re(kx,ky,kz)*zet_f(2) + vzk_re(kx,ky,kz)*zet_f(3)  
                vf_k_im(kx,ky,kz) = vxk_im(kx,ky,kz)*zet_f(1) + vyk_im(kx,ky,kz)*zet_f(2) + vzk_im(kx,ky,kz)*zet_f(3) 
                
                ! slow mode
                vs_k_re(kx,ky,kz) = vxk_re(kx,ky,kz)*zet_s(1) + vyk_re(kx,ky,kz)*zet_s(2) + vzk_re(kx,ky,kz)*zet_s(3)                 
                vs_k_im(kx,ky,kz) = vxk_im(kx,ky,kz)*zet_s(1) + vyk_im(kx,ky,kz)*zet_s(2) + vzk_im(kx,ky,kz)*zet_s(3)                 
               
            END DO
        END DO
    END DO

    ! loop over k space grid and deposit fourier amplitudes (squared) into power spectrum bins (i.e. shells in k-space)
     DO kz = 0, nmax/2-1
        DO ky = 0, nmax/2-1
            DO kx = 0, nmax/2-1
                
                IF((kx .EQ. 0) .AND. (ky .EQ. 0)  .AND. (kz .EQ. 0)) CYCLE

                k = k0 * SQRT(SNGL(kx**2 + ky**2 + kz**2))
                
                kvec(1) = k0*kx
                kvec(2) = k0*ky
                kvec(3) = k0*kz
                
                khat(:) = kvec(:)/k 

                ! paralell and perpendicular components of wave vector     
                kpar_hat(:) = B0(:) / bmag
                kpar(:) = (kvec(1)*kpar_hat(1)+kvec(2)*kpar_hat(2)+kvec(3)*kpar_hat(3)) * kpar_hat(:)
                kperp(:) = kvec(:) - kpar(:)
                
                kperp_hat(:) = 0.d0
                IF(SUM(kperp**2) .GT. 1.d-20) kperp_hat(:) = kperp(:) / SQRT(SUM(kperp**2))
               
                kpar_sqr  = SUM(kpar**2) 
                kpark_sqr  = SUM(kpar**2)/k**2 
                kperp_sqr = SUM(kperp**2)
                kperpk_sqr = SUM(kperp**2)/k**2              
              
                !*******************************************************
                ! check for parallel and perp propagation limiting cases
                !*******************************************************  
                IF(kperpk_sqr .LT. 1.d-6*kpark_sqr) THEN       ! parallel propagation
                    kperp(:) = 0.d0
                    kperp_hat(:) = 0.d0
                ELSE IF(kpark_sqr .LT. 1.d-6*kperpk_sqr) THEN  ! perp propagation
                    kpar(:) = 0.d0
                    kpar_hat(:) = 0.d0
                END IF
                
                ! compute displacement unit vectors for all three wave modes
                D = MAX(0.d0, (1.d0+alpha)**2 - 4.d0*alpha*kpar_sqr/(k**2))              
    
                zet_a(1) = kperp_hat(2)*B0_hat(3) - kperp_hat(3)*B0_hat(2)
                zet_a(2) = kperp_hat(3)*B0_hat(1) - kperp_hat(1)*B0_hat(3)
                zet_a(3) = kperp_hat(1)*B0_hat(2) - kperp_hat(2)*B0_hat(1)
                
                zet_f(:) = (-1.d0+alpha+SQRT(D))*kpar(:) + (1.d0+alpha+SQRT(D))*kperp(:) 
                IF(SUM(zet_f**2) .GT. 1.D-20) zet_f(:) = zet_f(:)/(SQRT(SUM(zet_f**2)))
                              
                zet_s(:) = (-1.d0+alpha-SQRT(D))*kpar(:) + (1.d0+alpha-SQRT(D))*kperp(:)  
                IF(SUM(zet_s**2) .GT. 1.D-20) zet_s(:) = zet_s(:)/(SQRT(SUM(zet_s**2)))
                
                
                tmp1 = khat(1)*zet_f(1) + khat(2)*zet_f(2) + khat(3)*zet_f(3) 
                tmp2 = khat(1)*zet_s(1) + khat(2)*zet_s(2) + khat(3)*zet_s(3) 
                delvk_f = 2.d0 * vf_k_im(kx,ky,kz)
                delvk_s = 2.d0 * vs_k_im(kx,ky,kz)
                
                ca = bmag                
                cf = SQRT(0.5d0 * (bmag**2) * (1.d0+alpha+SQRT(D)))
                cs = SQRT(MAX(0.d0, 0.5d0 * (bmag**2) * (1.d0+alpha-SQRT(D))))
                
                icf = 1.d0 / cf
                ics = 0.d0
                IF(cs .GT. 1.D-20) ics = 1.d0 / cs  ! for slow speed = 0, set ics to zero because magnetic, velocity and density amplitudes need to be zero 
                                                    ! delvk_s will be pretty close to zero anyways if that happens, so this is just for extra protection 
               
               
                !*******************************************************************
                GO TO 100
                IF(cf .LT. 1.D-25 .OR. cs .LT. 1.D-25) THEN
                    PRINT*,'WARNING!! Bad magnetosonic wave speed. cs, cf = ',cs,cf
                    PRINT*,' 1+alpha-SQRT(D)  = ',1.d0+alpha-SQRT(D)
                    PRINT*,'kpar = ',kpar
                    PRINT*,'mag(kpar/k)^2 = ',kpark_sqr
                    PRINT*,'kpar_hat = ',kpar_hat
                    PRINT*,'kperp = ',kperp
                    PRINT*,'mag(kperp/k)^2 = ',kperpk_sqr
                    PRINT*,'kperp_hat = ',kperp_hat
                    PRINT*,'kperpk_sqr/kpark_sqr = ',kperpk_sqr/kpark_sqr
                    PRINT*,'kpark_sqr/kperpk_sqr = ',kpark_sqr/kperpk_sqr
                    PRINT*,'mag(k)^2 = ',k**2
                    
                    PRINT*,'zet_a vector = ',zet_a
                    PRINT*,'mag(zet_a)^2 = ',SUM(zet_a**2)

                    PRINT*,'zet_f vector = ',zet_f                    
                    PRINT*,'mag(zet_f)^2 = ',SUM(zet_f**2)

                    PRINT*,'zet_s vector = ',zet_s
                    PRINT*,'mag(zet_s)^2 = ',SUM(zet_s**2)
         
                    PRINT*,'delvk_s = ',delvk_s
                    
                    PRINT*,'ics, icf = ',ics,icf
                    STOP
                END IF
                100 CONTINUE
                !*******************************************************************
                
                
                ! density fourier amplitudes (real part is zero)
                rhok_f(1) = 0.d0             
                rhok_s(1) = 0.d0
                rhok_f(2) = delvk_f * tmp1 * icf             
                rhok_s(2) = delvk_s * tmp2 * ics            
                
                
                tmp3 = SQRT((B0(2)*zet_f(3) - B0(3)*zet_f(2))**2 + (B0(3)*zet_f(1) - B0(1)*zet_f(3))**2 + &
                            (B0(1)*zet_f(2) - B0(2)*zet_f(1))**2)
                
                tmp4 = SQRT((B0(2)*zet_s(3) - B0(3)*zet_s(2))**2 + (B0(3)*zet_s(1) - B0(1)*zet_s(3))**2 + &
                            (B0(1)*zet_s(2) - B0(2)*zet_s(1))**2) 
                
                ! magnetic field fourier amplitudes
                bk_a(1) = va_k_re(kx,ky,kz) 
                bk_a(2) = va_k_im(kx,ky,kz)
                bk_f(1) = 0.d0
                bk_f(2) = delvk_f * tmp3 * icf 
                bk_s(1) = 0.d0 
                bk_s(2) = delvk_s * tmp4 * ics 
            
                ! compute amplitude squared            
                vsqr_a = 2.d0 *(va_k_re(kx,ky,kz)**2 + va_k_im(kx,ky,kz)**2)
                vsqr_f = 2.d0 *(vf_k_re(kx,ky,kz)**2 + vf_k_im(kx,ky,kz)**2)
                vsqr_s = 2.d0 *(vs_k_re(kx,ky,kz)**2 + vs_k_im(kx,ky,kz)**2)
                
                rhosqr_f =  2.d0 *(rhok_f(1)**2 + rhok_f(2)**2)
                rhosqr_s =  2.d0 *(rhok_s(1)**2 + rhok_s(2)**2)    
                
                bsqr_a =  2.d0 *(bk_a(1)**2 + bk_a(2)**2)
                bsqr_f =  2.d0 *(bk_f(1)**2 + bk_f(2)**2)
                bsqr_s =  2.d0 *(bk_s(1)**2 + bk_s(2)**2)
                
                ! deposit into k-bins
                ik = 1 + INT(k/dk)
                
                tmp1 = kmin + (ik-1) * dk 
                dv_shell = (FOURPI/3.d0)* ((tmp1+dk)**3 - tmp1**3)

                Pk_v(ik,1) = Pk_v(ik,1) + vsqr_a !* dv_shell
                Pk_v(ik,2) = Pk_v(ik,2) + vsqr_f !* dv_shell
                Pk_v(ik,3) = Pk_v(ik,3) + vsqr_s !* dv_shell
                Pk_b(ik,1) = Pk_b(ik,1) + bsqr_a !* dv_shell
                Pk_b(ik,2) = Pk_b(ik,2) + bsqr_f !* dv_shell
                Pk_b(ik,3) = Pk_b(ik,3) + bsqr_s !* dv_shell
                Pk_rho(ik,1) = Pk_rho(ik,1)  + rhosqr_f !* dv_shell
                Pk_rho(ik,2) = Pk_rho(ik,2)  + rhosqr_s !* dv_shell
    
            END DO
        END DO
    END DO    
        
    Pk_v = Pk_v * (dx**3)    
    Pk_b = Pk_b * (dx**3)    
    Pk_rho = Pk_rho * (dx**3)    
        
    ! dump power spectrum into file
    IF(.NOT. wavelet_reconstruction) Then
        filename = TRIM('Output/modes_Pk_dump=')//TRIM(uniti)//TRIM('.dat')
    ELSE        
        filename = TRIM('Output/modes_Pk_rec_dump=')//TRIM(uniti)//TRIM('.dat')
    END IF


    OPEN(UNIT=1, FILE=filename, FORM = 'UNFORMATTED', ACCESS = 'STREAM')
    
    DO ik = 1, kbins
        WRITE(1) (ik-1)*dk,Pk_v(ik,1),Pk_v(ik,2),Pk_v(ik,3),Pk_b(ik,1),Pk_b(ik,2),Pk_b(ik,3),Pk_rho(ik,1),Pk_rho(ik,2) 
    END DO
    
    CLOSE(UNIT=1)
    
    DEALLOCATE(Pk_v, Pk_b, Pk_rho)
 

END SUBROUTINE mode_decomposition_fourier



! MHD wave mode decomposition: First, we expand the velocity field in wavelets basis functions. Then we fourier transform each wavelet and decompose it
! into MHD modes. We sum up the contributions from all wavelets to obtain total power spectra for each MHD mode. 
SUBROUTINE mode_decomposition_wavelet(t_dump)

    INTEGER, INTENT(IN) :: t_dump 
    REAL*4 :: bmag, B0(3), alpha, D, zet_a(3), zet_f(3), zet_s(3), khat(3), &
              kdotB, kvec(3), kpar(3), kperp(3), kpar_hat(3), kperp_hat(3), tmp1, tmp2, &
              rhok_f(2), rhok_s(2), bk_a(2), bk_f(2), bk_s(2), &
              delvk_f, delvk_s, cf, cs, vsqr_a, vsqr_f, vsqr_s, rhosqr_f, rhosqr_s, bsqr_a, bsqr_f, bsqr_s  
    REAL*4 :: vwk_re(3), vwk_im(3), wcenter(3), vak_re, vak_im, vfk_re, vfk_im, vsk_re, vsk_im            
    INTEGER :: ix, iy, iz, kx, ky, kz, kbins, ik, counter
    REAL*4 :: kmin, kmax, dk, k, k0, dv_shell 
    CHARACTER(LEN=300) :: filename
    CHARACTER(LEN=6) :: uniti
    REAL*4, ALLOCATABLE :: Pk_v(:,:), Pk_b(:,:), Pk_rho(:,:) 


    IF(t_dump<10) THEN
        WRITE(uniti,'(I1.1)') t_dump
    ELSE IF(t_dump>=10 .and. t_dump<100) THEN
        WRITE(uniti,'(I2.2)') t_dump
    ELSE IF(t_dump>=100 .and. t_dump<1000) THEN
        WRITE (uniti,'(I3.3)') t_dump
    ELSE IF(t_dump>=1000 .and. t_dump<10000) THEN
        WRITE (uniti,'(I4.3)') t_dump
    ELSE IF(t_dump>=10000 .and. t_dump<100000) THEN
        WRITE (uniti,'(I5.3)') t_dump  
    END IF 
    
    ! set number of bins (i.e. k-shells)
    kbins = nmax
    
    ALLOCATE(Pk_v(1:kbins,3),Pk_b(1:kbins,3),Pk_rho(1:kbins,2))
    Pk_v = 0.d0
    Pk_b = 0.d0
    Pk_rho = 0.d0
    
    ! k band parameters
    k0 = TWOPI/Lx
    kmin = 0.d0
    kmax = SQRT(3.d0)*k0*(nmax/2)    
    dk = (kmax-kmin)/(kbins) ! uinform bin width


    !********************************************************
    ! Step 1: Compute wavelet co-efficients of velocity field
    !********************************************************
    PRINT*,'Computing wavelet transform of velocity field...'
    
    buffer(:,:,:) = fx(:,:,:,2) 
    CALL fwt_3d(buffer,vwx,1)
    buffer(:,:,:) = fx(:,:,:,3) 
    CALL fwt_3d(buffer,vwy,1)
    buffer(:,:,:) = fx(:,:,:,4) 
    CALL fwt_3d(buffer,vwz,1)

    filename = TRIM('Output/v_wt_dump=')//TRIM(uniti)//TRIM('.dat')
    OPEN(UNIT=1, FILE=filename, FORM = 'UNFORMATTED', ACCESS = 'STREAM')
    
    DO iz = 1, nranks_z*nz
        DO iy = 1, nranks_y*ny
            DO ix = 1, nranks_x*nx
                WRITE(1) vwx(ix,iy,iz),vwy(ix,iy,iz),vwz(ix,iy,iz)
            END DO
        END DO
    END DO
 
    CLOSE(UNIT=1)
    
    PRINT*,'Applying MHD mode decomposition to wavelet basis functions...'
   
    ! loop over wavelet basis functions
    counter = 0  
    
    !!$OMP PARALLEL PRIVATE(buffer,fwbk_re,fwbk_im,B0,bmag,alpha,k,kvec,kpar,kperp,khat,kpar_hat,kperp_hat,vwk_re,vwk_im,D,zet_a,zet_f,zet_s)  &   
    !!$OMP&  PRIVATE(vak_re,vak_im,vfk_re,vfk_im,vsk_re,vsk_im, tmp1,tmp2,delvk_f,delvk_s,cf,cs, rhok_f,rhok_s,bk_a,bk_f,bk_s, vsqr_a,vsqr_f,vsqr_s) &
    !!$OMP&  PRIVATE(rhosqr_f,rhosqr_s,bsqr_a,bsqr_f,bsqr_s, dv_shell, ix,iy,iz,kx,ky,kz,ik) 
    
    !!$OMP DO
    DO iz = 1, nranks_z*nz
        DO iy = 1, nranks_y*ny
            DO ix = 1, nranks_x*nx

                tmp1 = 100.0*SNGL(counter)/SNGL((nranks_x*nx)**3)
                IF(MOD(tmp1,1.0) .LT. 0.001) PRINT*,'% Complete = ',tmp1 
              
                !******************************************************************
                ! Step 2: Compute wavelet basis function and it's fourier transform 
                !******************************************************************
                
                buffer = 0.0
                CALL compute_basis_wavelet(ix,iy,iz,buffer)
                wcenter = MINLOC(buffer) ! find wavelet center
                CALL fft_3d(buffer, fwbk_re, fwbk_im)               

                !**********************************************************************************************************
                ! Step 3: Apply MHD mode decomposition to this wavelet and add up their contributions to the power spectrum
                !*********************************************************************************************************

                ! compute local mean magnetic field (i.e. near the wavelet center)
                ! Need to fiugure out a good way to compute this. Maybe average it over the wavelet function norm/magnitude?? 
                
                !B0(:) = 0.d0
                !DO k = 1, nranks_z*nz
                !    DO j = 1, nranks_y*ny
                !        DO i = 1, nranks_x*nx

                !          B0(1) = B0(1) + fx(i,j,k,5)
                !          B0(2) = B0(2) + fx(i,j,k,6)
                !          B0(3) = B0(3) + fx(i,j,k,7)
                        
                !        END DO
                !    END DO
                !END DO
                !B0(:) = B0(:) /(nmax**3) 
                
                ! mean B-field hard coded for now...
                B0 = (/ 1.0, 0.0, 0.0 /)
                bmag = SQRT(SUM(B0**2))
                alpha = 1.d0/(bmag**2)

                ! loop over k space grid
                DO kz = 0, nmax/2-1
                    DO ky = 0, nmax/2-1
                        DO kx = 0, nmax/2-1

                            IF((kx .EQ. 0) .AND. (ky .EQ. 0)  .AND. (kz .EQ. 0)) CYCLE

                            k = k0 * SQRT(SNGL(kx**2 + ky**2 + kz**2))
                            
                            kvec(1) = k0*kx
                            kvec(2) = k0*ky
                            kvec(3) = k0*kz                            
                            khat(:) = kvec(:)/k 

                            ! paralell and perpendicular components of wave vector     
                            kpar_hat(:) = B0(:) / bmag
                            kpar(:) = (kvec(1)*kpar_hat(1)+kvec(2)*kpar_hat(2)+kvec(3)*kpar_hat(3)) * kpar_hat(:)
                            kperp(:) = kvec(:) - kpar(:)
                            
                            kperp_hat(:) = 0.d0
                            IF(SUM(kperp**2) .GT. 0.d0) kperp_hat(:) = kperp(:) / SQRT(SUM(kperp**2))
                            
                            ! velocity fourier components of this wavelet
                            vwk_re(1) = vwx(ix,iy,iz) * fwbk_re(kx,ky,kz)
                            vwk_re(2) = vwy(ix,iy,iz) * fwbk_re(kx,ky,kz)
                            vwk_re(3) = vwz(ix,iy,iz) * fwbk_re(kx,ky,kz)
                           
                            vwk_im(1) = vwx(ix,iy,iz) * fwbk_im(kx,ky,kz)
                            vwk_im(2) = vwy(ix,iy,iz) * fwbk_im(kx,ky,kz)
                            vwk_im(3) = vwz(ix,iy,iz) * fwbk_im(kx,ky,kz)                  
                            
                       
                            ! compute displacement unit vectors for all three wave modes
                            D = MAX(0.d0, (1.d0+alpha)**2 - 4.d0*alpha*SUM(kpar**2)/(k**2))              
                
                            zet_a(1) = kperp_hat(2)*kpar_hat(3) - kperp_hat(3)*kpar_hat(2)
                            zet_a(2) = kperp_hat(3)*kpar_hat(1) - kperp_hat(1)*kpar_hat(3)
                            zet_a(3) = kperp_hat(1)*kpar_hat(2) - kperp_hat(2)*kpar_hat(1)
                            
                            zet_f(:) = (-1.d0+alpha+SQRT(D))*kpar(:) + (1.d0+alpha+SQRT(D))*kperp(:) 
                            IF(SUM(zet_f**2) .GT. 0) zet_f(:) = zet_f(:)/(SQRT(SUM(zet_f**2)))
                                          
                            zet_s(:) = (-1.d0+alpha-SQRT(D))*kpar(:) + (1.d0+alpha-SQRT(D))*kperp(:)  
                            IF(SUM(zet_s**2) .GT. 0) zet_s(:) = zet_s(:)/(SQRT(SUM(zet_s**2)))
                            
                            ! project velocity fourier components onto displacement unit vectors of each mode 
                            
                            ! alfven mode
                            vak_re = vwk_re(1)*zet_a(1) + vwk_re(2)*zet_a(2) + vwk_re(3)*zet_a(3)
                            vak_im = vwk_im(1)*zet_a(1) + vwk_im(2)*zet_a(2) + vwk_im(3)*zet_a(3)      
                            
                            ! fast mode
                            vfk_re = vwk_re(1)*zet_f(1) + vwk_re(2)*zet_f(2) + vwk_re(3)*zet_f(3)  
                            vfk_im = vwk_im(1)*zet_f(1) + vwk_im(2)*zet_f(2) + vwk_im(3)*zet_f(3) 
                            
                            ! slow mode
                            vsk_re = vwk_re(1)*zet_s(1) + vwk_re(2)*zet_s(2) + vwk_re(3)*zet_s(3)                 
                            vsk_im = vwk_im(1)*zet_s(1) + vwk_im(2)*zet_s(2) + vwk_im(3)*zet_s(3)                 
                           
                           
                            tmp1 = khat(1)*zet_f(1) + khat(2)*zet_f(2) + khat(3)*zet_f(3) 
                            tmp2 = khat(1)*zet_s(1) + khat(2)*zet_s(2) + khat(3)*zet_s(3) 
                            delvk_f = 2.d0 * vfk_im
                            delvk_s = 2.d0 * vsk_im
                            
                            cf = SQRT(0.5d0 * (bmag**2) * (1.d0+alpha+SQRT(D)))
                            cs = SQRT(MAX(0.d0, 0.5d0 * (bmag**2) * (1.d0+alpha-SQRT(D))))
                            
                            ! density fourier amplitudes (real part is zero)
                            rhok_f(1) = 0.d0             
                            rhok_s(1) = 0.d0
                            rhok_f(2) = delvk_f * tmp1 /cf             
                            rhok_s(2) = delvk_s * tmp2 /cs            
                            
                            
                            tmp1 = SQRT((B0(2)*zet_f(3) - B0(3)*zet_f(2))**2 + (B0(3)*zet_f(1) - B0(1)*zet_f(3))**2 + &
                                        (B0(1)*zet_f(2) - B0(2)*zet_f(1))**2)
                            
                            tmp2 = SQRT((B0(2)*zet_s(3) - B0(3)*zet_s(2))**2 + (B0(3)*zet_s(1) - B0(1)*zet_s(3))**2 + &
                                        (B0(1)*zet_s(2) - B0(2)*zet_s(1))**2) 
                            
                            ! magnetic field fourier amplitudes
                            bk_a(1) = vak_re
                            bk_a(2) = vak_im
                            bk_f(1) = 0.d0
                            bk_f(2) = delvk_f * tmp1 /cf 
                            bk_s(1) = 0.d0 
                            bk_s(2) = delvk_s * tmp2 /cs 
                        

                            ! compute amplitude squared            
                            vsqr_a = 2.d0 *(vak_re**2 + vak_im**2)
                            vsqr_f = 2.d0 *(vfk_re**2 + vfk_im**2)
                            vsqr_s = 2.d0 *(vsk_re**2 + vsk_im**2)
                            
                            rhosqr_f =  2.d0 *(rhok_f(1)**2 + rhok_f(2)**2)
                            rhosqr_s =  2.d0 *(rhok_s(1)**2 + rhok_s(2)**2)
                            
                            
                            bsqr_a =  2.d0 *(bk_a(1)**2 + bk_a(2)**2)
                            bsqr_f =  2.d0 *(bk_f(1)**2 + bk_f(2)**2)
                            bsqr_s =  2.d0 *(bk_s(1)**2 + bk_s(2)**2)
                            
                            ! deposit fourier amplitudes (squared) into power spectrum bins (i.e. shells in k-space)
                            ik = 1 + INT(k/dk)
                            
                            tmp1 = kmin + (ik-1) * dk 
                            dv_shell = (FOURPI/3.d0)* ((tmp1+dk)**3 - tmp1**3)

                            !!$OMP ATOMIC UPDATE
                            Pk_v(ik,1) = Pk_v(ik,1) + vsqr_a !* dv_shell
                            !!$OMP END ATOMIC
                            
                            !!$OMP ATOMIC UPDATE
                            Pk_v(ik,2) = Pk_v(ik,2) + vsqr_f !* dv_shell
                            !!$OMP END ATOMIC
                            
                            !!$OMP ATOMIC UPDATE
                            Pk_v(ik,3) = Pk_v(ik,3) + vsqr_s !* dv_shell
                            !!$OMP END ATOMIC
                            
                            !!$OMP ATOMIC UPDATE
                            Pk_b(ik,1) = Pk_b(ik,1) + bsqr_a !* dv_shell
                            !!$OMP END ATOMIC
                            
                            !!$OMP ATOMIC UPDATE                            
                            Pk_b(ik,2) = Pk_b(ik,2) + bsqr_f !* dv_shell
                            !!$OMP END ATOMIC
                            
                            !!$OMP ATOMIC UPDATE
                            Pk_b(ik,3) = Pk_b(ik,3) + bsqr_s !* dv_shell
                            !!$OMP END ATOMIC
                            
                            !!$OMP ATOMIC UPDATE
                            Pk_rho(ik,1) = Pk_rho(ik,1)  + rhosqr_f !* dv_shell
                            !!$OMP END ATOMIC
                            
                            !!$OMP ATOMIC UPDATE
                            Pk_rho(ik,2) = Pk_rho(ik,2)  + rhosqr_s !* dv_shell
                            !!$OMP END ATOMIC
                                                        
                        END DO
                    END DO
                END DO

                !!$OMP ATOMIC UPDATE
                counter = counter + 1
                !!$OMP END ATOMIC

            END DO
        END DO
    END DO
    !!$OMP END PARALLEL    

    ! apply normalization       
    Pk_v = Pk_v * (dx**3)    
    Pk_b = Pk_b * (dx**3)    
    Pk_rho = Pk_rho * (dx**3)    
        
    ! dump power spectrum into file
    filename = TRIM('Output/modes_Pk_dump=')//TRIM(uniti)//TRIM('.dat')

    OPEN(UNIT=1, FILE=filename, FORM = 'UNFORMATTED', ACCESS = 'STREAM')
    
    DO ik = 1, kbins
        WRITE(1) (ik-1)*dk,Pk_v(ik,1),Pk_v(ik,2),Pk_v(ik,3),Pk_b(ik,1),Pk_b(ik,2),Pk_b(ik,3),Pk_rho(ik,1),Pk_rho(ik,2) 
    END DO
    
    CLOSE(UNIT=1)
    
    DEALLOCATE(Pk_v, Pk_b, Pk_rho)
 

END SUBROUTINE mode_decomposition_wavelet



SUBROUTINE mode_decomposition_tracer_wavelets(t_dump)

    INTEGER, INTENT(IN) :: t_dump 
    REAL*4 :: bmag, B0(3), B0_hat(3), alpha, D, zet_a(3), zet_f(3), zet_s(3), khat(3), &
              kdotB, kvec(3), kpar(3), kperp(3), kpar_hat(3), kperp_hat(3), tmp1, tmp2, tmp3, tmp4, rhok_f(2), rhok_s(2), bk_a(2), bk_f(2), bk_s(2), &
              delvk_f, delvk_s, ca, cf, cs, icf, ics, vsqr_a, vsqr_f, vsqr_s, rhosqr_f, rhosqr_s, bsqr_a, bsqr_f, bsqr_s, sound_speed  
    INTEGER :: kx, ky, kz, kbins, ik, nranks_tot, irank, ri, rj, rk, i, j, ii, ix, iy, iz, indx
    REAL*4 :: kmin, kmax, dk, k, kpar_sqr, kperp_sqr, kpark_sqr, kperpk_sqr, k0, dv_shell, cmplt 
    CHARACTER(LEN=300) :: filename
    CHARACTER(LEN=6) :: uniti
    REAL*4, ALLOCATABLE :: Pk_v(:,:), Pk_b(:,:), Pk_rho(:,:) 
    REAL(4) :: vw(nx,ny,nz), buff_re(nx,ny,nz),buff_im(nx,ny,nz), rho0
    
    LOGICAL :: printd = .FALSE.

    IF(t_dump<10) THEN
        WRITE(uniti,'(I1.1)') t_dump
    ELSE IF(t_dump>=10 .and. t_dump<100) THEN
        WRITE(uniti,'(I2.2)') t_dump
    ELSE IF(t_dump>=100 .and. t_dump<1000) THEN
        WRITE (uniti,'(I3.3)') t_dump
    ELSE IF(t_dump>=1000 .and. t_dump<10000) THEN
        WRITE (uniti,'(I4.3)') t_dump
    ELSE IF(t_dump>=10000 .and. t_dump<100000) THEN
        WRITE (uniti,'(I5.3)') t_dump  
    END IF 
    
    ! set number of bins (i.e. k-shells)
    kbins = nx
    
    ALLOCATE(Pk_v(1:kbins,3),Pk_b(1:kbins,3),Pk_rho(1:kbins,2))
    Pk_v = 0.d0
    Pk_b = 0.d0
    Pk_rho = 0.d0
    
    ! k band parameters
    k0 = TWOPI/(dx*nx)
    kmin = 0.d0
    kmax = SQRT(3.d0)*k0*(nx/2)    
    dk = (kmax-kmin)/(kbins) ! uinform bin width

    ! set isothermal sound speed
    sound_speed = 1.D0
    
    ! resize FFT grid to patch size
    CALL fft_init(nx,ny,nz)

    ! resize fwt grid to patch size
    CALL fwt_init(nx,ny,nz)

    cmplt = 0.D0    
    nranks_tot = nranks_x * nranks_y * nranks_z
    
    ! loop over MPI ranks
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(irank,ri,rj,rk,i,j,ii,ix,iy,iz,vw,buff_re,buff_im,rho0,B0,B0_hat,bmag,kx,ky,kz) & 
    !$OMP& PRIVATE(k,kvec,khat,kpar,kpar_hat,kperp,kperp_hat,kpar_sqr,kperp_sqr,kpark_sqr,kperpk_sqr,alpha,D,zet_a,zet_f,zet_s) &
    !$OMP& PRIVATE(tmp1,tmp2,tmp3,tmp4,delvk_f,delvk_s,ca,cf,cs,icf,ics,bk_a,bk_f,bk_s) &
    !$OMP& PRIVATE(vsqr_a,vsqr_f,vsqr_s,rhosqr_f,rhosqr_s,bsqr_a,bsqr_f,bsqr_s,ik) &
    !$OMP& SHARED(Pk_v,Pk_b,Pk_rho) 
    
    DO irank = 1, nranks_tot
    
        CALL collapse_index(irank,ri,rj,rk)

        PRINT*,''
        PRINT*,'Performing wavelet mode decomp for rank: ',ri,rj,rk
  
        ! loop over tracer particles
        DO i = 1, dat_rec(ri,rj,rk)%nparticles
        
            !PRINT*,'Particle#, rank: ',i,ri,rj,rk
            
            ! loop over wavelets
            DO j = 1, nwavelets_per_scalar
      
                !PRINT*,'wavelet#, rank: ',j,ri,rj,rk

                IF(INT(dat_rec(ri,rj,rk)%tracer_wavelets(i,j,1)) .LE. 0) CYCLE
        
                ! flattened wavelet index
                ii = dat_rec(ri,rj,rk)%tracer_wavelets(i,j,1)
                
                ! collapse flattened wavelet index
                ix = 1 + MOD(ii-1,nx)
                iy = 1 + MOD(ii-1,nx*ny)/nx 
                iz = 1 + (ii-1) / (nx*ny) 
               
                ! generate the basis wavelet
                !PRINT*,'Computing basis wavelet, rank: ',j,ri,rj,rk
                CALL compute_basis_wavelet(ix,iy,iz,vw)
  
                ! get local average density and B field 
                rho0    = dat_rec(ri,rj,rk)%wavelet_avg_rhoB(dat_rec(ri,rj,rk)%wavelet_lookup_table(ii),2)   
                B0(1:3) = dat_rec(ri,rj,rk)%wavelet_avg_rhoB(dat_rec(ri,rj,rk)%wavelet_lookup_table(ii),3:5) 
                bmag = SQRT(SUM(B0**2))
         
                !IF(printd) PRINT*,'B0 = , rank: ',B0,ri,rj,rk
                !PRINT*,'B0, rho0 = ',B0,rho0
                !PRINT*,'bmag, rho0 = ',bmag,rho0


                B0_hat(:) = B0(:) / bmag
                alpha = 1.d0/SUM(B0**2)


                ! compute FFT of this wavelet
                !PRINT*,'Computing FFT, rank: ',ri,rj,rk
                CALL fft_3d(vw, buff_re, buff_im) 

                !PRINT*,'FFT DONE! rank: ',ri,rj,rk

                dat_rec(ri,rj,rk)%wavelet_vxk_re(:,:,:) = dat_rec(ri,rj,rk)%tracer_wavelets(i,j,2) * buff_re(:,:,:)
                dat_rec(ri,rj,rk)%wavelet_vxk_im(:,:,:) = dat_rec(ri,rj,rk)%tracer_wavelets(i,j,2) * buff_im(:,:,:)
                
                dat_rec(ri,rj,rk)%wavelet_vyk_re(:,:,:) = dat_rec(ri,rj,rk)%tracer_wavelets(i,j,3) * buff_re(:,:,:)
                dat_rec(ri,rj,rk)%wavelet_vyk_im(:,:,:) = dat_rec(ri,rj,rk)%tracer_wavelets(i,j,3) * buff_im(:,:,:)
                
                dat_rec(ri,rj,rk)%wavelet_vzk_re(:,:,:) = dat_rec(ri,rj,rk)%tracer_wavelets(i,j,4) * buff_re(:,:,:)
                dat_rec(ri,rj,rk)%wavelet_vzk_im(:,:,:) = dat_rec(ri,rj,rk)%tracer_wavelets(i,j,4) * buff_im(:,:,:)
    
                !**************************
                ! do the mode decomposition 
                !**************************

                !PRINT*,'Computing Modes... rank: ',ri,rj,rk

                ! loop over k space grid
                DO kz = 0, nz/2-1
                    DO ky = 0, ny/2-1
                        DO kx = 0, nx/2-1

                            !IF(kx .EQ. 1 .AND. ky .EQ. 0 .AND. kz .EQ. 0 .AND. i .EQ. 22 .AND. j .EQ. 27) THEN
                            !    printd = .TRUE.
                            !ELSE
                            !    printd = .FALSE.
                            !END IF    

                            IF(printd) PRINT*,'################'    
                            IF(printd) PRINT*,'kx, ky, kz, i, j = ',kx, ky, kz, i, j    
                            IF(printd) PRINT*,'################'    

                            IF((kx .EQ. 0) .AND. (ky .EQ. 0)  .AND. (kz .EQ. 0)) CYCLE

                            k = k0 * SQRT(SNGL(kx**2 + ky**2 + kz**2))

                            IF(printd) PRINT*,'k, rank: ',k,ri,rj,rk
                            IF(printd) PRINT*,'B0, rank: ',B0,ri,rj,rk
                            IF(printd) PRINT*,'bmag, rank: ',bmag,ri,rj,rk

                            
                            kvec(1) = k0*kx
                            kvec(2) = k0*ky
                            kvec(3) = k0*kz
                            khat(:) = kvec(:)/k 

                            IF(printd) PRINT*,'kvec, rank: ',kvec,ri,rj,rk
                            IF(printd) PRINT*,'khat, rank: ',khat,ri,rj,rk


                            ! paralell and perpendicular components of wave vector     
                            kpar_hat(:) = B0(:) / bmag
                            kpar(:) = (kvec(1)*kpar_hat(1)+kvec(2)*kpar_hat(2)+kvec(3)*kpar_hat(3)) * kpar_hat(:)
                            kperp(:) = kvec(:) - kpar(:)
                            
                            kperp_hat(:) = 0.d0
                            IF(SUM(kperp**2) .GT. 1.d-15) kperp_hat(:) = kperp(:) / SQRT(SUM(kperp**2))
                           
                            IF(printd) PRINT*,'kpar, rank: ',kpar,ri,rj,rk
                            IF(printd) PRINT*,'kperp, rank: ',kperp,ri,rj,rk
                            
                            kpar_sqr   = SUM(kpar**2) 
                            kpark_sqr  = SUM(kpar**2)/k**2 
                            kperp_sqr  = SUM(kperp**2)
                            kperpk_sqr = SUM(kperp**2)/k**2              
                          
                            IF(printd) PRINT*,'kpar_sqr, rank: ',kpar_sqr,ri,rj,rk
                            IF(printd) PRINT*,'kperp_sqr, rank: ',kperp_sqr,ri,rj,rk
                            IF(printd) PRINT*,'kpark_sqr, rank: ',kpark_sqr,ri,rj,rk
                            IF(printd) PRINT*,'kperpk_sqr, rank: ',kperpk_sqr,ri,rj,rk
                          
                            !*******************************************************
                            ! check for parallel and perp propagation limiting cases
                            !*******************************************************  
                            IF(kperpk_sqr .LT. 1.d-6*kpark_sqr) THEN       ! parallel propagation
                                kperp(:) = 0.d0
                                kperp_hat(:) = 0.d0
                            ELSE IF(kpark_sqr .LT. 1.d-6*kperpk_sqr) THEN  ! perp propagation
                                kpar(:) = 0.d0
                                kpar_hat(:) = 0.d0
                            END IF
                                                        
                            ! compute displacement unit vectors for all three wave modes
                            D = MAX(0.d0, (1.d0+alpha)**2 - 4.d0*alpha*SUM(kpar**2)/(k**2))              

                            zet_a(1) = kperp_hat(2)*B0_hat(3) - kperp_hat(3)*B0_hat(2)
                            zet_a(2) = kperp_hat(3)*B0_hat(1) - kperp_hat(1)*B0_hat(3)
                            zet_a(3) = kperp_hat(1)*B0_hat(2) - kperp_hat(2)*B0_hat(1)
                                            
                            zet_f(:) = (-1.d0+alpha+SQRT(D))*kpar(:) + (1.d0+alpha+SQRT(D))*kperp(:) 
                            IF(SUM(zet_f**2) .GT. 1.d-15) zet_f(:) = zet_f(:)/(SQRT(SUM(zet_f**2)))
                                          
                            zet_s(:) = (-1.d0+alpha-SQRT(D))*kpar(:) + (1.d0+alpha-SQRT(D))*kperp(:)  
                            IF(SUM(zet_s**2) .GT. 1.d-15) zet_s(:) = zet_s(:)/(SQRT(SUM(zet_s**2)))

                            IF(printd) PRINT*,'zet_a, rank: ',zet_a,ri,rj,rk
                            IF(printd) PRINT*,'zet_f, rank: ',zet_f,ri,rj,rk
                            IF(printd) PRINT*,'zet_s, rank: ',zet_s,ri,rj,rk

                            !********************************************************************************                
                            ! project velocity fourier components onto displacement unit vectors of each mode
                            !******************************************************************************** 

                            ! alfven mode
                            dat_rec(ri,rj,rk)%wavelet_va_k_re(kx,ky,kz) = dat_rec(ri,rj,rk)%wavelet_vxk_re(kx,ky,kz)*zet_a(1) + &
                                                                          dat_rec(ri,rj,rk)%wavelet_vyk_re(kx,ky,kz)*zet_a(2) + &
                                                                          dat_rec(ri,rj,rk)%wavelet_vzk_re(kx,ky,kz)*zet_a(3)
                                                                       
                            dat_rec(ri,rj,rk)%wavelet_va_k_im(kx,ky,kz) = dat_rec(ri,rj,rk)%wavelet_vxk_im(kx,ky,kz)*zet_a(1) + &
                                                                          dat_rec(ri,rj,rk)%wavelet_vyk_im(kx,ky,kz)*zet_a(2) + &
                                                                          dat_rec(ri,rj,rk)%wavelet_vzk_im(kx,ky,kz)*zet_a(3)      
                            
                            ! fast mode
                            dat_rec(ri,rj,rk)%wavelet_vf_k_re(kx,ky,kz) = dat_rec(ri,rj,rk)%wavelet_vxk_re(kx,ky,kz)*zet_f(1) + &
                                                                          dat_rec(ri,rj,rk)%wavelet_vyk_re(kx,ky,kz)*zet_f(2) + &
                                                                          dat_rec(ri,rj,rk)%wavelet_vzk_re(kx,ky,kz)*zet_f(3) 
                                                                       
                            dat_rec(ri,rj,rk)%wavelet_vf_k_im(kx,ky,kz) = dat_rec(ri,rj,rk)%wavelet_vxk_im(kx,ky,kz)*zet_f(1) + &
                                                                          dat_rec(ri,rj,rk)%wavelet_vyk_im(kx,ky,kz)*zet_f(2) + &
                                                                          dat_rec(ri,rj,rk)%wavelet_vzk_im(kx,ky,kz)*zet_f(3) 
                            
                            ! slow mode
                            dat_rec(ri,rj,rk)%wavelet_vs_k_re(kx,ky,kz) = dat_rec(ri,rj,rk)%wavelet_vxk_re(kx,ky,kz)*zet_s(1) + &
                                                                          dat_rec(ri,rj,rk)%wavelet_vyk_re(kx,ky,kz)*zet_s(2) + &
                                                                          dat_rec(ri,rj,rk)%wavelet_vzk_re(kx,ky,kz)*zet_s(3)  
                                                                       
                            dat_rec(ri,rj,rk)%wavelet_vs_k_im(kx,ky,kz) = dat_rec(ri,rj,rk)%wavelet_vxk_im(kx,ky,kz)*zet_s(1) + &
                                                                          dat_rec(ri,rj,rk)%wavelet_vyk_im(kx,ky,kz)*zet_s(2) + &
                                                                          dat_rec(ri,rj,rk)%wavelet_vzk_im(kx,ky,kz)*zet_s(3)                 
                           
                           
                            tmp1 = khat(1)*zet_f(1) + khat(2)*zet_f(2) + khat(3)*zet_f(3) 
                            tmp2 = khat(1)*zet_s(1) + khat(2)*zet_s(2) + khat(3)*zet_s(3) 
                            delvk_f = 2.d0 * dat_rec(ri,rj,rk)%wavelet_vf_k_im(kx,ky,kz)
                            delvk_s = 2.d0 * dat_rec(ri,rj,rk)%wavelet_vs_k_im(kx,ky,kz)
                            
                            
                            IF(printd) PRINT*,'tmp1, rank: ',tmp1,ri,rj,rk
                            IF(printd) PRINT*,'tmp2, rank: ',tmp2,ri,rj,rk
                            IF(printd) PRINT*,'delvk_f, rank: ',delvk_f,ri,rj,rk
                            IF(printd) PRINT*,'delvk_s, rank: ',delvk_s,ri,rj,rk

                            ca = bmag                
                            cf = SQRT(0.5d0 * (bmag**2) * (1.d0+alpha+SQRT(D)))
                            cs = SQRT(MAX(0.d0, 0.5d0 * (bmag**2) * (1.d0+alpha-SQRT(D))))


                            IF(printd) PRINT*,'ca, rank: ',ca,ri,rj,rk
                            IF(printd) PRINT*,'cf, rank: ',cf,ri,rj,rk
                            IF(printd) PRINT*,'cs, rank: ',cs,ri,rj,rk
                            
                            icf = 1.d0 / cf
                            ics = 0.d0
                            IF(cs .GT. 1.D-15) ics = 1.d0 / cs  ! for slow speed = 0, set ics to zero because magnetic, velocity and density amplitudes need to be zero 
                                                                ! delvk_s will be pretty close to zero anyways if that happens, so this is just for extra protection                 
              
                            IF(printd) PRINT*,'icf, rank: ',icf,ri,rj,rk
                            IF(printd) PRINT*,'ics, rank: ',ics,ri,rj,rk
              
                            ! density fourier amplitudes (real part is zero)
                            rhok_f(1) = 0.d0             
                            rhok_s(1) = 0.d0
                            rhok_f(2) = delvk_f * tmp1 * icf             
                            rhok_s(2) = delvk_s * tmp2 * ics            
                            
                            
                            IF(printd) PRINT*,'rhok_f, rank: ',rhok_f,ri,rj,rk
                            IF(printd) PRINT*,'rhok_s, rank: ',rhok_s,ri,rj,rk
                            
                            
                            tmp3 = SQRT((B0(2)*zet_f(3) - B0(3)*zet_f(2))**2 + (B0(3)*zet_f(1) - B0(1)*zet_f(3))**2 + &
                                        (B0(1)*zet_f(2) - B0(2)*zet_f(1))**2)
                            
                            tmp4 = SQRT((B0(2)*zet_s(3) - B0(3)*zet_s(2))**2 + (B0(3)*zet_s(1) - B0(1)*zet_s(3))**2 + &
                                        (B0(1)*zet_s(2) - B0(2)*zet_s(1))**2) 
                            
                            
                            IF(printd) PRINT*,'tmp3, rank: ',tmp3,ri,rj,rk
                            IF(printd) PRINT*,'tmp4, rank: ',tmp4,ri,rj,rk
                            
                            ! magnetic field fourier amplitudes
                            bk_a(1) = dat_rec(ri,rj,rk)%wavelet_va_k_re(kx,ky,kz) 
                            bk_a(2) = dat_rec(ri,rj,rk)%wavelet_va_k_im(kx,ky,kz)
                            bk_f(1) = 0.d0
                            bk_f(2) = delvk_f * tmp3 * icf 
                            bk_s(1) = 0.d0 
                            bk_s(2) = delvk_s * tmp4 * ics 
           
                            IF(printd) PRINT*,'bk_a, rank: ',bk_a,ri,rj,rk
                            IF(printd) PRINT*,'bk_f, rank: ',bk_f,ri,rj,rk
                            IF(printd) PRINT*,'bk_s, rank: ',bk_s,ri,rj,rk
           
                            ! compute amplitude squared            
                            vsqr_a = 2.d0 *(dat_rec(ri,rj,rk)%wavelet_va_k_re(kx,ky,kz)**2 + dat_rec(ri,rj,rk)%wavelet_va_k_im(kx,ky,kz)**2)
                            vsqr_f = 2.d0 *(dat_rec(ri,rj,rk)%wavelet_vf_k_re(kx,ky,kz)**2 + dat_rec(ri,rj,rk)%wavelet_vf_k_im(kx,ky,kz)**2)
                            vsqr_s = 2.d0 *(dat_rec(ri,rj,rk)%wavelet_vs_k_re(kx,ky,kz)**2 + dat_rec(ri,rj,rk)%wavelet_vs_k_im(kx,ky,kz)**2)
                            
                            IF(printd) PRINT*,'vsqr_a, rank: ',vsqr_a,ri,rj,rk
                            IF(printd) PRINT*,'vsqr_f, rank: ',vsqr_f,ri,rj,rk
                            IF(printd) PRINT*,'vsqr_s, rank: ',vsqr_s,ri,rj,rk


                            rhosqr_f =  2.d0 *(rhok_f(1)**2 + rhok_f(2)**2)
                            rhosqr_s =  2.d0 *(rhok_s(1)**2 + rhok_s(2)**2)    
                            
                            IF(printd) PRINT*,'rhosqr_f, rank: ',rhosqr_f,ri,rj,rk
                            IF(printd) PRINT*,'rhosqr_s, rank: ',rhosqr_s,ri,rj,rk


                            bsqr_a =  2.d0 *(bk_a(1)**2 + bk_a(2)**2)
                            bsqr_f =  2.d0 *(bk_f(1)**2 + bk_f(2)**2)
                            bsqr_s =  2.d0 *(bk_s(1)**2 + bk_s(2)**2)

                            IF(printd) PRINT*,'bsqr_a, rank: ',bsqr_a,ri,rj,rk
                            IF(printd) PRINT*,'bsqr_f, rank: ',bsqr_f,ri,rj,rk
                            IF(printd) PRINT*,'bsqr_s, rank: ',bsqr_s,ri,rj,rk

                            
                            ! deposit into k-bins
                            ik = 1 + INT(k/dk)
                            
                            IF(printd) PRINT*,'ik, rank: ',ik,ri,rj,rk
                            

                            !tmp1 = kmin + (ik-1) * dk 
                            !dv_shell = (FOURPI/3.d0)* ((tmp1+dk)**3 - tmp1**3)

                            dat_rec(ri,rj,rk)%wavelet_Pk_v(ik,1) = dat_rec(ri,rj,rk)%wavelet_Pk_v(ik,1) + vsqr_a !* dv_shell
                            dat_rec(ri,rj,rk)%wavelet_Pk_v(ik,2) = dat_rec(ri,rj,rk)%wavelet_Pk_v(ik,2) + vsqr_f !* dv_shell
                            dat_rec(ri,rj,rk)%wavelet_Pk_v(ik,3) = dat_rec(ri,rj,rk)%wavelet_Pk_v(ik,3) + vsqr_s !* dv_shell
                            dat_rec(ri,rj,rk)%wavelet_Pk_b(ik,1) = dat_rec(ri,rj,rk)%wavelet_Pk_b(ik,1) + bsqr_a !* dv_shell
                            dat_rec(ri,rj,rk)%wavelet_Pk_b(ik,2) = dat_rec(ri,rj,rk)%wavelet_Pk_b(ik,2) + bsqr_f !* dv_shell
                            dat_rec(ri,rj,rk)%wavelet_Pk_b(ik,3) = dat_rec(ri,rj,rk)%wavelet_Pk_b(ik,3) + bsqr_s !* dv_shell
                            dat_rec(ri,rj,rk)%wavelet_Pk_rho(ik,1) = dat_rec(ri,rj,rk)%wavelet_Pk_rho(ik,1)  + rhosqr_f !* dv_shell
                            dat_rec(ri,rj,rk)%wavelet_Pk_rho(ik,2) = dat_rec(ri,rj,rk)%wavelet_Pk_rho(ik,2)  + rhosqr_s !* dv_shell
                           
                        END DO
                    END DO
                END DO   
                !PRINT*,'Done computing Modes. rank: ',ri,rj,rk
                
            END DO   
        END DO       
  
        !$OMP ATOMIC
        cmplt = cmplt + 1
   
        PRINT*,'Done wavelet mode decomp for rank: ',ri,rj,rk
        PRINT*,'Percent complete = ',100.0*REAL(cmplt)/REAL(nranks_tot)
                
    END DO
    !$OMP END PARALLEL DO
   
    ! sum up contirubutions from all mpi ranks
    DO irank = 1, nranks_tot
        CALL collapse_index(irank,ri,rj,rk)
        Pk_v = Pk_v + dat_rec(ri,rj,rk)%wavelet_Pk_v
        Pk_b = Pk_b + dat_rec(ri,rj,rk)%wavelet_Pk_b
        Pk_rho = Pk_rho + dat_rec(ri,rj,rk)%wavelet_Pk_rho
    END DO
   
   
    Pk_v = Pk_v * (dx**3)    
    Pk_b = Pk_b * (dx**3)    
    Pk_rho = Pk_rho * (dx**3)    
        
    ! dump power spectrum into file
    filename = TRIM('Output/modes_wavelet_Pk_dump=')//TRIM(uniti)//TRIM('.dat')

    PRINT*,'Saving to file...'
    OPEN(UNIT=11, FILE=filename, FORM = 'UNFORMATTED', ACCESS = 'STREAM')
    
    DO ik = 1, kbins
        WRITE(11) (ik-1)*dk,Pk_v(ik,1),Pk_v(ik,2),Pk_v(ik,3),Pk_b(ik,1),Pk_b(ik,2),Pk_b(ik,3),Pk_rho(ik,1),Pk_rho(ik,2) 
    END DO
    
    CLOSE(UNIT=11)
    

    DEALLOCATE(Pk_v, Pk_b, Pk_rho)
 

END SUBROUTINE mode_decomposition_tracer_wavelets



SUBROUTINE compute_pk(t_dump, fieldtype, fxk_re, fxk_im, fyk_re, fyk_im, fzk_re, fzk_im)

    INTEGER, INTENT(IN) :: t_dump 
    INTEGER, INTENT(IN) :: fieldtype 
    REAL*4, INTENT(IN) :: fxk_re(0:nmax-1,0:nmax-1,0:nmax-1), fxk_im(0:nmax-1,0:nmax-1,0:nmax-1), &
                          fyk_re(0:nmax-1,0:nmax-1,0:nmax-1), fyk_im(0:nmax-1,0:nmax-1,0:nmax-1), &
                          fzk_re(0:nmax-1,0:nmax-1,0:nmax-1), fzk_im(0:nmax-1,0:nmax-1,0:nmax-1)
    INTEGER :: kx, ky, kz, kbins, ik
    REAL*4 :: kmin, kmax, dk, k, k0, fk_sqr
    REAL*4, ALLOCATABLE :: Pk(:) 
    CHARACTER(LEN=300) :: filename
    CHARACTER(LEN=6) :: uniti


    IF(t_dump<10) THEN
        WRITE(uniti,'(I1.1)') t_dump
    ELSE IF(t_dump>=10 .and. t_dump<100) THEN
        WRITE(uniti,'(I2.2)') t_dump
    ELSE IF(t_dump>=100 .and. t_dump<1000) THEN
        WRITE (uniti,'(I3.3)') t_dump
    ELSE IF(t_dump>=1000 .and. t_dump<10000) THEN
        WRITE (uniti,'(I4.3)') t_dump
    ELSE IF(t_dump>=10000 .and. t_dump<100000) THEN
        WRITE (uniti,'(I5.3)') t_dump  
    END IF   
    
    ! set number of bins 
    kbins = nmax
    
    ALLOCATE(Pk(1:kbins))
    Pk = 0.d0
    
    ! k band parameters
    k0 = TWOPI/Lx
    kmin = 0.d0
    kmax = SQRT(3.d0)*k0*(nmax/2)    
    dk = (kmax-kmin)/(kbins) ! uinform bin width
    
    ! loop over k space grid and deposit velocity fourier amplitudes (squared) into power spectrum bins (i.e. shells in k-space)
     DO kz = 0, nmax/2-1
        DO ky = 0, nmax/2-1
            DO kx = 0, nmax/2-1
                
                k = k0 * SQRT(SNGL(kx**2 + ky**2 + kz**2))
            
                fk_sqr = 2.d0 *(fxk_re(kx,ky,kz)**2 + fxk_im(kx,ky,kz)**2 + &
                         fyk_re(kx,ky,kz)**2 + fyk_im(kx,ky,kz)**2 + &
                         fzk_re(kx,ky,kz)**2 + fzk_im(kx,ky,kz)**2) 
                
                ! zero frequency component
                IF((kx .EQ. 0) .AND. (ky .EQ. 0) .AND. (kz .EQ. 0)) fk_sqr = 0.5d0 * fk_sqr
                
                ik = 1 + INT(k/dk)

                Pk(ik) = Pk(ik) + fk_sqr
    
            END DO
        END DO
    END DO
    
    
    Pk = Pk * (dx**3)    
    
    ! dump power spectrum into file
    IF(fieldtype .EQ. 0) THEN
        filename = TRIM('Output/density_Pk_dump=')//TRIM(uniti)//TRIM('.dat')
    ELSE IF(fieldtype .EQ. 1) THEN
        filename = TRIM('Output/velocity_Pk_dump=')//TRIM(uniti)//TRIM('.dat')
    ELSE IF(fieldtype .EQ. 2) THEN
        filename = TRIM('Output/Bfield_Pk_dump=')//TRIM(uniti)//TRIM('.dat')    
    ELSE IF(fieldtype .EQ. 3) THEN
        filename = TRIM('Output/v_wavelet_Pk_dump=')//TRIM(uniti)//TRIM('.dat')    
    END IF    


    OPEN(UNIT=1, FILE=filename, FORM = 'UNFORMATTED', ACCESS = 'STREAM')
    
    DO ik = 1, kbins
        WRITE(1) (ik-0.5)*dk,Pk(ik)
    END DO
    
    CLOSE(UNIT=1)
    
    DEALLOCATE(Pk)


END SUBROUTINE compute_pk



SUBROUTINE compute_density_pdf(t_dump)

    INTEGER, INTENT(IN) :: t_dump 
    INTEGER :: i, j, k, nbins, ibin
    REAL*4 :: rhomin, rhomax, drho, pdfnorm
    REAL*4, ALLOCATABLE :: PDF(:) 
    CHARACTER(LEN=300) :: filename
    CHARACTER(LEN=6) :: uniti


    IF(t_dump<10) THEN
        WRITE(uniti,'(I1.1)') t_dump
    ELSE IF(t_dump>=10 .and. t_dump<100) THEN
        WRITE(uniti,'(I2.2)') t_dump
    ELSE IF(t_dump>=100 .and. t_dump<1000) THEN
        WRITE (uniti,'(I3.3)') t_dump
    ELSE IF(t_dump>=1000 .and. t_dump<10000) THEN
        WRITE (uniti,'(I4.3)') t_dump
    ELSE IF(t_dump>=10000 .and. t_dump<100000) THEN
        WRITE (uniti,'(I5.3)') t_dump  
    END IF   
    
    ! set number of bins 
    nbins = nmax
    
    ALLOCATE(PDF(1:nbins))
    PDF = 0.d0
    
    ! set upper/lower bounds and bin size
    rhomin = 0.d0
    rhomax = 4.d0    
    drho = (rhomax-rhomin)/(nbins) ! uinform bin width 
    
    
    DO k = 1, nranks_z*nz
        DO j = 1, nranks_y*ny
            DO i = 1, nranks_x*nx

                ibin = 1 + INT(fx(i,j,k,1)/drho)
                PDF(ibin) = PDF(ibin) + 1.d0
                
            END DO
        END DO
    END DO

    ! normalization
    pdfnorm = drho*SUM(PDF)
    PDF = PDF / pdfnorm

    ! dump density pdf into file
    filename = TRIM('Output/density_PDF_dump=')//TRIM(uniti)//TRIM('.dat')
    OPEN(UNIT=1, FILE=filename, FORM = 'UNFORMATTED', ACCESS = 'STREAM')
    
    DO ibin = 1, nbins
        WRITE(1) (ibin-0.5)*drho,PDF(ibin)
    END DO
    
    CLOSE(UNIT=1)
    
    DEALLOCATE(PDF)

END SUBROUTINE compute_density_pdf


!SUBROUTINE enstrophy()
!
!END SUBROUTINE enstrophy



SUBROUTINE shift_fft(fxk_re, fxk_im)

    REAL*4, INTENT(INOUT) :: fxk_re(-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2), &
                             fxk_im(-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2)  
    INTEGER :: kx, ky, kz, kx_shft, ky_shft, kz_shft
    REAL*4, ALLOCATABLE :: fs_re(:,:,:), fs_im(:,:,:)
    
    ALLOCATE(fs_re(-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2))
    ALLOCATE(fs_im(-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2))
    
    DO kz = -nmax/2, -1+nmax/2
        DO ky = -nmax/2, -1+nmax/2
            DO kx = -nmax/2, -1+nmax/2

                kx_shft = kx - nmax * (-1 + (SIGN(1,kx)))/2 
                ky_shft = ky - nmax * (-1 + (SIGN(1,ky)))/2 
                kz_shft = kz - nmax * (-1 + (SIGN(1,kz)))/2          

                fs_re(kx,ky,kz) = fxk_re(kx_shft,ky_shft,kz_shft)
                fs_im(kx,ky,kz) = fxk_im(kx_shft,ky_shft,kz_shft)
    
            END DO
        END DO
    END DO
    
    DO kz = -nmax/2, -1+nmax/2
        DO ky = -nmax/2, -1+nmax/2
            DO kx = -nmax/2, -1+nmax/2

                fxk_re(kx+nmax/2,ky+nmax/2,kz+nmax/2) = fs_re(kx,ky,kz)
                fxk_im(kx+nmax/2,ky+nmax/2,kz+nmax/2) = fs_im(kx,ky,kz)
    
            END DO
        END DO
    END DO
    
    

    DEALLOCATE(fs_re, fs_im)

END SUBROUTINE shift_fft



SUBROUTINE compute_ke()

    INTEGER :: i, j, k
    REAL(4) :: rho 

    total_ke = 0.d0

    DO k = 1, nranks_z*nz
    DO j = 1, nranks_y*ny
    DO i = 1, nranks_x*nx

        rho = fx(i,j,k,1) 
        total_ke = total_ke + (fx(i,j,k,2)**2 + fx(i,j,k,3)**2 + fx(i,j,k,4)**2) /rho  

    END DO
    END DO
    END DO

    total_ke = 0.5d0 * total_ke * (dx**3)

    PRINT*,'Kinetic Energy = ',total_ke 


END SUBROUTINE compute_ke


SUBROUTINE compute_vrms()

    INTEGER :: i,j,k
    REAL(4) :: rho 

    v_rms = 0.d0

    DO k = 1, nranks_z*nz
    DO j = 1, nranks_y*ny
    DO i = 1, nranks_x*nx

        rho = fx(i,j,k,1) 
        v_rms = v_rms + (fx(i,j,k,2)**2 + fx(i,j,k,3)**2 + fx(i,j,k,4)**2) /(rho**2)  

    END DO
    END DO
    END DO

    v_rms = SQRT(v_rms*(dx**3))  

    PRINT*,'rms velocity = ',v_rms 


END SUBROUTINE compute_vrms


SUBROUTINE compute_brms()

    INTEGER :: i,j,k

    b_rms = 0.d0

    DO k = 1, nranks_z*nz
    DO j = 1, nranks_y*ny
    DO i = 1, nranks_x*nx

        b_rms = b_rms + (fx(i,j,k,5)**2 + fx(i,j,k,6)**2 + fx(i,j,k,7)**2)   

    END DO
    END DO
    END DO

    b_rms = SQRT(b_rms*(dx**3)) 

    PRINT*,'rms magnetic field = ',b_rms 


END SUBROUTINE compute_brms



SUBROUTINE save_fft_to_file(t_dump)

    INTEGER, INTENT(IN) :: t_dump 
    INTEGER :: kx, ky, kz, kx_shft, ky_shft, kz_shft
    CHARACTER(LEN=300) :: filename
    CHARACTER(LEN=6) :: uniti


    PRINT*,'Saving FFT to file..'

    IF(t_dump<10) THEN
        WRITE(uniti,'(I1.1)') t_dump
    ELSE IF(t_dump>=10 .and. t_dump<100) THEN
        WRITE(uniti,'(I2.2)') t_dump
    ELSE IF(t_dump>=100 .and. t_dump<1000) THEN
        WRITE (uniti,'(I3.3)') t_dump
    ELSE IF(t_dump>=1000 .and. t_dump<10000) THEN
        WRITE (uniti,'(I4.3)') t_dump
    ELSE IF(t_dump>=10000 .and. t_dump<100000) THEN
        WRITE (uniti,'(I5.3)') t_dump  
    END IF
    


    filename = TRIM('Output/velocity_vxk_dump=')//TRIM(uniti)//TRIM('.dat')

    OPEN(UNIT=1, FILE=filename, FORM = 'UNFORMATTED', ACCESS = 'STREAM')
    
    DO kz = -nmax/2, -1+nmax/2
        DO ky = -nmax/2, -1+nmax/2
            DO kx = -nmax/2, -1+nmax/2

                kx_shft = kx - nmax * (-1 + (SIGN(1,kx)))/2 
                ky_shft = ky - nmax * (-1 + (SIGN(1,ky)))/2 
                kz_shft = kz - nmax * (-1 + (SIGN(1,kz)))/2          

                !WRITE(1) fxk_re(kx_shft,ky_shft,kz_shft), fxk_im(kx_shft,ky_shft,kz_shft)
                !WRITE(1) fxk_re(kx+nmax/2,ky+nmax/2,kz+nmax/2), fxk_im(kx+nmax/2,ky+nmax/2,kz+nmax/2)
    
            END DO
        END DO
    END DO

    
    CLOSE(UNIT=1)

END SUBROUTINE save_fft_to_file


END PROGRAM analysis