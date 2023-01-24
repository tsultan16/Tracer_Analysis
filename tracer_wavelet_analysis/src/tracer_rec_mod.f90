MODULE tracer_reconstruction_mod

!USE MPI
USE OMP_LIB
USE constants_mod
USE grid_arrays_mod
USE fwt_mod

IMPLICIT NONE

CONTAINS


SUBROUTINE tracer_reconstruct()

    INTEGER :: ii, i, j, k, nranks_tot
    REAL :: time, dt, cmplt
    


    ! set wavelet transform grid-size
    CALL fwt_init(nwx,nwy,nwz)
    
    cmplt = 0.D0    
    
    nranks_tot = nranks_x * nranks_y * nranks_z
    
    ! loop over MPI ranks
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i,j,k)
    DO ii = 1, nranks_tot
    
        CALL collapse_index(ii,i,j,k)

        PRINT*,''
        PRINT*,'Reconstructing for rank: ',i,j,k
  
        CALL tracer_wavelet_reconstruction(i, j, k)     
   
        !$OMP ATOMIC
        cmplt = cmplt + 1
   
        PRINT*,'Done reconstructing for rank: ',i,j,k
        PRINT*,'Percent complete = ',100.0*REAL(cmplt)/REAL(nranks_tot)
                
    END DO
   !$OMP END PARALLEL DO

    CALL fwt_destroy()

    ! save to file
    IF(tracer_rec_file) CALL file_output()
   

END SUBROUTINE tracer_reconstruct


! partial reconstruction using wavelets picked up by the tracers
SUBROUTINE tracer_wavelet_reconstruction(ri, rj, rk)

    INTEGER, INTENT(IN) :: ri, rj, rk
    REAL(4) :: fbasis(nwx,nwy,nwz), amp
    INTEGER :: i, j, k, ii, ix, iy, iz, isc, nn(nscalars), ilow, ihi, jlow, jhi, klow, khi
   
    nn = 0


        ! loop over particles
        DO i = 1, dat_rec(ri,rj,rk)%nparticles
        
            !PRINT*,'Particle# ',i
            
            ! loop over wavelets
            DO j = 1, nwavelets_per_scalar
                
                IF(INT(dat_rec(ri,rj,rk)%tracer_wavelets(i,j,1)) .LE. 0) CYCLE
        
                ! flattened wavelet index
                ii = dat_rec(ri,rj,rk)%tracer_wavelets(i,j,1)
                
                ! collapse flattened wavelet index
                ix = 1 + MOD(ii-1,nwx)
                iy = 1 + MOD(ii-1,nwx*nwy)/nwx 
                iz = 1 + (ii-1) / (nwx*nwy) 
                                                                                    
                IF(ix .LT. 1 .OR. ix .GT. nwx .OR. iy .LT. 1 .OR. iy .GT. nwy .OR. iz .LT. 1 .OR. iz .GT. nwz) THEN
                    PRINT*,'ERROR!! Invalid wavelet index: ',ix,iy,iz
                    STOP
                END IF       
                                               
                ! compute basis wavelet
                CALL compute_basis_wavelet(ix,iy,iz,fbasis)
       
                ! add wavelet contribution to the reconstructed function
                DO isc = 1, nscalars
                    IF(keep_bndry) THEN
                        dat_rec(ri,rj,rk)%patch_array(:,:,:,isc) = dat_rec(ri,rj,rk)%patch_array(:,:,:,isc) + dat_rec(ri,rj,rk)%tracer_wavelets(i,j,1+isc)*fbasis(:,:,:) 
                    ELSE
                        dat_rec(ri,rj,rk)%patch_array(1:nx,1:ny,1:nz,isc) = dat_rec(ri,rj,rk)%patch_array(1:nx,1:ny,1:nz,isc) + dat_rec(ri,rj,rk)%tracer_wavelets(i,j,1+isc)*fbasis(:,:,:) 
                    END IF
                    
                    nn(isc) = nn(isc) + 1 
                END DO
            
            END DO
        
        END DO
     
 
    ! append reconstructed field to full-grid array 
    ilow = file_meta(ri,rj,rk)%origin_zone(1) 
    ihi = ilow + nx - 1
    jlow = file_meta(ri,rj,rk)%origin_zone(2) 
    jhi = jlow + ny - 1 
    klow = file_meta(ri,rj,rk)%origin_zone(3) 
    khi = klow + nz - 1    
    
    
    fx(ilow:ihi,jlow:jhi,klow:khi,wombat_dump_nvars+1) = dat_rec(ri,rj,rk)%patch_array(1:nx,1:ny,1:nz,1)
    fx(ilow:ihi,jlow:jhi,klow:khi,wombat_dump_nvars+2) = dat_rec(ri,rj,rk)%patch_array(1:nx,1:ny,1:nz,2)
    fx(ilow:ihi,jlow:jhi,klow:khi,wombat_dump_nvars+3) = dat_rec(ri,rj,rk)%patch_array(1:nx,1:ny,1:nz,3)
    
 
    PRINT*,'Total number of wavelets used in the reconstruction (per scalars) =',nn(1)
 
 
END SUBROUTINE tracer_wavelet_reconstruction


SUBROUTINE collapse_index(ii, i, j, k)

    INTEGER, INTENT(IN) :: ii
    INTEGER, INTENT(OUT) :: i, j, k

    
    i = MOD(ii-1, nranks_x)
    j = MOD(ii-1, nranks_x*nranks_y) / nranks_x 
    k = (ii-1) / (nranks_x*nranks_y) 
              
END SUBROUTINE collapse_index


SUBROUTINE file_output()

    INTEGER :: i, j, k, isc, ix, iy, iz
    CHARACTER(LEN=100) :: filename
    CHARACTER(LEN=6) :: uniti
   
    IF(dump_num<10) THEN
        WRITE(uniti,'(I1.1)') dump_num
    ELSE IF(dump_num>=10 .and. dump_num<100) THEN
        WRITE(uniti,'(I2.2)') dump_num
    ELSE IF(dump_num>=100 .and. dump_num<1000) THEN
        WRITE (uniti,'(I3.3)') dump_num
    ELSE IF(dump_num>=1000 .and. dump_num<10000) THEN
        WRITE (uniti,'(I4.3)') dump_num
    ELSE IF(dump_num>=10000 .and. dump_num<100000) THEN
        WRITE (uniti,'(I5.3)') dump_num  
    END IF
    
    PRINT*,''
    PRINT*,'Writing to file...'
    
    
    filename =TRIM('Output/reconstruction')//TRIM('_dump=')//TRIM(uniti)//TRIM('.dat')                           
    OPEN(UNIT = 10, FILE = filename, FORM = 'UNFORMATTED', ACCESS='STREAM')

    DO i = 0, nranks_x-1
        DO j = 0, nranks_y-1
            DO k = 0, nranks_z-1
    
     
                DO isc = 1, nscalars
                DO iz = 1, nz
                DO iy = 1, ny
                DO ix = 1, nx
    
                WRITE(10) dat_rec(i,j,k)%patch_array(ix,iy,iz,isc)
        
                END DO
                END DO
                END DO
                END DO
        
            END DO
        END DO
    END DO
   
    CLOSE(UNIT=10) 

   
END SUBROUTINE file_output


END MODULE tracer_reconstruction_mod