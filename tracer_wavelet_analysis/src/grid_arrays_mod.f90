MODULE grid_arrays_mod

USE constants_mod
IMPLICIT NONE


INTEGER :: nxtot, nytot, nztot, nwx, nwy, nwz
INTEGER ::  nwavelets_per_scalar, nwavelets_max

LOGICAL :: particle_presence(0:nranks_x*nranks_y*nranks_z -1 ,1:192) = .FALSE.

TYPE :: patch_scalar

    INTEGER :: nparticles, nwavelet_indices
    REAL(4), ALLOCATABLE :: patch_array(:,:,:,:)
    REAL(4), ALLOCATABLE :: tracer_wavelets(:,:,:)
    REAL(4), ALLOCATABLE :: wavelet_avg_rhoB(:,:)
    
    REAL(4), ALLOCATABLE :: wavelet_vxk_re(:,:,:)
    REAL(4), ALLOCATABLE :: wavelet_vxk_im(:,:,:)
    REAL(4), ALLOCATABLE :: wavelet_vyk_re(:,:,:)
    REAL(4), ALLOCATABLE :: wavelet_vyk_im(:,:,:)
    REAL(4), ALLOCATABLE :: wavelet_vzk_re(:,:,:)
    REAL(4), ALLOCATABLE :: wavelet_vzk_im(:,:,:)

    REAL(4), ALLOCATABLE :: wavelet_va_k_re(:,:,:)
    REAL(4), ALLOCATABLE :: wavelet_va_k_im(:,:,:)
    REAL(4), ALLOCATABLE :: wavelet_vf_k_re(:,:,:)
    REAL(4), ALLOCATABLE :: wavelet_vf_k_im(:,:,:)
    REAL(4), ALLOCATABLE :: wavelet_vs_k_re(:,:,:)
    REAL(4), ALLOCATABLE :: wavelet_vs_k_im(:,:,:)
    
    REAL(4), ALLOCATABLE :: wavelet_Pk_v(:,:)
    REAL(4), ALLOCATABLE :: wavelet_Pk_b(:,:)
    REAL(4), ALLOCATABLE :: wavelet_Pk_rho(:,:)

    INTEGER, ALLOCATABLE :: wavelet_lookup_table(:)


END TYPE patch_scalar


TYPE :: patch_metadata

    REAL(4) :: dt, t
    INTEGER :: worldsize(3), patchsize(3), ndomains(3), mpi_rank, rank_coordinates(3), origin_zone(3), nparticles, nwavelet_indices
   
END TYPE patch_metadata


TYPE(patch_scalar), ALLOCATABLE :: dat_rec(:,:,:)
TYPE(patch_metadata), ALLOCATABLE :: file_meta(:,:,:)


REAL*4, ALLOCATABLE :: fx(:,:,:,:), buffer(:,:,:), vp(:,:,:,:)
REAL*4, ALLOCATABLE :: vxk_re(:,:,:), vxk_im(:,:,:),vyk_re(:,:,:), vyk_im(:,:,:),vzk_re(:,:,:), vzk_im(:,:,:), &
                       bxk_re(:,:,:), bxk_im(:,:,:),byk_re(:,:,:), byk_im(:,:,:),bzk_re(:,:,:), bzk_im(:,:,:), &
                       rhok_re(:,:,:), rhok_im(:,:,:), va_k_re(:,:,:), va_k_im(:,:,:), vf_k_re(:,:,:), vf_k_im(:,:,:), &
                       vs_k_re(:,:,:), vs_k_im(:,:,:), &
                       vwx(:,:,:), vwy(:,:,:), vwz(:,:,:), fwbk_re(:,:,:), fwbk_im(:,:,:)
                       



INTEGER :: nmax, mem_bytes


CONTAINS


SUBROUTINE create_patch_arrays()


    INTEGER :: i, j, k, nparticles  
    LOGICAL :: found
    
    found = .FALSE.

    ! world grid size
    nxtot = nranks_x*nx
    nytot = nranks_y*ny
    nztot = nranks_z*nz

    ! wavelet enlarged patch size
    IF(keep_bndry) THEN
        nwx = nx  + 2*nb 
        nwy = ny  + 2*nb
        nwz = nz  + 2*nb
    ELSE
        nwx = nx   
        nwy = ny  
        nwz = nz  
    END IF

    nmax = MIN(nxtot, nytot, nztot)

    ! max number of wavelets per tracer
    nwavelets_per_scalar =  nwavelets_per_bin * (1+2*bin_neighborhood)**3
    nwavelets_max = nwavelets_per_scalar

    PRINT*,''
    PRINT*,'---> Constructing patch arrays and reading tracer data for dump # ',dump_num
    PRINT*,''

    ALLOCATE(dat_rec(0:nranks_x-1,0:nranks_y-1,0:nranks_z-1))
    ALLOCATE(file_meta(0:nranks_x-1,0:nranks_y-1,0:nranks_z-1))
    
  
    
    ! obtain file meta data for all patches and allocate tracer data arrays
    CALL read_metadata()             
        
    DO k = 0, nranks_z-1
        DO j = 0, nranks_y-1
            DO i = 0, nranks_x-1
    
                dat_rec(i,j,k)%nparticles = file_meta(i,j,k)%nparticles
                dat_rec(i,j,k)%nwavelet_indices = file_meta(i,j,k)%nwavelet_indices

                ALLOCATE(dat_rec(i,j,k)%patch_array(1:nx+nb,1-nb:ny+nb,1-nb:nz+nb,1:nscalars))
                dat_rec(i,j,k)%patch_array = 0.D0
    
                ! allocate memory for tracer data
                ALLOCATE(dat_rec(i,j,k)%tracer_wavelets(1:dat_rec(i,j,k)%nparticles,nwavelets_max+1,4))  ! added an extra slot on the second index, will use it for storing div(V) value at tracer location
                dat_rec(i,j,k)%tracer_wavelets = 0.0
    
                ! allocate memory for wavelet averages 
                ALLOCATE(dat_rec(i,j,k)%wavelet_avg_rhoB(1:dat_rec(i,j,k)%nwavelet_indices,5))
                dat_rec(i,j,k)%wavelet_avg_rhoB = 0.0
    
                ! allocate memory for wavelet FFT arrays
                ALLOCATE(dat_rec(i,j,k)%wavelet_vxk_re(0:nx-1,0:ny-1,0:nz-1))
                ALLOCATE(dat_rec(i,j,k)%wavelet_vxk_im(0:nx-1,0:ny-1,0:nz-1))
                ALLOCATE(dat_rec(i,j,k)%wavelet_vyk_re(0:nx-1,0:ny-1,0:nz-1))
                ALLOCATE(dat_rec(i,j,k)%wavelet_vyk_im(0:nx-1,0:ny-1,0:nz-1))
                ALLOCATE(dat_rec(i,j,k)%wavelet_vzk_re(0:nx-1,0:ny-1,0:nz-1))
                ALLOCATE(dat_rec(i,j,k)%wavelet_vzk_im(0:nx-1,0:ny-1,0:nz-1))
                ALLOCATE(dat_rec(i,j,k)%wavelet_va_k_re(0:nx-1,0:ny-1,0:nz-1))
                ALLOCATE(dat_rec(i,j,k)%wavelet_va_k_im(0:nx-1,0:ny-1,0:nz-1))
                ALLOCATE(dat_rec(i,j,k)%wavelet_vf_k_re(0:nx-1,0:ny-1,0:nz-1))
                ALLOCATE(dat_rec(i,j,k)%wavelet_vf_k_im(0:nx-1,0:ny-1,0:nz-1))
                ALLOCATE(dat_rec(i,j,k)%wavelet_vs_k_re(0:nx-1,0:ny-1,0:nz-1))
                ALLOCATE(dat_rec(i,j,k)%wavelet_vs_k_im(0:nx-1,0:ny-1,0:nz-1))
                ALLOCATE(dat_rec(i,j,k)%wavelet_Pk_v(1:nx,3))
                ALLOCATE(dat_rec(i,j,k)%wavelet_Pk_b(1:nx,3))
                ALLOCATE(dat_rec(i,j,k)%wavelet_Pk_rho(1:nx,2))
                ALLOCATE(dat_rec(i,j,k)%wavelet_lookup_table(1:nx*ny*nz))
                dat_rec(i,j,k)%wavelet_vxk_re = 0.d0
                dat_rec(i,j,k)%wavelet_vxk_im = 0.d0
                dat_rec(i,j,k)%wavelet_vyk_re = 0.d0
                dat_rec(i,j,k)%wavelet_vyk_im = 0.d0
                dat_rec(i,j,k)%wavelet_vzk_re = 0.d0
                dat_rec(i,j,k)%wavelet_vzk_im = 0.d0
                dat_rec(i,j,k)%wavelet_va_k_re = 0.d0
                dat_rec(i,j,k)%wavelet_va_k_im = 0.d0
                dat_rec(i,j,k)%wavelet_vf_k_re = 0.d0
                dat_rec(i,j,k)%wavelet_vf_k_im = 0.d0
                dat_rec(i,j,k)%wavelet_vs_k_re = 0.d0
                dat_rec(i,j,k)%wavelet_vs_k_im = 0.d0
                
                dat_rec(i,j,k)%wavelet_Pk_v = 0.d0
                dat_rec(i,j,k)%wavelet_Pk_b = 0.d0
                dat_rec(i,j,k)%wavelet_Pk_rho = 0.d0
                dat_rec(i,j,k)%wavelet_lookup_table = 0
    
                ! now read tracer data      
                WRITE(*,FMT='("Reading tracer data for rank=",i3,", coords = ",i2,i2,i2, ", origin_zone = ", i3,i3,i3)') file_meta(i,j,k)%mpi_rank,file_meta(i,j,k)%rank_coordinates(1), &
                      file_meta(i,j,k)%rank_coordinates(2),file_meta(i,j,k)%rank_coordinates(3),file_meta(i,j,k)%origin_zone(1),file_meta(i,j,k)%origin_zone(2),file_meta(i,j,k)%origin_zone(3)     
                CALL read_tracers(i,j,k,file_meta(i,j,k)%mpi_rank, found)
    
            END DO
        END DO
    END DO
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    PRINT*,''
    PRINT*,'ALL(particle_presence) = ',ALL(particle_presence)
    PRINT*,''
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ALLOCATE(vxk_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(vxk_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(vyk_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(vyk_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(vzk_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(vzk_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(bxk_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(bxk_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(byk_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(byk_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(bzk_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(bzk_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(rhok_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(rhok_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(va_k_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(va_k_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(vf_k_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(vf_k_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(vs_k_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(vs_k_im(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(fwbk_re(0:nmax-1,0:nmax-1,0:nmax-1))
    ALLOCATE(fwbk_im(0:nmax-1,0:nmax-1,0:nmax-1))

    ALLOCATE(fx(1:nxtot,1:nytot,1:nztot,1:wombat_dump_nvars+3))
    ALLOCATE(vp(1:nwx,1:nwy,1:nwz,3))
    ALLOCATE(buffer(1:nxtot,1:nytot,1:nztot))
    ALLOCATE(vwx(1:nwx,1:nwy,1:nwz))
    ALLOCATE(vwy(1:nwx,1:nwy,1:nwz))
    ALLOCATE(vwz(1:nwx,1:nwy,1:nwz))

    fx = 0.D0

    
    mem_bytes = 22*SIZEOF(vxk_re) + SIZEOF(fx) + SIZEOF(vp) + 4*SIZEOF(buffer) + SIZEOF(dat_rec)
   
    PRINT*,''
    PRINT*,'Memory allocated for work arrays (Mb) = ', mem_bytes*1.e-6
    PRINT*,''


END SUBROUTINE create_patch_arrays


SUBROUTINE destroy_patch_arrays()

    INTEGER :: i, j, k
    
    DO k = 0, nranks_z-1
        DO j = 0, nranks_y-1
            DO i = 0, nranks_x-1
    
                DEALLOCATE(dat_rec(i,j,k)%patch_array)
                DEALLOCATE(dat_rec(i,j,k)%tracer_wavelets)
                DEALLOCATE(dat_rec(i,j,k)%wavelet_vxk_re)
                DEALLOCATE(dat_rec(i,j,k)%wavelet_vxk_im)
                DEALLOCATE(dat_rec(i,j,k)%wavelet_vyk_re)
                DEALLOCATE(dat_rec(i,j,k)%wavelet_vyk_im)
                DEALLOCATE(dat_rec(i,j,k)%wavelet_vzk_re)
                DEALLOCATE(dat_rec(i,j,k)%wavelet_vzk_im)
                DEALLOCATE(dat_rec(i,j,k)%wavelet_va_k_re)
                DEALLOCATE(dat_rec(i,j,k)%wavelet_va_k_im)
                DEALLOCATE(dat_rec(i,j,k)%wavelet_vf_k_re)
                DEALLOCATE(dat_rec(i,j,k)%wavelet_vf_k_im)
                DEALLOCATE(dat_rec(i,j,k)%wavelet_vs_k_re)
                DEALLOCATE(dat_rec(i,j,k)%wavelet_vs_k_im)
                DEALLOCATE(dat_rec(i,j,k)%wavelet_Pk_v)
                DEALLOCATE(dat_rec(i,j,k)%wavelet_Pk_b)
                DEALLOCATE(dat_rec(i,j,k)%wavelet_Pk_rho)
                DEALLOCATE(dat_rec(i,j,k)%wavelet_lookup_table)
                
            END DO
        END DO
    END DO
    
    DEALLOCATE(dat_rec)
    DEALLOCATE(file_meta)
    
    DEALLOCATE(fx, vp, buffer)
    DEALLOCATE(vxk_re, vxk_im, vyk_re, vyk_im, vzk_re, vzk_im)
    DEALLOCATE(bxk_re, bxk_im, byk_re, byk_im, bzk_re, bzk_im)
    DEALLOCATE(rhok_re, rhok_im)
    DEALLOCATE(va_k_re, va_k_im, vf_k_re, vf_k_im, vs_k_re, vs_k_im)
    DEALLOCATE(vwx, vwy, vwz, fwbk_re, fwbk_im)
    

END SUBROUTINE destroy_patch_arrays



SUBROUTINE read_metadata()

    CHARACTER(LEN=200) :: filename
    CHARACTER(LEN=6) :: uniti
    INTEGER :: i, j, k, byte_offset, item_bytes, ntot
    REAL(4) :: meta_buff(20)
    
    ntot = 0
    
    ! byte size per data item
    item_bytes = 4
   
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

    filename = TRIM(output_filepath)//TRIM('TRACER/tracer_parallel_meta_dump=')//TRIM(uniti)//TRIM('.dat')        
           
    OPEN(UNIT=10, FILE = filename, FORM = 'UNFORMATTED', STATUS = 'OLD', ACCESS = 'STREAM')
   
    byte_offset = 1
  
    DO k = 0, nranks_z-1
        DO j = 0, nranks_y-1
            DO i = 0, nranks_x-1

                READ(10, POS = byte_offset) meta_buff
                byte_offset = byte_offset + 20*item_bytes 
    
                file_meta(i,j,k)%t = meta_buff(1)
                file_meta(i,j,k)%dt = meta_buff(2)
                file_meta(i,j,k)%worldsize(1:3) = meta_buff(3:5)
                file_meta(i,j,k)%patchsize(1:3) = meta_buff(6:8)
                file_meta(i,j,k)%ndomains(1:3) = meta_buff(9:11)
                file_meta(i,j,k)%mpi_rank = meta_buff(12)
                file_meta(i,j,k)%rank_coordinates(1:3) = meta_buff(13:15)
                file_meta(i,j,k)%origin_zone(1:3) = meta_buff(16:18)
                file_meta(i,j,k)%nparticles = meta_buff(19)
                file_meta(i,j,k)%nwavelet_indices = meta_buff(20)

                ntot = ntot + file_meta(i,j,k)%nparticles
    
                PRINT*,''
                PRINT*,'mpi_rank = ',file_meta(i,j,k)%mpi_rank
                PRINT*,'Rank Coords = ',file_meta(i,j,k)%rank_coordinates
                PRINT*,'origin zone = ',file_meta(i,j,k)%origin_zone
                PRINT*,'nparticles = ',file_meta(i,j,k)%nparticles
                PRINT*,'nwavelet_indices = ',file_meta(i,j,k)%nwavelet_indices
                PRINT*,''               

            END DO
        END DO
    END DO

    PRINT*,'#############################'
    PRINT*,'Nparticles_total = ',ntot
    PRINT*,'#############################'
    PRINT*,''
    
    CLOSE(UNIT=10)
   
END SUBROUTINE read_metadata



SUBROUTINE read_tracers(ri, rj, rk, myrank, found)

    INTEGER, INTENT(IN) :: ri, rj, rk, myrank
    CHARACTER(LEN=200) :: filename
    CHARACTER(LEN=6) :: uniti
    INTEGER :: i, j, k, isc, item_bytes, byte_offset, ix, iy, iz
    REAL(4) :: pid, init_rank, xflat, divV, tmp
    LOGICAL, INTENT(INOUT) :: found
    LOGICAL :: fexist
    
   
    ! byte size per tracer data item
    item_bytes = 4
        
    IF(dump_num<10) THEN
        WRITE(uniti,'(I1.1)') dump_num
    ELSE IF(dump_num>=10 .and. dump_num < 100) THEN
        WRITE(uniti,'(I2.2)') dump_num
    ELSE IF(dump_num>=100 .and. dump_num < 1000) THEN
        WRITE (uniti,'(I3.3)') dump_num
    ELSE IF(dump_num>=1000 .and. dump_num < 10000) THEN
        WRITE (uniti,'(I4.3)') dump_num
    ELSE IF(dump_num>=10000 .and. dump_num < 100000) THEN
        WRITE (uniti,'(I5.3)') dump_num  
    END IF

    filename = TRIM(output_filepath)//TRIM('TRACER/tracer_parallel_dump=')//TRIM(uniti)//TRIM('.dat')        
   
    OPEN(UNIT=10, FILE = filename, FORM = 'UNFORMATTED', STATUS = 'OLD', ACCESS = 'STREAM')
   
    ! get the byte offset for this mpi rank
    byte_offset = 1 + myrank * item_bytes
    
    READ(10, POS = byte_offset) tmp  
    byte_offset = 1 + INT(tmp) 
   
    ! read tracer particle data
    DO i = 1, dat_rec(ri,rj,rk)%nparticles  
    
        
        READ(10, POS = byte_offset) pid 
        READ(10, POS = byte_offset + item_bytes) init_rank 
        READ(10, POS = byte_offset + 2*item_bytes) xflat
        READ(10, POS = byte_offset + 3*item_bytes) divV
        byte_offset  = byte_offset + 4*item_bytes

        ! collapse flattened particle position index
        ix = MOD(INT(xflat)-1,nx+2)
        iy = MOD(INT(xflat)-1,(nx+2)*(ny+2)) / (nx+2) 
        iz = (INT(xflat)-1) / ((nx+2)*(ny+2))


        !PRINT*,'PARTICLE# ',i
        !PRINT*,'Particle id, init_rank = ',INT(pid), INT(init_rank)
        !PRINT*,'Particle position (flat, i,j,k) = ',INT(xflat),ix,iy,iz
        !PRINT*,'Div(V) = ',divV            

        IF(ix .LT. 0 .OR. ix .GT. nx+1 .OR. iy .LT. 0 .OR. iy .GT. ny+1 .OR. iz .LT. 0 .OR. iz .GT. nz+1) THEN
            PRINT*,'ERROR!! Invalid particle coordinates!! '
            PRINT*,'Particle id, init_rank = ',pid, init_rank
            PRINT*,'Particle position (flat, i,j,k) = ',INT(xflat),ix,iy,iz
            STOP
        END IF
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(.NOT. particle_presence(INT(init_rank),INT(pid))) THEN
            particle_presence(INT(init_rank),INT(pid)) = .TRUE.
        ELSE
            PRINT*,'ERROR! Found duplicate particle, id, init_rank = ',INT(pid),INT(init_rank)
            STOP
        END IF    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
     
        DO j = 1, nwavelets_per_scalar
                READ(10, POS = byte_offset) dat_rec(ri,rj,rk)%tracer_wavelets(i,j,1)                           ! wavelet index
                READ(10, POS = byte_offset + item_bytes) dat_rec(ri,rj,rk)%tracer_wavelets(i,j,2:2+nscalars-1) ! wavelet amplitudes for x,y,z components
                byte_offset = byte_offset + (1+nscalars) * item_bytes            
                !PRINT*,'Wavelet #, Amplitude = ',dat_rec(ri,rj,rk)%tracer_wavelets(i,j,:)            
        END DO
        
        dat_rec(ri,rj,rk)%tracer_wavelets(i,nwavelets_per_scalar+1,1) = divV  ! velocity divergence at tracer location
        
        
        IF(INT(pid) .EQ. 10 .AND. INT(init_rank) .EQ. 1)THEN
           
            PRINT*,''
            PRINT*,'%%%%%%%%%%%%%%%'
            PRINT*,'Particle found.'
            PRINT*,'%%%%%%%%%%%%%%%'
            PRINT*,''

            PRINT*,'Particle id, init_rank = ',INT(pid), INT(init_rank)
            PRINT*,'Particle position (flat, i,j,k) = ',INT(xflat),ix,iy,iz
            PRINT*,'Div(V) = ',divV       
        
            IF(.NOT. found) THEN
                ! append tracer data to a unique file
                !INQUIRE(FILE="tracer_100_0.dat", EXIST=fexist)
                
                !PRINT*,'fexist = ',fexist
                
                !IF (fexist) THEN
                    OPEN(UNIT=12, FILE="Output/tracer_6_1.dat", FORM = 'UNFORMATTED', STATUS = 'UNKNOWN', ACCESS = 'STREAM', POSITION='APPEND', ACTION='WRITE')
                    WRITE(12) file_meta(ri,rj,rk)%t,file_meta(ri,rj,rk)%dt,pid,init_rank,ri,rj,rk,ix,iy,iz
                !ELSE
                !    OPEN(UNIT=12, FILE="Output/tracer_1_0.dat", FORM = 'UNFORMATTED', STATUS = 'NEW', ACCESS = 'STREAM', ACTION='WRITE')
                !    WRITE(12) file_meta(ri,rj,rk)%t,file_meta(ri,rj,rk)%dt,pid,init_rank,ri,rj,rk,ix,iy,iz
                !END IF

                CLOSE(12)   
                found = .TRUE.
                
            ELSE
                PRINT*,'ERROR!!! Multiple copies found for the same particle!!!'
                !STOP
            END IF        
                    
        END IF
        
        
    END DO
    
    ! read wavelet average data
    DO i = 1, dat_rec(ri,rj,rk)%nwavelet_indices
        READ(10, POS = byte_offset) dat_rec(ri,rj,rk)%wavelet_avg_rhoB(i,1:5) 
        byte_offset = byte_offset + 5 * item_bytes 
        
        !IF(ri .EQ. 0 .AND. rj .EQ. 2 .AND. rk .EQ. 1 .AND. i .GE. 6609) PRINT*,'Wavelet Index, Avg = ',dat_rec(ri,rj,rk)%wavelet_avg_rhoB(i,1:5)
        
        ! record array index in lookup table
        dat_rec(ri,rj,rk)%wavelet_lookup_table(INT(dat_rec(ri,rj,rk)%wavelet_avg_rhoB(i,1))) = i
    END DO    
    
    CLOSE(10)
   
    !IF(found)THEN
    !    PRINT*,'Particle found.'
    !ELSE
    !    PRINT*,'Particle not found.'
    !END IF
    
   
END SUBROUTINE read_tracers



END MODULE grid_arrays_mod