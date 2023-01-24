MODULE grid_arrays_mod

USE constants_mod
IMPLICIT NONE


INTEGER :: nxtot, nytot, nztot, nwx, nwy, nwz, &
           nwavelets_per_scalar, nwavelets_max

!LOGICAL :: particle_presence(1:192*nranks_x*nranks_y*nranks_z) = .FALSE. ! assuming 192 particles per mpi rank at the beginning of the simulation
!LOGICAL :: particle_presence(1:12288*nranks_x*nranks_y*nranks_z) = .FALSE. ! assuming 12288 particles per mpi rank at the beginning of the simulation
LOGICAL :: particle_presence(1:1536*nranks_x*nranks_y*nranks_z) = .FALSE. ! assuming 12288 particles per mpi rank at the beginning of the simulation

TYPE :: patch_scalar

    INTEGER :: nparticles, nwavelet_indices
    REAL(4), ALLOCATABLE :: patch_array(:,:,:,:)
    REAL(4), ALLOCATABLE :: tracer_wavelets(:,:,:)
    REAL(4), ALLOCATABLE :: wavelet_avg_rhoB(:,:)
    INTEGER, ALLOCATABLE :: wavelet_lookup_table(:)
     

END TYPE patch_scalar


TYPE :: patch_metadata

    REAL(4) :: dt, t
    INTEGER :: worldsize(3), patchsize(3), ndomains(3), mpi_rank, rank_coordinates(3), origin_zone(3), &
               nparticles, nwavelet_indices, start_ID, end_ID, init_rank
   
END TYPE patch_metadata


TYPE :: outfile

    INTEGER :: nparticles
    INTEGER :: start_ID, end_ID
    CHARACTER(LEN=200) :: filename

    REAL(4), ALLOCATABLE :: buffer(:,:)

END TYPE outfile


TYPE(patch_scalar), ALLOCATABLE :: dat_rec(:,:,:)
TYPE(patch_metadata), ALLOCATABLE :: file_meta(:,:,:)                    
TYPE(outfile) :: file_pars(1:nfiles)

REAL(4) :: particle_trajectory(0:ndumps_tot,3)

INTEGER :: nmax, mem_bytes, nparticles_tot, nparticles_per_rank, nparticles_per_file


CONTAINS


SUBROUTINE reorganize_files()


    INTEGER :: ndump, i, j, k, start_id, end_id, npfile, xflat, ix, iy, iz
    LOGICAL :: found
    CHARACTER(LEN=6) :: uniti, uintj
    

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
    
      
    ! obtain file meta data for all patches for zeroth time dump
    nparticles_tot = 0   
    CALL read_metadata(start_dump_num, nparticles_tot)             
    
    nparticles_per_rank = nparticles_tot/(nranks_x*nranks_y*nranks_z)
    nparticles_per_file = nparticles_tot/nfiles
    PRINT*,'# of particles per rank = ',nparticles_per_rank

    start_id = 1
    DO i = 1, nfiles        
        file_pars(i)%nparticles = nparticles_per_file
        file_pars(i)%start_ID = start_id        
        end_id   = start_id +  file_pars(i)%nparticles - 1
        start_id = end_id + 1
        file_pars(i)%end_ID = end_id  

        ALLOCATE(file_pars(i)%buffer(file_pars(i)%nparticles,particle_nvars+nwavelets_max*4))  
        
    END DO
    
    
    PRINT*,'# of files = ',nfiles
    PRINT*,''
    DO i = 1, nfiles
        WRITE(*,FMT='("File ",i3,", # of particles = ",i5)') i,file_pars(i)%nparticles  
        PRINT*,'start_id, end_id =',file_pars(i)%start_ID, file_pars(i)%end_ID
    END DO
    PRINT*,''     
    
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Loop over all time dump files (to do..)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO ndump = start_dump_num, start_dump_num + ndumps_tot
   
        PRINT*,'###################################'
        PRINT*,'REORGANIZING FOR TIME DUMP #',ndump
        PRINT*,'% Complete = ',REAL(ndump-start_dump_num)*100.D0/REAL(ndumps_tot)
        PRINT*,'###################################'
   
        WRITE(uintj,FMT='(I0.5)') ndump

        DO i = 1, nfiles
            WRITE(uniti,'(I1.1)') i
            file_pars(i)%filename = TRIM('Output/reorganized_file_')//TRIM(uniti)//TRIM("_dump=")//TRIM(uintj)//TRIM('.dat')
            !PRINT*,'File # ',i        
            !PRINT*,'File Name = ',file_pars(i)%filename        
        END DO

        ! read tracer metadata file
        nparticles_tot = 0
        CALL read_metadata(ndump, nparticles_tot)             
        
        ! allocate tracer data arrays    
        CALL create_patch_arrays()    
        
        found = .FALSE.
        particle_presence = .FALSE.
        
        ! now read tracer data for this time dunmp              
        DO i = 0, nranks_x-1
            DO j = 0, nranks_y-1
                DO k = 0, nranks_z-1
        
                    IF(print_debug) THEN
                        PRINT*,''
                        PRINT*,''
                        WRITE(*,FMT='("Reading tracer data for rank=",i3,", coords = ",i2,i2,i2, ", origin_zone = ", i3,i3,i3)') file_meta(i,j,k)%mpi_rank,file_meta(i,j,k)%rank_coordinates(1), &
                              file_meta(i,j,k)%rank_coordinates(2),file_meta(i,j,k)%rank_coordinates(3),file_meta(i,j,k)%origin_zone(1),file_meta(i,j,k)%origin_zone(2),file_meta(i,j,k)%origin_zone(3)     
                    END IF
                    CALL read_tracers(ndump, i,j,k,file_meta(i,j,k)%mpi_rank, found)
        
                END DO
            END DO
        END DO
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        PRINT*,''
        PRINT*,'ALL(particle_presence) = ',ALL(particle_presence)
        PRINT*,''
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        GO TO 111
        
        PRINT*,''
        DO i = 1, nfiles        

            PRINT*,'#############'        
            PRINT*,'FILE #', i
            PRINT*,'#############'
            
            DO j = 1, file_pars(i)%nparticles
                PRINT*,'Particle ID_flat, grid_pos_flat, Bsrq, divV = ',INT(file_pars(i)%buffer(j,1)),INT(file_pars(i)%buffer(j,2)),file_pars(i)%buffer(j,3),file_pars(i)%buffer(j,4)
                PRINT*,'Wavelet data = ',file_pars(i)%buffer(j,5:particle_nvars+nwavelets_max*4)
            END DO
            PRINT*,''        
            
        END DO
        

        xflat =  file_pars(1)%buffer(1,2)  ! only for particle id#1

        ! collapse flattened particle position index
        ix = MOD(INT(xflat)-1,nxtot+2)
        iy = MOD(INT(xflat)-1,(nxtot+2)*(nytot+2)) / (nxtot+2) 
        iz = (INT(xflat)-1) / ((nxtot+2)*(nytot+2))

        !PRINT*,'PARTICLE # 1, FLATTENED WORLDGRID COORDINATE   = ',file_pars(1)%buffer(1,2)
        !PRINT*,'PARTICLE # 1, UNFLATTENED WORLDGRID COORDINATE = ',ix,iy,iz

        IF((ix .LT. 0 .OR. ix .GT. nxtot+1) .OR. (iy .LT. 0 .OR. iy .GT. nytot+1) .OR. &
           (iz .LT. 0 .OR. iz .GT. nztot+1)) THEN
            PRINT*,'ERROR!!! INVALID PARTICLE WORLD-GRID-COORDINATE!!! ABORTING...'
            STOP
        END IF

        111 CONTINUE    


        ! write to reorganized file        
        DO i = 1, nfiles        
        
            OPEN(UNIT=12, FILE=file_pars(i)%filename, FORM = 'UNFORMATTED', STATUS = 'UNKNOWN', ACCESS = 'STREAM', ACTION='WRITE')
            WRITE(12) file_pars(i)%nparticles
            DO j = 1, file_pars(i)%nparticles
                WRITE(12) file_pars(i)%buffer(j,:)
            END DO
            CLOSE(UNIT=12)
            
        END DO
        
        !#################################################################
 
        xflat =  file_pars(1)%buffer(1,2)  ! only for particle id#1

        ! collapse flattened particle position index
        ix = MOD(INT(xflat)-1,nxtot+2)
        iy = MOD(INT(xflat)-1,(nxtot+2)*(nytot+2)) / (nxtot+2) 
        iz = (INT(xflat)-1) / ((nxtot+2)*(nytot+2))
        
        particle_trajectory(ndump,:) = (/ REAL(ix), REAL(iy), REAL(iz) /)
        
        !#################################################################
            
        mem_bytes = SIZEOF(file_meta) +  SIZEOF(dat_rec)
       
        PRINT*,''
        PRINT*,'Memory allocated for work arrays (Mb) = ', mem_bytes*1.e-6
        PRINT*,''



        ! deallocate tracer data arrays
        CALL destroy_patch_arrays()
        
        IF(single_file) EXIT
        
    END DO


            
    RETURN

    ! write particle trajectory to a file
    OPEN(UNIT=13, FILE='Output/particle_trajectory_id=1.dat', FORM = 'UNFORMATTED', STATUS = 'UNKNOWN', ACCESS = 'STREAM', ACTION='WRITE')
    DO j = 0, ndumps_tot
        WRITE(13) particle_trajectory(j,1),particle_trajectory(j,2),particle_trajectory(j,3)
        !PRINT*,'dump, x, y, z = ',j,particle_trajectory(j,1),particle_trajectory(j,2),particle_trajectory(j,3) 
    END DO
    CLOSE(UNIT=13)
    PRINT*,'Particle id#1 trajectory written to file.'

END SUBROUTINE reorganize_files



SUBROUTINE create_patch_arrays()

    INTEGER :: i, j, k

    ! allocate tracer data arrays    
        DO i = 0, nranks_x-1
            DO j = 0, nranks_y-1
                DO k = 0, nranks_z-1
        
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
        
                    ALLOCATE(dat_rec(i,j,k)%wavelet_lookup_table(1:nx*ny*nz))
                    dat_rec(i,j,k)%wavelet_lookup_table = 0

                END DO
            END DO
        END DO


END SUBROUTINE create_patch_arrays



SUBROUTINE destroy_patch_arrays()

    INTEGER :: i, j, k
    
    DO k = 0, nranks_z-1
        DO j = 0, nranks_y-1
            DO i = 0, nranks_x-1
    
                DEALLOCATE(dat_rec(i,j,k)%patch_array)
                DEALLOCATE(dat_rec(i,j,k)%tracer_wavelets)
                DEALLOCATE(dat_rec(i,j,k)%wavelet_avg_rhoB)
                DEALLOCATE(dat_rec(i,j,k)%wavelet_lookup_table)
                
            END DO
        END DO
    END DO
    
    
    !DO i = 1, nfiles         
    !    DEALLOCATE(file_pars(i)%buffer)  
    !END DO
    
    
    !DEALLOCATE(dat_rec)
    !DEALLOCATE(file_meta)
    
END SUBROUTINE destroy_patch_arrays



SUBROUTINE read_metadata(ndump, nparticles_tot)

    INTEGER, INTENT(IN) :: ndump
    INTEGER, INTENT(OUT) :: nparticles_tot
    CHARACTER(LEN=200) :: filename
    CHARACTER(LEN=6) :: uniti
    INTEGER :: i, j, k, byte_offset, item_bytes, ntot
    REAL(4) :: meta_buff(20)
    
    ntot = 0
    
    ! byte size per data item
    item_bytes = 4
   
    IF(ndump<10) THEN
        WRITE(uniti,'(I1.1)') ndump
    ELSE IF(ndump>=10 .and. ndump<100) THEN
        WRITE(uniti,'(I2.2)') ndump
    ELSE IF(ndump>=100 .and. ndump<1000) THEN
        WRITE (uniti,'(I3.3)') ndump
    ELSE IF(ndump>=1000 .and. ndump<10000) THEN
        WRITE (uniti,'(I4.3)') ndump
    ELSE IF(ndump>=10000 .and. ndump<100000) THEN
        WRITE (uniti,'(I5.3)') ndump  
    END IF

    filename = TRIM(output_filepath)//TRIM('TRACER/tracer_parallel_meta_dump=')//TRIM(uniti)//TRIM('.dat')        
           
    OPEN(UNIT=10, FILE = filename, FORM = 'UNFORMATTED', STATUS = 'OLD', ACCESS = 'STREAM')
   
    byte_offset = 1
  
    DO i = 0, nranks_x-1
        DO j = 0, nranks_y-1
            DO k = 0, nranks_z-1

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

                file_meta(i,j,k)%init_rank = file_meta(i,j,k)%rank_coordinates(1) + nranks_x * file_meta(i,j,k)%rank_coordinates(2) + &
                                             (nranks_x * nranks_y) * file_meta(i,j,k)%rank_coordinates(3) 
                file_meta(i,j,k)%start_ID = 1
                file_meta(i,j,k)%end_ID = file_meta(i,j,k)%nparticles

                nparticles_tot = nparticles_tot + file_meta(i,j,k)%nparticles
    
                IF(print_debug) THEN
                    PRINT*,''
                    PRINT*,'mpi_rank = ',file_meta(i,j,k)%mpi_rank
                    PRINT*,'Rank Coords = ',file_meta(i,j,k)%rank_coordinates
                    PRINT*,'origin zone = ',file_meta(i,j,k)%origin_zone
                    PRINT*,'nparticles = ',file_meta(i,j,k)%nparticles
                    PRINT*,'nwavelet_indices = ',file_meta(i,j,k)%nwavelet_indices
                    PRINT*,'init_rank = ',file_meta(i,j,k)%init_rank                
                    PRINT*,''               
                END IF
            END DO
        END DO
    END DO

    PRINT*,'#############################'
    PRINT*,'Nparticles_total = ',nparticles_tot
    PRINT*,'#############################'
    PRINT*,''
    
    CLOSE(UNIT=10)
   
END SUBROUTINE read_metadata



SUBROUTINE read_tracers(ndump, ri, rj, rk, myrank, found)

    INTEGER, INTENT(IN) :: ndump, ri, rj, rk, myrank
    CHARACTER(LEN=200) :: filename
    CHARACTER(LEN=6) :: uniti
    INTEGER :: i, j, k, isc, item_bytes, byte_offset, ix, iy, iz, pid_flat, grid_pos(3), grid_pos_flat, file_index, ib
    REAL(4) :: pid, init_rank, xflat, divV, Bsqr, tmp  !, buffer(200)
    LOGICAL, INTENT(INOUT) :: found
    LOGICAL :: fexist
    
   
    ! byte size per tracer data item
    item_bytes = 4
        
    IF(ndump<10) THEN
        WRITE(uniti,'(I1.1)') ndump
    ELSE IF(ndump>=10 .and. ndump < 100) THEN
        WRITE(uniti,'(I2.2)') ndump
    ELSE IF(ndump>=100 .and. ndump < 1000) THEN
        WRITE (uniti,'(I3.3)') ndump
    ELSE IF(ndump>=1000 .and. ndump < 10000) THEN
        WRITE (uniti,'(I4.3)') ndump
    ELSE IF(ndump>=10000 .and. ndump < 100000) THEN
        WRITE (uniti,'(I5.3)') ndump  
    END IF

    filename = TRIM(output_filepath)//TRIM('TRACER/tracer_parallel_dump=')//TRIM(uniti)//TRIM('.dat')        
   
    OPEN(UNIT=10, FILE = filename, FORM = 'UNFORMATTED', STATUS = 'OLD', ACCESS = 'STREAM')
   
    
    !READ(10, POS = 1) buffer 
    !PRINT*,''
    !PRINT*,'buffer item#  |  item value '
    !DO i = 1, 200
    !    PRINT*,i,buffer(i)
    !END DO
    !PRINT*,''

   
    ! get the byte offset for this mpi rank
    byte_offset = 1 + myrank * item_bytes
    
    READ(10, POS = byte_offset) tmp  
    byte_offset = 1 + INT(tmp) 
   
    ! read tracer particle data
    DO i = 1, dat_rec(ri,rj,rk)%nparticles  
    
        
        READ(10, POS = byte_offset) pid 
        READ(10, POS = byte_offset + 1*item_bytes) xflat
        READ(10, POS = byte_offset + 2*item_bytes) Bsqr
        READ(10, POS = byte_offset + 3*item_bytes) divV
        byte_offset  = byte_offset + 4*item_bytes

        ! collapse flattened particle position index
        ix = MOD(INT(xflat)-1,nx+2)
        iy = MOD(INT(xflat)-1,(nx+2)*(ny+2)) / (nx+2) 
        iz = (INT(xflat)-1) / ((nx+2)*(ny+2))
  

        !IF(i .EQ. 1) THEN
            !PRINT*,''
            !PRINT*,'PARTICLE# ',i
            !PRINT*,'Particle id = ',INT(pid)
            !PRINT*,'Particle position (flat, i,j,k) = ',INT(xflat),ix,iy,iz
            !PRINT*,'Div(V) = ',divV
            !PRINT*,''    
        !END IF

        ! check for particle position errors
        IF(ix .LT. 0 .OR. ix .GT. nx+1 .OR. iy .LT. 0 .OR. iy .GT. ny+1 .OR. iz .LT. 0 .OR. iz .GT. nz+1) THEN
            PRINT*,'ERROR!! Invalid particle coordinates!! '
            PRINT*,'Particle id, init_rank = ',pid, init_rank
            PRINT*,'Particle position (flat, i,j,k) = ',INT(xflat),ix,iy,iz
            STOP
        END IF
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(.NOT. particle_presence(INT(pid))) THEN
            particle_presence(INT(pid)) = .TRUE.
        ELSE
            PRINT*,'ERROR! Found duplicate particle with global_id = ',INT(pid)
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
        
        
        
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
        ! flattened particle id
        !pid_flat = pid + init_rank * nparticles_per_rank
        
        ! flattened world grid position
        grid_pos(1) = ix + file_meta(ri,rj,rk)%origin_zone(1) - 1
        grid_pos(2) = iy + file_meta(ri,rj,rk)%origin_zone(2) - 1
        grid_pos(3) = iz + file_meta(ri,rj,rk)%origin_zone(3) - 1
        
        grid_pos_flat = (ix + file_meta(ri,rj,rk)%origin_zone(1) - 1) + nxtot * (iy + file_meta(ri,rj,rk)%origin_zone(2) - 1) + &
                        (nxtot * nytot) * (iz + file_meta(ri,rj,rk)%origin_zone(3) - 1)


        !PRINT*,''   
        !PRINT*,'Particle# ',i
        !PRINT*,'Patch position      = ',ix,iy,iz
        !PRINT*,'Rank origin         = ',file_meta(ri,rj,rk)%origin_zone(1),file_meta(ri,rj,rk)%origin_zone(2),file_meta(ri,rj,rk)%origin_zone(3) 
        !PRINT*,'World grid position = ',grid_pos(1),grid_pos(2),grid_pos(3)
        !PRINT*,'Bsqr, divV          = ',Bsqr, divV 
        !PRINT*,''
        !PRINT*,''

        ! add particle data to reorganized file buffer
        IF(pid .LE. file_pars(1)%end_ID) THEN
            file_pars(1)%buffer(pid,1) = pid
            file_pars(1)%buffer(pid,2) = grid_pos_flat
            file_pars(1)%buffer(pid,3) = Bsqr
            file_pars(1)%buffer(pid,4) = divV
            ib = 5
            DO j = 1, nwavelets_per_scalar
                file_pars(1)%buffer(pid,ib:ib+nscalars) = dat_rec(ri,rj,rk)%tracer_wavelets(i,j,1:2+nscalars-1)
                ib = ib + 1 + nscalars
            END DO    
        ELSE
            file_pars(2)%buffer(pid-file_pars(1)%nparticles,1) = pid
            file_pars(2)%buffer(pid-file_pars(1)%nparticles,2) = grid_pos_flat
            file_pars(2)%buffer(pid-file_pars(1)%nparticles,3) = Bsqr
            file_pars(2)%buffer(pid-file_pars(1)%nparticles,4) = divV
            ib = 5
            DO j = 1, nwavelets_per_scalar
                file_pars(2)%buffer(pid-file_pars(1)%nparticles,ib:ib+nscalars) = dat_rec(ri,rj,rk)%tracer_wavelets(i,j,1:2+nscalars-1)
                ib = ib + 1 + nscalars
            END DO             
        END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        IF(INT(pid) .EQ. 10)THEN
           
            !PRINT*,''
            !PRINT*,'%%%%%%%%%%%%%%%'
            !PRINT*,'Particle found.'
            !PRINT*,'%%%%%%%%%%%%%%%'
            !PRINT*,''

            !PRINT*,'Particle global_id = ',INT(pid)
            !PRINT*,'Particle position (flat, i,j,k) = ',INT(xflat),ix,iy,iz
            !PRINT*,'Div(V) = ',divV       
        
            IF(.NOT. found) THEN
                ! append tracer data to a unique file
                !INQUIRE(FILE="tracer_100_0.dat", EXIST=fexist)
                
                !PRINT*,'fexist = ',fexist
                
                !IF (fexist) THEN
                    !OPEN(UNIT=12, FILE="Output/tracer_6_1.dat", FORM = 'UNFORMATTED', STATUS = 'UNKNOWN', ACCESS = 'STREAM', POSITION='APPEND', ACTION='WRITE')
                    !WRITE(12) file_meta(ri,rj,rk)%t,file_meta(ri,rj,rk)%dt,pid,init_rank,ri,rj,rk,ix,iy,iz
                !ELSE
                !    OPEN(UNIT=12, FILE="Output/tracer_1_0.dat", FORM = 'UNFORMATTED', STATUS = 'NEW', ACCESS = 'STREAM', ACTION='WRITE')
                !    WRITE(12) file_meta(ri,rj,rk)%t,file_meta(ri,rj,rk)%dt,pid,init_rank,ri,rj,rk,ix,iy,iz
                !END IF

                !CLOSE(12)   
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