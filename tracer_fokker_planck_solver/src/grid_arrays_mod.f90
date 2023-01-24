MODULE grid_arrays_mod

USE constants_mod
IMPLICIT NONE


INTEGER :: nxtot, nytot, nztot, nwx, nwy, nwz, np_tot


TYPE :: tracer
    
    INTEGER :: tid, x(3)  ! tracer ID, origin mpi rank, tracer cell indices 
    REAL(4) :: divV, Bsqr                  ! velocity divergence
    REAL(4) :: Bavg(3), rhoavg       ! average density and magnetic field
    
    REAL(8), ALLOCATABLE :: fp2_cc(:) ! distribution function and it's first and second moments
    
END TYPE tracer



TYPE :: outfile

    INTEGER :: nparticles
    INTEGER :: start_ID, end_ID
    CHARACTER(LEN=200) :: filename

    REAL(4), ALLOCATABLE :: buffer(:,:)
    TYPE(tracer), ALLOCATABLE :: particle(:)


END TYPE outfile                   


TYPE(outfile) :: file_pars(1:nfiles)


INTEGER :: nmax, mem_bytes

    ! max number of wavelets per tracer
INTEGER, PARAMETER :: nwavelets_per_scalar =  nwavelets_per_bin * (1+2*bin_neighborhood)**3
INTEGER, PARAMETER :: nwavelets_max = nwavelets_per_scalar
    

CONTAINS


SUBROUTINE create_patch_arrays()


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


    !CALL read_tracer_data(0)

    
    !mem_bytes = 22*SIZEOF(vxk_re) + SIZEOF(fx) + SIZEOF(vp) + 4*SIZEOF(buffer) + SIZEOF(dat_rec)
   
    !PRINT*,''
    !PRINT*,'Memory allocated for work arrays (Mb) = ', mem_bytes*1.e-6
    !PRINT*,''


END SUBROUTINE create_patch_arrays


SUBROUTINE destroy_patch_arrays()

    INTEGER :: i, j, k

    DO i = 1, nfiles
        DEALLOCATE(file_pars(i)%buffer)
        DO j = 1, file_pars(i)%nparticles
            DEALLOCATE(file_pars(i)%particle(j)%fp2_cc)    
        END DO
        DEALLOCATE(file_pars(i)%particle)
    END DO
  

END SUBROUTINE destroy_patch_arrays



SUBROUTINE read_tracer_data(ndump)


    INTEGER :: i, j, k, ndump, xflat, ix, iy, iz  
    CHARACTER(LEN=6) :: uniti, uintj
    REAL(4) :: particle_buffer(1:particle_nvars+nwavelets_max*4)
    
    PRINT*,'#################################################'
    PRINT*,'READING REORGANIZED FILES FOR TIME DUMP #',ndump
    PRINT*,'#################################################'

    WRITE(uintj,FMT='(I0.5)') ndump


    ! loop over redistributed files
    DO i = 1, nfiles
    
        WRITE(uniti,'(I1.1)') i
        file_pars(i)%filename = TRIM(reorganized_filepath)//TRIM('reorganized_file_')//TRIM(uniti)//TRIM("_dump=")//TRIM(uintj)//TRIM('.dat')

        OPEN(UNIT=12, FILE=file_pars(i)%filename, FORM = 'UNFORMATTED', STATUS = 'OLD', ACCESS = 'STREAM', ACTION='READ')             
        READ(12) file_pars(i)%nparticles

        PRINT*,''
        PRINT*,'File # ',i
        PRINT*,'Nparticles = ',file_pars(i)%nparticles
        PRINT*,''


        IF(ndump .EQ. start_dump_num) THEN 
            ALLOCATE(file_pars(i)%buffer(file_pars(i)%nparticles,particle_nvars+nwavelets_max*4)) 
            ALLOCATE(file_pars(i)%particle(file_pars(i)%nparticles))

            DO j = 1, file_pars(i)%nparticles
                ALLOCATE(file_pars(i)%particle(j)%fp2_cc(1:npbins))    
            END DO

        END IF    

        DO j = 1, file_pars(i)%nparticles
            
            READ(12) particle_buffer(:)
            
            file_pars(i)%buffer(j,:) = particle_buffer(:)
            
            xflat = file_pars(i)%buffer(j,2)
            
            ! collapse flattened particle position index
            ix = MOD(INT(xflat)-1,nxtot+2)
            iy = MOD(INT(xflat)-1,(nxtot+2)*(nytot+2)) / (nxtot+2) 
            iz = (INT(xflat)-1) / ((nxtot+2)*(nytot+2))
            
            !IF(i .EQ. 1 .AND. j .EQ. 1) THEN
            !    PRINT*,'PARTICLE # 1, FLATTENED WORLDGRID COORDINATE = ',xflat !file_pars(i)%buffer(j,2)
            !    PRINT*,'PARTICLE # 1, UNFLATTENED WORLDGRID COORDINATE = ',ix,iy,iz
            !END IF


            IF((ix .LT. 0 .OR. ix .GT. nxtot+1) .OR. (iy .LT. 0 .OR. iy .GT. nytot+1) .OR. &
               (iz .LT. 0 .OR. iz .GT. nztot+1)) THEN
                PRINT*,'ERROR!!! INVALID PARTICLE WORLD-GRID-COORDINATE!!! ABORTING...'
                PRINT*,'PARTICLE # ',j
                PRINT*,'FLATTENED WORLDGRID COORDINATE = ',xflat 
                PRINT*,'UNFLATTENED WORLDGRID COORDINATE = ',ix,iy,iz
            
                STOP
            END IF

            !PRINT*,''
            !PRINT*,'Particle # ',j
            !PRINT*,'Particle ID_flat, grid_pos_flat = ',INT(file_pars(i)%buffer(j,1)),INT(file_pars(i)%buffer(j,2))
            !WRITE(*,FMT='(" Grid Position unflattened       = ",i5, i5, i5)') ix, iy, iz
            !PRINT*,'Bsqr, divV                      = ',file_pars(i)%buffer(j,3),file_pars(i)%buffer(j,4)
            !PRINT*,'particle buffer = ',particle_buffer
            !PRINT*,''
            
            file_pars(i)%particle(j)%x    = (/ ix, iy, iz /)            
            file_pars(i)%particle(j)%tid  = file_pars(i)%buffer(j,1)
            file_pars(i)%particle(j)%Bsqr = file_pars(i)%buffer(j,3)
            file_pars(i)%particle(j)%divV = file_pars(i)%buffer(j,4)

        END DO
        
        CLOSE(UNIT=12)

        PRINT*,'Done reading from file # ',i

    END DO
    

END SUBROUTINE read_tracer_data




END MODULE grid_arrays_mod