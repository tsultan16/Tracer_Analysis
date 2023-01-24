PROGRAM solver_test


USE constants_mod
USE fp_solver_mod
USE grid_arrays_mod

IMPLICIT NONE


INTEGER, PARAMETER :: tskip = 1           ! file output interval
INTEGER, PARAMETER :: nsteps = ndumps_tot ! (t_sim/tracer_dump_interval) - 1


REAL(8), ALLOCATABLE :: pp(:), fp2_cc(:), fp2_sum(:)
INTEGER :: i, j, k, ntot
REAL(8) :: my_dt, sim_t, time(0:nsteps), pmean_cc(0:nsteps), pmean_exact(0:nsteps), Bsqr, divV 
LOGICAL :: solve_status

ALLOCATE(pp(0:npbins+1), fp2_cc(1:npbins), fp2_sum(1:npbins))


! read in tracer particle data from reorganized files
CALL create_patch_arrays()


! initilize fp solver
CALL init_fp_solver()


! initialize distribution function and moments
pp = p 
CALL adiabatic_expansion_test_init(pp, fp2_cc)
!CALL diffusion_test_init(pp, fp2_cc)

CALL file_io(0,fp2_cc)


! evolve distribution function
sim_t       = 0.d0
pmean_cc    = 0.D0
pmean_exact = 0.D0
time = 0.D0

!CALL compute_mean_momentum(0,fp2_cc)

DO i = start_dump_num + 1, start_dump_num + nsteps

    my_dt = tracer_dump_interval  ! uniform time step size (i.e. tracer IO dump interval)
    solve_status = .FALSE.
    
    PRINT*,'###############################'   
    PRINT*,'Solving for time dump #',i
    PRINT*,'###############################'   
    
    !**************************
    ! read tracer particle data
    !**************************
    CALL read_tracer_data(i-1)
    
    !*******************************
    ! loop over all tracer particles
    !*******************************

    fp2_sum = 0.D0
    ntot = 0

    !PRINT*,'Begin Fokker-Planck solve...'

    ! just one file for now
    DO j = 1, nfiles

        WRITE(*,FMT ='(" Evolving distribution function for particles in file #",I6)') j
        
        ! just one particle for now
        DO k = 1, file_pars(j)%nparticles
 
            !WRITE(*,FMT ='(" Evolving distribution function for FILE # ",I2, ", Tracer #",I6)') j,k
            !WRITE(*,FMT ='(" Evolving distribution function for Tracer #",I6)') k


            !***********************************************
            ! evolve CR distribution function on each tracer
            !***********************************************
            
            ! get the relevant information saved on the tracer (just divV and Bsqr for now)
            IF((i-1) .EQ. start_dump_num) THEN
                Bsqr = 1.D-6
                divV = 0.D0
                file_pars(j)%particle(k)%fp2_cc = fp2_cc
            ELSE
                Bsqr = file_pars(j)%particle(k)%Bsqr
                divV = file_pars(j)%particle(k)%divV
            END IF
                
            CALL solve(my_dt, file_pars(j)%particle(k)%fp2_cc, Bsqr, divV, solve_status)

            fp2_sum = fp2_sum + file_pars(j)%particle(k)%fp2_cc


        END DO
        ntot = ntot + file_pars(j)%nparticles 
    END DO
        
    ! fp2_sum normalization
    fp2_sum = fp2_sum / ntot

    !PRINT*,'fp2_sum = ',fp2_sum


    ! file IO dump
    !IF(MOD(i, tskip) .EQ. 0) CALL file_io(i)
    
    PRINT*,'Saving to file...'
    CALL file_io(i-start_dump_num,fp2_sum)
    CALL file_io_2(i-start_dump_num,file_pars(1)%particle(1)%x,file_pars(1)%particle(1)%fp2_cc)
    !PRINT*,'Done saving to file.'

    sim_t = sim_t + my_dt
    time(i-start_dump_num) = sim_t
   
    CALL compute_mean_momentum(i-start_dump_num,fp2_sum)
    CALL compute_mean_momentum(i-start_dump_num,file_pars(1)%particle(1)%fp2_cc)
   
    PRINT*,''
    PRINT*,'t, dt = ',i,sim_t, my_dt
    PRINT*,''
   
END DO


! destroy solver
CALL destroy_patch_arrays()

CALL destroy_fp_solver()

DEALLOCATE(pp,fp2_cc,fp2_sum)

PRINT*,''
PRINT*,'Test completed!'
PRINT*,''



CONTAINS


SUBROUTINE file_io(tstep, fp2)

    INTEGER, INTENT(IN) :: tstep
    REAL(8), INTENT(IN) :: fp2(:)
    INTEGER :: i
    CHARACTER(LEN=300) :: filename
    CHARACTER(LEN=6) :: uniti

    WRITE(uniti,FMT='(I5.5)') tstep
    filename = TRIM('Output/fp_tracer_dump=')//TRIM(uniti)//TRIM(".dat")        
    OPEN(UNIT = 17,FILE = filename, FORM = 'UNFORMATTED', ACCESS = 'STREAM')


    !PRINT*,'|----p----|----fp2----|'
    DO i = 1, npbins 
        WRITE(17) pp(i),fp2(i)
        !PRINT*,pp(i), fp2_cc(i)
    END DO

    CLOSE(UNIT=17)


    PRINT*,'File write completed.'

  
END SUBROUTINE file_io

    
SUBROUTINE file_io_2(tstep, x, fp2)

    INTEGER, INTENT(IN) :: tstep, x(3)
    REAL(8), INTENT(IN) :: fp2(:)
    INTEGER :: i
    CHARACTER(LEN=300) :: filename
    CHARACTER(LEN=6) :: uniti

    WRITE(uniti,FMT='(I5.5)') tstep
    filename = TRIM('Output/fp_tracer_id=1_dump=')//TRIM(uniti)//TRIM(".dat")        
    OPEN(UNIT = 17,FILE = filename, FORM = 'UNFORMATTED', ACCESS = 'STREAM')


    !PRINT*,'|----p----|----fp2----|'
    WRITE(17) DBLE(x(:))
    DO i = 1, npbins 
        WRITE(17) pp(i),fp2(i)
        !PRINT*,pp(i), fp2_cc(i)
    END DO

    CLOSE(UNIT=17)


    PRINT*,'File write completed.'

  
END SUBROUTINE file_io_2


  
SUBROUTINE compute_mean_momentum(tstep, fp2)

    INTEGER, INTENT(IN) :: tstep
    REAL(8), INTENT(IN) :: fp2(1:npbins)
    REAL(8) :: p_cc, n_cc, f0_cc, dp, p
    INTEGER :: i
    
    n_cc  = 0.D0
    p_cc  = 0.D0
    f0_cc = 0.D0
    
    DO i = 1, npbins 

        dp = pp(i+1) - pp(i)
        p = 0.5D0 *(pp(i) + pp(i+1))

        f0_cc = 0.5D0 * (fp2(i) + fp2(i+1))

        p_cc = p_cc +  (p**(1)) * f0_cc * dp   
        
        n_cc = n_cc + fp2(i) * delp(i)
        
    END DO
    
    !p_cgmv = p_cgmv/1.D8
    !p_cc = p_cc/1.D8
    
    p_cc = p_cc      !/n_cc
    
    
    pmean_cc(tstep)    = p_cc
    pmean_exact(tstep) = pmean_cc(0) * EXP(4.D0 * sim_t)
 
     
    PRINT*,''
    PRINT*,'N_tot CC   = ',n_cc
    PRINT*,''
    
    
    PRINT*,''
    PRINT*,'Mean momentum CC   = ',p_cc
    PRINT*,''
 
    IF(tstep .EQ. nsteps) THEN
    OPEN(UNIT=10,FILE='Output/pmean.dat', FORM = 'UNFORMATTED', ACCESS = 'STREAM')

    DO i = 0, nsteps 
        WRITE(10) time(i),pmean_cc(i),pmean_exact(i)
    END DO
    
    CLOSE(UNIT=10)        
    
    END IF
    
    
END SUBROUTINE compute_mean_momentum 



END PROGRAM solver_test

