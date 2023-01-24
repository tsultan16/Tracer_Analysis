MODULE constants_mod

IMPLICIT NONE


!********************************
! simulation and grid parameters
!********************************
INTEGER, PARAMETER :: nx       = 32    ! grid size per MPI rank (or in this case, per patch)
INTEGER, PARAMETER :: ny       = 32    	     
INTEGER, PARAMETER :: nz       = 32  	     
INTEGER, PARAMETER :: nb       = 6     ! number of boundary cells
INTEGER, PARAMETER :: dump_num = 0 !50


INTEGER, PARAMETER :: nwavelets_per_bin = 30 ! # of wavelets per tracer per bin
INTEGER, PARAMETER :: bin_neighborhood = 0
INTEGER, PARAMETER :: nscalars = 3   ! # of scalars (wavelet amplitudes for vx,vy,vz)

LOGICAL, PARAMETER :: print_debug     = .FALSE.
LOGICAL, PARAMETER :: override_checks = .TRUE.

LOGICAL, PARAMETER :: output_basis = .FALSE.
LOGICAL, PARAMETER :: keep_bndry   = .FALSE.

REAL, PARAMETER :: FOURPI = 16.d0*ATAN(1.d0)
REAL, PARAMETER :: TWOPI  = 08.d0*ATAN(1.d0)

!*************************************
! MPI domain decomposition parameters
!*************************************
INTEGER, PARAMETER :: nranks_x = 3  ! #  of ranks along x direction
INTEGER, PARAMETER :: nranks_y = 3  ! # of ranks along y direction
INTEGER, PARAMETER :: nranks_z = 3  ! # of ranks along z direction

!********************
! physics parameters
!********************
REAL :: sound_speed = 1.d0  ! constant isothermal sound speed

!************************
! File output Parameters
!************************
CHARACTER(LEN=200), PARAMETER :: output_filepath = '/data/uchu/tanzid/Wombat_Sprint_5/WAVELETS/wombat/build/Snapshots/'     
REAL, PARAMETER :: dump_frequency = 0.5    
INTEGER, PARAMETER :: wombat_dump_nvars = 7
LOGICAL, PARAMETER :: tracer_rec_file = .TRUE.

END MODULE constants_mod