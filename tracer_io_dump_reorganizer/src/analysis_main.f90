PROGRAM analysis

USE constants_mod
USE grid_arrays_mod
USE OMP_LIB

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


! initialize grid arrays and reorganize tracer data
CALL reorganize_files()




! grid size
Lx = 1.d0
Ly = 1.d0
Lz = 1.d0
dx = Lx/DBLE(nranks_x*nx)


t1 = OMP_GET_WTIME()


! perform tracer wavelet reconstruction of the velocity field

t2 =  OMP_GET_WTIME()

!CALL tracer_reconstruct()

t3 = OMP_GET_WTIME()



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


!CALL destroy_patch_arrays()

    
PRINT*,'Done!'
PRINT*,''


!##########################################################################################


CONTAINS



END PROGRAM analysis