MODULE readfile_mod

USE constants_mod
IMPLICIT NONE

CONTAINS



SUBROUTINE readfile_wombat(t_dump, field_array)


    INTEGER, INTENT(IN) :: t_dump
    REAL(4), INTENT(INOUT) :: field_array(:,:,:,:) ! field array single precision by default
    INTEGER :: i, j, k, iv, ib, ri, rj, rk, &
               ilow, ihi, jlow, jhi, klow, khi, izone_low, izone_hi, jzone_low, jzone_hi, kzone_low, kzone_hi
    INTEGER :: item_bytesize, byte_offset, buff_size
    CHARACTER(LEN=300) :: filename
    CHARACTER(LEN=6) :: uniti
    REAL(4), ALLOCATABLE :: file_buff(:) 
    
    
    buff_size = (nranks_x * nranks_y * nranks_z) * ((nx+2*nb) * (ny+2*nb) * (nz+2*nb) * wombat_dump_nvars)
    ALLOCATE(file_buff(buff_size))

    item_bytesize = 4    
    
    
    WRITE(uniti,FMT='(I4.4)') t_dump

    filename = TRIM(output_filepath)//TRIM('MHD/compr-')//TRIM(uniti)//TRIM('-part00-000')        

    PRINT*,''
    PRINT*,'Reading wombat dump # ',t_dump
    PRINT*,'Filename: ',filename
    PRINT*,''
       
    
    OPEN(UNIT=17,FILE=filename, FORM = 'UNFORMATTED', STATUS = 'OLD', ACCESS = 'STREAM')
            
    byte_offset = 1
    
    READ(17, POS = 1) file_buff
    CLOSE(UNIT=17)

    ib = 1
    
    DO ri = 0, nranks_x-1 
    DO rj = 0, nranks_y-1
    DO rk = 0, nranks_z-1
    
        
        ilow = 1 - nb + nx * ri 
        ihi  = ilow + nx + 2*nb - 1
        jlow = 1 - nb + ny * rj
        jhi  = jlow + ny + 2*nb - 1
        klow = 1 - nb + nz * rk
        khi  = klow + nz + 2*nb - 1
                   
        izone_low = 1 + nx * ri
        jzone_low = 1 + ny * rj
        kzone_low = 1 + nz * rk
        
        izone_hi = izone_low + nx - 1
        jzone_hi = jzone_low + ny - 1
        kzone_hi = kzone_low + nz - 1


        DO iv = 1, wombat_dump_nvars        
        DO k = klow, khi 
        DO j = jlow, jhi
        DO i = ilow, ihi
           
                IF(i .GE. izone_low .AND. i.LE. izone_hi .AND. j .GE. jzone_low .AND. j .LE. jzone_hi .AND.  k .GE. kzone_low .AND. k .LE. kzone_hi) THEN   
                    field_array(i,j,k,iv) = file_buff(ib)
                END IF
                ib = ib + 1        
                 
        END DO       
        END DO    
        END DO
        END DO

    
    END DO    
    END DO
    END DO    
    
    
    DEALLOCATE(file_buff)
    
    
    
END SUBROUTINE readfile_wombat



END MODULE readfile_mod 