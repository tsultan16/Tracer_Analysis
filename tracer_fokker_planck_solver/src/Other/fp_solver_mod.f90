MODULE fp_solver_mod

! FILENAME: cosmicrays.f90
!
! DESCRIPTION: This module contains an implementation of the CGMV scheme for solving the Lagrangiam cosmic ray 
!              transport equation (i.e. momentum space evolution of the isotropic CR distribution function) 
!
!
!              This implementation only considers a passive population of cosmic ray electrons or protons 
!              and does NOT include spatial diffusion.
!
!              Momentum Diffusion Co-effeicient: Typically Dpp ~ p^2/tacc, where tacc is a characrteristic acceleration time-scale
!                                                which may or may not be momentum dependent. 
!
!              Slope reconstruction: The momentum-fluxes for updating the n and g moments require intra-bin slopes of the piecewise 
!                                    power law distribution function. The intra-bin slope can be computed by solving for the root of the
!                                    non-linear function:  gnp(q) = (q-3)*(1-d^(4-q)) / (q-4)*(1-d^(3-q))  
!                                    using an iterative root finder algorithm (such as Newton-Raphson).  
!                                      
!                                    Alternatively, the smoothness and monotonicity of this function make it advantageous to use an   
!                                    approximate non-iterative root solving approach via interpolation. We can sample this function                 
!                                    within a physically relevant domain, e.g. gnp in [gnp(q=-15),gnp(q=15)], then interpolate between
!                                    these samples. Piecewise cubic spline interpolation, using as few as only 3 samples (a.k.a. "knot points"), 
!                                    yields excellent approximations. The other added benefit of computing the slopes this way is that it can be 
!                                    done inside a vector loop. I've included subroutines for computing the slope using both approaches.                                                                               
!
!
! TODO: - Add protection floors inside root solver slope calculation (doesn't matter if we end up using the non-iterative root solver)
!       
!
! DEPENDENCIES: constants_mod
!
! REFERENCES: Jones and Kang (2005), Chang and Cooper (1970)
!
! HISTORY:
!    5/?/2021  - Tanzid Sultan  
!

USE constants_mod


! local variables (scratch arrays)

REAL(8), ALLOCATABLE :: work_np(:), work_gp(:), work_fp(:) , &
                        p(:), q(:), Fn(:), Fg(:), slopes(:), &
                        dn(:), dg(:), delp(:), delq(:)

REAL(8) :: gnp_knots(3), q_knots(3), qp_knots(3) ! arrays for hermite cubic spline knot points
REAL(8) :: dx, dt, dlp, dp, pmin, pmax, dgp, gnpmin, gnpmax, np_sum, gp_sum 

! Solver Parameters
REAL(8),    PARAMETER :: third   = 1.D0 / 3.D0
REAL(8),    PARAMETER :: ifourpi = 1.D0 / fourpi
REAL(8),    PARAMETER :: ZERO    = 1.D-3
REAL(8),    PARAMETER :: COUR_p  = 0.3D0            ! momentum convection courant number
INTEGER,    PARAMETER :: maxIter = 50               ! maximum number of iterations allowed for Newton-Raphson
REAL(8),    PARAMETER :: alamb   = 2.D0             ! constant greater than unity neded for injection momentum (JK05) FF injection fraction (JK05), default set this to 2.0
REAL(8),    PARAMETER :: ceps1   = 0.D0             ! flux fraction epsilon parameter (JK05)
REAL(8),    PARAMETER :: shkmin  = 1.15D0           ! minimum shock mach number for FF injection (JK05)
REAL(8),    PARAMETER :: qmax    = 15.D0            ! slope cutoff for iterative root solver




REAL(8) :: gamma, c_light, beta0, m_e, m_p, sigmaT, tsc, itacc, &
             sigmaSB, redshift, Tcmb, Ucmb, Urad, nmin

LOGICAL :: CREloss, injOn 



CONTAINS


SUBROUTINE init_fp_solver()


    !gamma = gamma_in

    ! set speed of light 
    !c_light = simunits%c_light ! speed of light in cm/s

    ! set inverse speed of light in code units
    !beta0 = simunits%velocity_cgs / c_light
   
    ! set electron and proton masses 
    !m_e = simunits%m_e
    !m_p = simunits%m_p 

    ! set Thomson cross-section 
    !sigmaT = simunits%sigmaT
    
    ! other useful constants
    !sigmaSB = simunits%sigmaSB ! Stephan-Boltzmann constant
    !redshift = 0.D0            ! redshift of computational domain 
    !Tcmb = simunits%Tcmb       ! WMAP5 CMB temperature
    
    ! set syncrhotron characteristic lifetime (in code units) of a particle of momentum (m_p*c) in a magnetic field of strength 1 code units
    !tsc =  (3.D0 / 4.D0) * (m_e / m_p)**2 *  &
    !            (m_p * c_light / sigmaT) * 8.d0 * simunits%pi / simunits%magnetic_cgs**2  ! in seconds
   
    !tsc = tsc / simunits%time_cgs ! convert into code units

    ! set inverse Fermi-II acceleration characteristic time-scale (in code units)
    itacc = 0.D0
    
    ! set CMB and other background magnetic field energy densities into code units  
    !Ucmb = 4.D0 * sigmaSB * ((1.d0 + redshift) * Tcmb)**4 / (c_light * simunits%energy_density_cgs)
    !Urad = Urad / simunits%energy_density_cgs


    ! radiative loss setting
    CREloss = .FALSE.  ! hard coded for now...

    ! particle injection setting
    injOn = .FALSE.

    ! protection floor for n moment
    nmin = 1.D-24
 
    ! p-space bounds (hardcoded for now)
    pmin  = 1.D-1
    pmax  = 1.D10

    dlp = (LOG(pmax)-LOG(pmin))/npbins  ! uniform logarithmic spacing 
    dp  = EXP(dlp)  

    ! allocate memory for all work arrays
    CALL create_workpool()
    
    ! prepare array of cubic spline knot points
    CALL set_spline_knot_points()
  
    PRINT*,''
    PRINT*,"CGMV Parameters:"
    PRINT*,"    nbins     =  ",npbins
    PRINT*,"    pmin      =  ",pmin
    PRINT*,"    pmax      =  ",pmax
    PRINT*,"    dp        =  ",dp
    PRINT*,"    d logp    =  ",dlp
    PRINT*,''

    PRINT*,''
    PRINT*,'Fokker-Planck solver initilized..'
    PRINT*,''


END SUBROUTINE init_fp_solver


SUBROUTINE destroy_fp_solver()

    CALL destroy_workpool()

    PRINT*,''
    PRINT*,'Fokker-Planck solver destroyed..'
    PRINT*,''    

END SUBROUTINE destroy_fp_solver


! Allocate memory for all scratch arrays and initialize
SUBROUTINE create_workpool()

    INTEGER :: i
          
    ALLOCATE(work_np(1:npbins),work_gp(1:npbins),work_fp(1:npbins),p(0:npbins+1),Fn(1:npbins+1),Fg(1:npbins+1), slopes(0:npbins+1), &
             dn(1:npbins), dg(1:npbins), delp(0:npbins+1), delq(0:npbins+1),q(0:npbins+1))        
  
    work_np = 0.D0
    work_gp = 0.D0
    work_fp = 0.D0
    p       = 0.D0
    q       = 0.D0
    Fn      = 0.D0
    Fg      = 0.D0
    delp    = 0.D0
    delq    = 0.D0
    
    DO i = 0, npbins+1
        p(i) = pmin * dp**(i)  ! value of p at left-edge of ith bin
    END DO

    DO i = 0, npbins
        q(i) = 0.5D0 * (p(i+1) + p(i))  ! q_i = p_i+1/2 = 0.5*(p_i+1 + p_i)
    END DO
    
    DO i = 1, npbins
        delp(i) = 0.5D0 * (p(i+1) - p(i-1)) ! delp_i = 0.5* (p_i+1 - p_i-1)
        delq(i) = p(i+1) - p(i)           ! delq_i = p_i+1 - p_i
    END DO
    
    
END SUBROUTINE create_workpool



! Deallocate memory for all scratch arrays
SUBROUTINE destroy_workpool()

    DEALLOCATE(work_np, work_gp, work_fp, p, q, Fn, Fg, slopes, dn, dg, delp, delq)
    

END SUBROUTINE destroy_workpool



! top-level interface for performing a full CGMV update 
SUBROUTINE solve(dt_in, fp, np, gp, fp2_cc, update_cmplt)

    REAL(8)  :: dt_in, fp(:), np(:), gp(:), fp2_cc(:)
    LOGICAL, INTENT(OUT) :: update_cmplt
    REAL(8) :: alpha, N0
    
    dt = dt_in
    
    ! update particle distribution for adiabatic changes and losses
    !CALL momentum_convection_3d(fp, np, gp)
    
    CALL chang_cooper(fp2_cc)
    
    CALL diagnostics()
    
    
    update_cmplt = .TRUE.    


END SUBROUTINE solve



SUBROUTINE diagnostics()

    INTEGER :: i
    REAL(8) :: del_np_sum, del_gp_sum

    del_np_sum = np_sum
    del_gp_sum = gp_sum
    
    np_sum = 0.D0
    gp_sum = 0.D0

    DO i = 1, npbins
        np_sum = np_sum  + work_np(i)    
        gp_sum = gp_sum  + work_gp(i)     
    END DO

    del_np_sum = ABS(np_sum - del_np_sum)/np_sum
    del_gp_sum = ABS(gp_sum - del_gp_sum)/gp_sum


    PRINT*,''
    WRITE(*,FMT='("np_sum = ",E11.4,", fractional del_np_sum",E11.4)') np_sum, del_np_sum
    WRITE(*,FMT='("gp_sum = ",E11.4,", fractional del_gp_sum",E11.4)') gp_sum, del_gp_sum
    PRINT*,''

END SUBROUTINE diagnostics


!****************************
!    Chang Cooper Scheme    !
!****************************


! Updates momentum convection part using Chang-Cooper Scheme
SUBROUTINE chang_cooper(fp2)

    REAL(8), INTENT(INOUT) :: fp2(:)
    REAL(8) :: pdot, pdot_max, pdot_loss, subdt, Bsqr
    INTEGER :: i, nsub, n, ip, j
    REAL(8) :: H(0:npbins+1), Wp(0:npbins+1), Wm(0:npbins+1)
    REAL(8) :: Dpp(1:npbins), Qp(1:npbins)
    REAL(8) :: A(0:npbins+1), B(0:npbins+1), C(0:npbins+1), r(0:npbins+1)
    REAL(8) :: gam(0:npbins+1)
    REAL(8) :: wi, expwi, bet, divV
    

    divV = 0.D0

    ! compute moomentum velocity (pdot) from radiative losses
    Bsqr = 0.D0 !SUM(workp%grid1d(i,5:7)**2)  ! local magnetic field strength squared
    pdot_loss = 0.D0 ! self%radiative_losses(Bsqr)
    
    ! compute required time step size with CFL number = 0.8 
    subdt = dt !COUR_p * pdot_max / tls%dlp  
    
    ! compute number of subcycle steps
    nsub = 1      
    
    
    ! Dpp_i+1/2
    DO i = 1, npbins
        Dpp(i) = 0.5D0 * (p(i)**2 + p(i+1)**2)  
    END DO
    Dpp(:) = Dpp * 1.D-4
    
    
    ! no injection
    Qp  = 0.D0   

    
    ! loop over subcycles
    DO n = 1, nsub       
                  
        !***************************************************
        ! Compute H_i+1/2, Wp_i+1/2, Wm_i+1/2, A_i, B_i, C_i
        !***************************************************        
        
        H = 0.D0
        DO ip = 1, npbins-1            
            H(ip) = -2.D0 * Dpp(ip) / q(ip) + third * divV * q(ip)  ! only adiabatic and diffusion terms, not including radiative losses
        END DO 
        
        Wp = 0.D0
        Wm = 0.D0
        DO ip = 2, npbins-1
                        
            wi =  MIN(700.D0, SIGN(1.D0, H(ip) * delq(ip) / Dpp(ip)) * &
                  MAX(1.D-8 , ABS(H(ip) * delq(ip) / Dpp(ip))))   ! wi restricted between (10^-8, 700)
                                                                  
            expwi = EXP(wi)     
                                
            Wp(ip) = 1.D0 / (1.D0 - 1.D0 / expwi)                                     
            Wm(ip) = 1.D0 / (expwi - 1.D0)
         
        END DO
        
        A = 0.D0
        B = 0.D0
        C = 0.D0
        r = 0.D0
        DO ip = 1, npbins           
            A(ip) = - subdt * H(ip)   * Wp(ip)   / delp(ip)          
            C(ip) = - subdt * H(ip-1) * Wm(ip-1) / delp(ip)          
            B(ip) = 1.D0 + subdt * (H(ip-1) * Wp(ip-1) + H(ip) * Wm(ip)) / delp(ip)  ! + particle loss term goes here        
            r(ip) = fp2(ip) + subdt * Qp(ip)
        END DO
                  
        
        !********************************
        ! Update fp2 := p^2 f   
        !********************************    
        
        ! Tridiagonal matrix inversion (using back substitution)
        bet = B(1)
        fp2(1) = r(1) / bet

        DO ip = 2, npbins
            gam(ip) = A(ip-1) / bet
            bet = B(ip)- C(ip) * gam(ip)
            fp2(ip) = (r(ip) - C(ip) * fp2(ip-1)) / bet
        END DO
        
        DO ip = npbins-1, 1, -1 
            fp2(ip) = fp2(ip) - gam(ip+1) * fp2(ip+1) 
        END DO
        
    END DO

       

END SUBROUTINE chang_cooper


!****************************
!        CGMV Scheme        !
!****************************

SUBROUTINE momentum_convection_3d(fp, np, gp)

    REAL(8), INTENT(INOUT) :: np(:), gp(:), fp(:)
    REAL(8) :: pdot, pdot_max, pdot_loss, subdt, Bsqrl, divV, Bsqr
    INTEGER :: i, j, k, nsub, n

    ! get velocity divergence
    divV = 0.D0
    
    ! compute moomentum velocity (pdot) from radiative losses
    Bsqr = 0.D0   ! local magnetic field strength squared
    pdot_loss = 0.D0 !radiative_losses(Bsqr)

    !**************************************************************************
    ! Check if time step satisfies courant condition for CR momentum convection
    ! If not, we'll need to subcycle for stability   
    !**************************************************************************

    ! compute maximum momentum "velocity"
    pdot_max = third * ABS(divV) + p(npbins) * ABS(pdot_loss) + ABS(MAXVAL(ABS(slopes)) * 0.D0 ) 

    ! compute required time step size with the given CFL number  
    subdt = MIN(COUR_p * dlp / pdot_max, dt)  


    !PRINT*,'******************************************'
    !PRINT*,'  dt, subdt =  ',dt,subdt
    !PRINT*,'******************************************'


    ! compute number of subcycle steps
    nsub = CEILING(dt / subdt)         

    ! set sub step size
    subdt = dt/nsub

    ! copy into work array 
    work_np = np
    work_gp = gp
    work_fp = fp

    PRINT*,'# of subcycles = ',nsub

    ! loop over subcycles
    DO n = 1, nsub       

        !PRINT*,'Subcycle step #',n

        !*************************
        ! Compute intra-bin slopes
        !*************************
        !PRINT*,'Computing slopes ...'        
        CALL compute_slopes()
        !CALL compute_slopes_approx()

        !***************************************
        ! Compute upwinded momentum space fluxes
        !***************************************
        !PRINT*,'Computing momentum fluxes ...'
        CALL compute_momentum_fluxes(divV, pdot_loss)


        !***************
        ! Update n and g    
        !***************    
        !PRINT*,'Updating distribution function ...'
        CALL update_ng(subdt, divV, pdot_loss)                      
                                              
    END DO

    ! copyback from work arrays
    np = work_np
    gp = work_gp
    fp = work_fp


END SUBROUTINE momentum_convection_3d


! this subroutine advances n and g by one time-step
SUBROUTINE update_ng(subdt, divV, pdot_loss)

    REAL(8), INTENT(IN) :: subdt, divV, pdot_loss
    INTEGER :: ip, ipg
    REAL(8) :: lossterm, diffterm, Dpp, q3, q4, q5, d4, d5, num, den

    ! n and g updates done on separate loops (apparently no vectorization if done in the same loop)
    DO ip = 1, npbins
        
        work_np(ip) = work_np(ip) - subdt*(Fn(ip+1) - Fn(ip))  
        
        ! apply protection floor
        work_np(ip) = MAX(nmin * p(ip)**(2.1D0), work_np(ip))
        
    END DO
            
    
    DO ip = 1, npbins   

        ! Compute the extra loss term: Integral_p(i)^p(i+1) dp * pdot_loss * p^2 * f
        ! Since pdot_loss ~ p^2 * sigmaT * (UB+ Ucmb)/(m_e^2 c^2) => we need to evaluate an integral of p^4*f 
        ! i.e. the next higher moment after g, which = fi * pi^5 *(1- di^(5-qi)) / (qi-5)  = p5f *(1-d5) / q5
        !                                            = p * (gi*q4 / (1-d4) ) * (1-d5) / q5 = p * gi * q4* (1-d5) / (q5 * (1-d4))
        
        lossterm = 0.D0 
        IF(CREloss .AND. pdot_loss .NE. 0.D0)THEN
            
            q4  = slopes(ip) - 4.D0
            q5  = slopes(ip) - 5.D0
            d4  = EXP(-q4*dlp)
            d5  = dp * d4           
                    
            num = q4 * (1.D0 - d5)  
            dem = q5 * (1.D0 - d4)  

            ! When q equal to 5, loss terms becomes singular. If that happens, replace num and den with appropriate limits.
            IF(ABS(q5) .LT. 1.D-3)THEN
                num = (1.D0 - d5) + q4 * q5 * d4
                den = (1.D0 - d4) + q4 * q5 * EXP((3.D0 - slopes(ip)) * dlp)
            END IF

            lossterm = pdot_loss * p(ip) * num / den

        END IF

        ! Compute the weighted diffusion term: Integral_p(i)^p(i+1) dp * D * q * p * f  
        ! D = itacc*p^2 for now, so that this integral ~ g moment. In general, this will depend on the precise functional form of D(p)
        
        
        ! set momentum diffusion co-efficient D (here Dpp = D / p^2)
        Dpp = p(ip)**2 ! 1.D0 ! itacc  ! D is hard-coded to be itacc*p^2 for now        
        
        diffterm = slopes(ip) * work_gp(ip)  ! diff term is the integral of (q_i * p * Dpp * f(p)) over the momentum bin


        work_gp(ip) = work_gp(ip) - subdt * &
                      (Fg(ip+1) - Fg(ip) + (third * divV - lossterm) * work_gp(ip) - diffterm )
                
        ! apply protection floor
        work_gp(ip) = MAX(nmin * p(ip)**(1.1D0), work_gp(ip))
                
    END DO
           
    ! Now update f(p)       
    DO ip = 1, npbins              

        ! this seems like a good place to update the grid cosmic ray distribution function..
        q3 = slopes(ip) - 3.D0
        work_fp(ip) = work_np(ip) * q3 / (p(ip)**3 * (1.D0 - EXP(-q3*dlp)))                 

    END DO
    

END SUBROUTINE update_ng


! This subroutine computes momentum fluxes through bin edges over a time-step
SUBROUTINE compute_momentum_fluxes(divV, pdot_loss)

    REAL(8), INTENT(IN) :: divV
    REAL(8), INTENT(IN) :: pdot_loss
    REAL(8) :: pdot, pdot_max, Dpp
    REAL(8) :: q3, q4, d3, d4, p3f, p4f, p5f, avg_slope
    INTEGER :: ip

    ! apply no-flux boundary condition at lowest bin
    Fn(1) = 0.D0
    Fg(1) = 0.D0

    ! Loop over momentum bins and compute momentum fluxes for n and g                
    DO ip = 2, npbins+1         

        ! set momentum diffusion co-efficient
        Dpp = p(ip)**2 
            
        ! compute total pdot term with a 'p' factored out from the numerator   
        !avg_slope = 0.5D0 * (slopes(ip-1) + slopes(ip))        
        avg_slope = 0.5D0 * (SIGN(1.D0, slopes(ip-1)) + SIGN(1.D0,slopes(ip))) *  MIN(ABS(slopes(ip-1)), ABS(slopes(ip)))      
        pdot = -third* divV  + avg_slope * Dpp / (p(ip)**2)    
         
        ! compute upwind flux        
        IF(pdot .LT. 0.D0) THEN
                
            ! need p^3*f and p^4*f for computing fluxes, can get these from n and g 
            q3  = slopes(ip) - 3.D0
            q4  = slopes(ip) - 4.D0
            d3  = EXP(-q3*dlp)
            d4  = dp * d3           
                    
            IF(ip .LE. npbins) THEN 
                p3f = work_np(ip) * q3 / (1.D0 - d3) 
                p4f = work_gp(ip) * q4 / (1.D0 - d4) 
            ELSE
                ! assume n = g = 0 at highest momentum bin (i.e. no-inflow)
                p3f = 0.D0
                p4f = 0.D0
            END IF
                    
        ELSE
                
            ! need p^3*f and p^4*f for computing fluxes, can get these from n and g 
            q3  = slopes(ip-1) - 3.D0
            q4  = slopes(ip-1) - 4.D0
            d3  = EXP(-q3*dlp)
            d4  = dp * d3           
                          
            p3f = work_np(ip-1) * q3 * d3 / (1.D0 - d3) 
            p4f = work_gp(ip-1) * q4 * d4 / (1.D0 - d4) 
 
            !p3f = work_np(ip-1) * q3 / (1.D0 - d3) 
            !p4f = work_gp(ip-1) * q4 / (1.D0 - d4) 
  
        END IF
                
        ! compute bin edge fluxes
        Fn(ip) = pdot * p3f  
        Fg(ip) = pdot * p4f 
        
        
        !#################
        CYCLE
        !#################
        
        
        
        ! If radiative losses are off, then just cycle
        IF(.NOT. CREloss) CYCLE
        
        ! add non-zero radiative loss term to momentum flux
        IF(pdot_loss .NE. 0.D0) THEN

            ! since pdot_loss < 0, loss term is upwinded to higher momentum 

            IF(ip .LE. npbins) THEN 

                q3  = slopes(ip) - 3.D0
                q4  = slopes(ip) - 4.D0
                d3  = EXP(-q3*dlp)
                d4  = dp * d3           
            
                ! set p4f and p5f from n and g
                p4f = work_np(ip-1) * p(ip) * q3 / (1.D0 - d3)  ! need to check these formulae p4f and p5f calculations
                p5f = work_gp(ip-1) * p(ip) * q4 / (1.D0 - d4) 

            ! assume n = g = 0 at highest momentum bin (i.e. no-inflow)
            ELSE
            
                p4f = 0.D0
                p5f = 0.D0
            
            END IF

        END IF         

        ! add radiative loss contribution to bin edge fluxes
        Fn(ip) = Fn(ip) + pdot_loss * p4f 
        Fg(ip) = Fg(ip) + pdot_loss * p5f
                
    END DO

    !PRINT*,'############################################'
    !PRINT*,'Fn = ',Fn
    !PRINT*,'############################################'


END SUBROUTINE compute_momentum_fluxes


! This subroutine reconstructs the intra-bin slopes (denoted here by 'q') from the n and g moments by 
! solving for the root of the equation [psi(q) = 0],  where psi(q) = g / n*p - q3*(1-d^-q4) / q4*(1-d^-q3)     
! Instead of computing the exact root using a root finder, we approximate it using a cubic hermite spline interpolation of pre-computed values.
SUBROUTINE compute_slopes_approx()

    INTEGER :: i, ix
    REAL(8) :: gnp, idgnp1, q3
    REAL(8) :: xx, xx2, xx3, dx, x1, x2, y1, y2, yp1, yp2
    REAL(8) :: h00, h10, h01, h11


    idgnp1 = 1.D0 / (gnp_knots(2) - gnp_knots(1)) 

    ! loop over momentum bins and compute slope
    DO i= 1, npbins
    
        gnp = work_gp(i) / (work_np(i) * p(i))   
        
        !*******************************************************************************        
        ! obtain slope in ith bin via cubic spline interpolation of pre-computed samples
        !*******************************************************************************
             
        ! first apply protection for out of bounds gnp value
        gnp = MAX(gnp_knots(1), MIN(gnp, gnp_knots(3)))    
 
        ! find out which knot point interval contains this gnp
        ix = 1 + (gnp - gnp_knots(1)) * idgnp1
    
        dx = gnp_knots(ix+1) - gnp_knots(ix) 
        xx = (gnp -  gnp_knots(ix)) / dx
        xx2 = xx * xx
        xx3 = xx2 * xx
        
        ! compute cubic spline interpolation co-efficients
        h00 = 2.D0 * xx3 - 3.D0 * xx2 + 1.D0 
        h10 =  xx3 - 2.0 * xx2 + xx  
        h01 = -2.D0 * xx3 + 3.D0 * xx2
        h11 = xx3 - xx2
   
        ! compute the interpolant 
        slopes(i) = h00 * q_knots(ix)   + h01 * q_knots(ix+1) + &
                       (h10 * qp_knots(ix)  + h11 * qp_knots(ix+1)) * dx
        
                
    END DO
        
    ! set slope in boundary bins (extrapolation from interior)
    slopes(0) = 2.D0 * slopes(1) - slopes(2)
    slopes(npbins+1) = 2.D0 * slopes(npbins) - slopes(npbins-1)


END SUBROUTINE compute_slopes_approx


! This subroutine reconstructs the intra-bin slopes (denoted here by 'q') from the n and g moments by 
! solving for the root of the equation [psi(q) = 0],  where psi(q) = g / n*p - q3*(1-d^-q4) / q4*(1-d^-q3)     
! with q3:= q-3 and q4:= q-4 using the Newton-Raphson root finder. 
SUBROUTINE compute_slopes()

    INTEGER :: i
    REAL(8) :: gnp, q3
    
    ! loop over momentum bins
    DO i= 1, npbins
    
        gnp = work_gp(i) / (work_np(i) * p(i)) 

        ! obtain slope in ith bin
        CALL newton_raphson(gnp,slopes(i))    
        
        
    END DO
        
    ! set slope in boundary bins (extrapolation from interior bins)
    slopes(0) = 2.D0 * slopes(1) - slopes(2)
    slopes(npbins+1) = 2.D0 * slopes(npbins) - slopes(npbins-1)


END SUBROUTINE compute_slopes

    
! Netwon-Raphson root solver
SUBROUTINE newton_raphson(gnp, q)

    REAL(8), INTENT(IN) :: gnp
    REAL(8), INTENT(INOUT) :: q
    REAL(8) :: psi, psip, q3, q4, d3, d4, q0, den, num, denp, nump, iden
    LOGICAL :: converged
    INTEGER :: iters
    converged = .FALSE.

    !initial guess value is set to 4.5
    q0 = 4.5
    
    ! clear iteration counter
    iters = 0

    DO WHILE(.NOT. converged .AND. iters .LT. maxIter)
            
        q3 = q0 - 3.D0
        q4 = q0 - 4.D0
        d3 = EXP(-q3*dlp)
        d4 = dp * d3 

        num  = q3 * (1.D0 - d4)
        den  = q4 * (1.D0 - d3)
        nump = (1.D0 - d4) + dlp * q3 * d4 
        denp = (1.D0 - d3) + dlp * q4 * d3

        
        ! Need protection from den = zero (which happens if q->3 or q-> 4). 
        ! To prevent numerical over/under flow, replace with linearized taylor 
        ! expansion (i.e. linear in q-4 if q->4 or linear in q-3 if q->3) and/or impose cut-off.

        IF(ABS(q4) .LT. ZERO) THEN
        
            ! need to figure out how to do this part correctly...
        
        ELSE IF(ABS(q3) .LT. ZERO) THEN

        END IF          
        
        psi = gnp - num / den       
        
        psip = -nump / den + num * denp / (den**2)

        ! compute next value in the iteration
        q = q0 - psi/psip 

        ! apply protection to q (q bounded in [-15,15] )
        q = MIN(qmax, MAX(-qmax, q))
        
        ! check for convergence (tolerance is set to 10^-3)
        IF(ABS(q-q0) .LT. 10.d-3) converged = .TRUE.

        iters = iters + 1

        ! update initial guess
        q0 = q
       
    END DO


END SUBROUTINE newton_raphson


! This function computes the flux of particles in a momentum bin due to radiative losses
! via synchrotron emissions and inverse compton scattering with CMB and other background radiation fields.
FUNCTION radiative_losses(B_sqr) RESULT(pdot_loss)

    REAL(8), INTENT(IN) :: B_sqr
    REAL(8) :: pdot
    REAL(8) :: UB, itsc
        
    ! compute local magnetic field energy density 
    UB = 0.5D0 * B_sqr
    
    ! compute inverse of characteristic synchrotron lifetime in code units
    itsc = 1.D0 / tsc 

    ! compute total momentum flux from all radiative losses
    pdot_loss = - itsc * (UB + Ucmb + Urad)
    
        
END FUNCTION radiative_losses


! This subroutine pre-computes 3 knot points which can be used to construct piecewise cubic hermite splines passing 
! through the region between them:  (gnp1, q1),  (gnp2, q2) and (gnp3, q3)
SUBROUTINE set_spline_knot_points()

    INTEGER :: i    
    REAL(8) :: q(3)
    REAL(8) :: psi, psip, q3, q4, d3, d4, q0, den, num
    REAL(8) :: q1, q2, gnp1, gnp2


    ! set 3 knot points q values (uniformly spaced in the interval [-15,15])
    q_knots = (/ 15.D0, 0.D0, -15.D0/)
    
    DO i = 1, 3
    
        q0 = q_knots(i)
        q3 = q0 - 3.D0
        q4 = q0 - 4.D0
        d3 = EXP(-q3*dlp)
        d4 = dp * d3 

        num  = q3 * (1.D0 - d4)
        den  = q4 * (1.D0 - d3)
              
        ! compute gnp value at this knot point
        gnp_knots(i) = num / den
    
        ! now compute the derivative q'(gnp) at this knot point
        q0 = q_knots(i) - 1.D-8
        q3 = q0 - 3.D0
        q4 = q0 - 4.D0
        d3 = EXP(-q3*dlp)
        d4 = dp * d3 

        num  = q3 * (1.D0 - d4)
        den  = q4 * (1.D0 - d3)
    
        q1 = q0
        gnp1 = num / den
    
        q0 = q_knots(i) + 1.D-8
        q3 = q0 - 3.D0
        q4 = q0 - 4.D0
        d3 = EXP(-q3*dlp)
        d4 = dp * d3 

        num  = q3 * (1.D0 - d4)
        den  = q4 * (1.D0 - d3)
        
        q2 = q0
        gnp2 = num/den
    
        qp_knots(i) = (q2 - q1) / (gnp2 - gnp1) 
            
    END DO     
             
         
END SUBROUTINE set_spline_knot_points


! initial conditions for adiabatic expansion test (i.e. momentum space advection)
SUBROUTINE adiabatic_expansion_test_init(fp, np, gp, p_in, fp2_cc)

    REAL(8), INTENT(INOUT) :: fp(:), np(:), gp(:), p_in(0:npbins+1), fp2_cc(:)
    INTEGER :: i
    REAL(8) :: pLC, pHC, N0, q3, q4, d3, d4, p3, p4
    
    
    p_in = p
    
    pLC = 1.D3 !100.D0
    pHC = 1.D9 !300.D0
    N0  = 1.D5 !.D3
    
    DO i = 1, npbins
   
        ! set f    
        fp(i) = MAX(1.D-30, N0 * p(i)**(-2.D0) * (1.D0 - p(i) / pHC) * (1.D0 - pLC / p(i)) )    
    
        ! set fp2_cc
        fp2_cc(i) = fp(i) * (p(i)**2)
    
        q3 = q0 - 3.D0
        q4 = q0 - 4.D0
        d3 = EXP(-q3 * dlp)
        d4 = dp * d3 
        p3 = p(i)**3
        p4 = p3 * p(i)

        ! set n
        np(i) = fp(i) * p3 *(1.D0 - d3) / q3 

        ! set g
        gp(i) = fp(i) * p4 *(1.D0 - d4) / q4       

    END DO


    PRINT*,''
    PRINT*,'Initialization for adiabatic expansion test completed..'
    PRINT*,''



END SUBROUTINE adiabatic_expansion_test_init



! initial conditions for adiabatic expansion test (i.e. momentum space advection)
SUBROUTINE diffusion_test_init(fp, np, gp, p_in, fp2_cc)

    REAL(8), INTENT(INOUT) :: fp(:), np(:), gp(:), p_in(0:npbins+1), fp2_cc(:)
    INTEGER :: i
    REAL(8) :: q3, q4, d3, d4, p3, p4, pC, sig, f0
    
    
    p_in = p
    pC = 1.D4
    sig = 1.D3
    f0 = 1.D15 !1.D-85
 
    PRINT*,''
    
    DO i = 1, npbins
   
        ! set f    
        !fp(i) = MAX(1.D-30, f0 * EXP(-(p(i)-pC)/sig)**2)    
        fp(i) = MAX(1.D-30, f0 * EXP(-((p(i)-pC)/sig)**2))    
    
        PRINT*,'p, f(p) = ',p(i), fp(i)
    
        ! set fp2_cc
        fp2_cc(i) = fp(i) * (p(i)**2)
    
        q3 = q0 - 3.D0
        q4 = q0 - 4.D0
        d3 = EXP(-q3 * dlp)
        d4 = dp * d3 
        p3 = p(i)**3
        p4 = p3 * p(i)

        ! set n
        np(i) = fp(i) * p3 *(1.D0 - d3) / q3 

        ! set g
        gp(i) = fp(i) * p4 *(1.D0 - d4) / q4       

    END DO


    

    PRINT*,''
    PRINT*,'Initialization for diffusion expansion test completed..'
    PRINT*,''



END SUBROUTINE diffusion_test_init



END MODULE fp_solver_mod