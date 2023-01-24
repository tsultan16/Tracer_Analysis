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

REAL(8), ALLOCATABLE :: p(:), q(:), &
                        delp(:), delq(:)

REAL(8) :: dx, dt, dlp, dp, pmin, pmax, dgp

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
    pmax  = 1.D15

    dlp = (LOG(pmax)-LOG(pmin))/npbins  ! uniform logarithmic spacing 
    dp  = EXP(dlp)  

    ! allocate memory for all work arrays
    CALL create_workpool()

  
    PRINT*,''
    PRINT*,"Fokker-Planck Solver Parameters:"
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
          
    ALLOCATE(p(0:npbins+1), delp(0:npbins+1), delq(0:npbins+1),q(0:npbins+1))        
  
    p       = 0.D0
    q       = 0.D0
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

    DEALLOCATE(p,q,delp,delq)
    

END SUBROUTINE destroy_workpool



! top-level interface for performing a full CGMV update 
SUBROUTINE solve(dt_in, fp2_cc, Bsqr, divV, update_cmplt)

    REAL(8), INTENT(IN)  :: dt_in, Bsqr, divV
    REAL(8), INTENT(INOUT)  :: fp2_cc(:)
    LOGICAL, INTENT(OUT) :: update_cmplt
    
    
    dt = dt_in
    
    ! update particle distribution for adiabatic changes and losses

    CALL chang_cooper(fp2_cc, Bsqr, divV)
    
    CALL diagnostics()
    
    
    update_cmplt = .TRUE.    


END SUBROUTINE solve



SUBROUTINE diagnostics()

    INTEGER :: i
    REAL(8) :: del_np_sum, del_gp_sum


    !PRINT*,''
    !WRITE(*,FMT='("np_sum = ",E11.4,", fractional del_np_sum",E11.4)') np_sum, del_np_sum
    !WRITE(*,FMT='("gp_sum = ",E11.4,", fractional del_gp_sum",E11.4)') gp_sum, del_gp_sum
    !PRINT*,''

END SUBROUTINE diagnostics


!*****************************
!    Chang Cooper Scheme     
!*****************************


! Updates momentum convection part using Chang-Cooper Scheme
SUBROUTINE chang_cooper(fp2, Bsqr, divV)

    REAL(8), INTENT(INOUT) :: fp2(:)
    REAL(8), INTENT(IN) :: Bsqr, divV
    REAL(8) :: pdot, pdot_max, pdot_loss, subdt
    INTEGER :: i, nsub, n, ip, j
    REAL(8) :: H(0:npbins+1), Wp(0:npbins+1), Wm(0:npbins+1)
    REAL(8) :: Dpp(1:npbins), Qp(1:npbins)
    REAL(8) :: A(0:npbins+1), B(0:npbins+1), C(0:npbins+1), r(0:npbins+1)
    REAL(8) :: gam(0:npbins+1)
    REAL(8) :: wi, expwi, bet
    
    ! compute moomentum velocity (pdot) from radiative losses
    pdot_loss = 0.5D0 * Bqsr ! self%radiative_losses(Bsqr)    

    


    ! compute required time step size with CFL number = 0.8 
    subdt =       !  dt / DBLE(nsub) !COUR_p * pdot_max / tls%dlp  
        
    ! compute number of subcycle steps
    nsub = 2000 !20 !60 !120      
    
    
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
            H(ip) = -2.D0 * (Dpp(ip) / q(ip)) + (third * divV + pdot_loss) * q(ip)  ! adiabatic + diffusion + radiative loss
            
                    !-2.D0 * Dpp(ip) / q(ip) + third * divV * q(ip)  ! only adiabatic and diffusion terms, not including radiative losses
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



! initial conditions for adiabatic expansion test (i.e. momentum space advection)
SUBROUTINE adiabatic_expansion_test_init(p_in, fp2_cc)

    REAL(8), INTENT(INOUT) :: p_in(0:npbins+1), fp2_cc(:)
    INTEGER :: i
    REAL(8) :: pLC, pHC, N0
    
    
    p_in = p
    
    pLC = 1.D6 !1.D3 !100.D0
    pHC = 1.D8 !300.D0
    N0  = 1.D5 !.D3
    
    DO i = 1, npbins
   
        ! set fp2_cc  
        fp2_cc(i) = MAX(1.D-30, N0 * p(i)**(-4.D0) * (1.D0 - p(i) / pHC) * (1.D0 - pLC / p(i)) * p(i)**(2.D0) )  

    END DO


    PRINT*,''
    PRINT*,'Initialization for adiabatic expansion test completed..'
    PRINT*,''



END SUBROUTINE adiabatic_expansion_test_init



! initial conditions for adiabatic expansion test (i.e. momentum space advection)
SUBROUTINE diffusion_test_init(p_in, fp2_cc)

    REAL(8), INTENT(INOUT) :: p_in(0:npbins+1), fp2_cc(:)
    INTEGER :: i
    REAL(8) :: pC, sig, N0
    
    
    p_in = p
    pC = 1.D4
    sig = 1.D5
    N0 = 1.D5 !1.D-85
 
    PRINT*,''
    
    DO i = 1, npbins
   
        ! set fp2_cc
        fp2_cc(i) = MAX(1.D-30, N0 * EXP(-((p(i)-pC)/sig)**2))  * (p(i)**2)
    
    END DO
    

    PRINT*,''
    PRINT*,'Initialization for diffusion expansion test completed..'
    PRINT*,''



END SUBROUTINE diffusion_test_init



END MODULE fp_solver_mod