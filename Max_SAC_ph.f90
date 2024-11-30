     Program Test

       Use MaxEnt_stoch_mod
  
       Implicit Real (KIND=8) (A-G,O-Z)
       Implicit Integer (H-N)
       
       Real (Kind=8), Dimension(:), allocatable :: XQMC, XTAU, Alpha_tot
       Real (Kind=8), Dimension(:,:), allocatable :: XCOV, Xn
       Real (Kind=8), External :: XKER, Back_trans_Aom

       
       Open(Unit=50,file='info', status = 'unknown')

       open(unit=30,file='paramSAC_PH',status='old',action='read', iostat=io_error) 
       if (io_error.eq.0) then
          write(50,*) 'Reading in from paramSAC_PH'
          read(30,*)  Ngamma, OM_st, OM_en, Ndis, NBins, NSweeps, Nwarm
          read(30,*)  N_alpha, alpha_st, R
          read(30,*)  rel_err
          write(50,*) 'Frequency range: ', OM_st, ' to ', OM_en
          write(50,*) 'Tolerated relative error ',  rel_err
       else
          write(6,*) 'No file paramSAC_PH! ' 
          stop
       endif

       !X0: this possible offset is equal to unity for den den-den correlations
       !    and zero for the spin-spin correlations
       ! If < n (tau) n (0) > and not < n(tau) n(0) > - < n(tau)> < n(0)>
       ! is computed, then need this offset

       ! Differs from Max_SAC_p.f90
       ! IMPORTANT: Here assuming symmetry over beta/2
       ! this symmetry can be absent. May change to beta
       ! then the same as Max_SAC_p.f90, except kernal and Back_trans_Aom
       open (unit=10,File="g_dat", status="unknown") 
       read(10,*) beta, ntau, X0
       Write(50,*) 'Beta is:  ', beta
       Write(50,*) 'Using the data from 0 to Beta/2 '
!       ntau = ntau/2
       Allocate ( XCOV(NTAU,NTAU), XQMC(NTAU),XTAU(NTAU) )
       XCOV  = 0.d0
       Do nt = 1,NTAU
          read(10,*) xtau(nt), xqmc(nt), err 
          xqmc(nt) = Xqmc(nt) - X0
          xcov(nt,nt) = err*err
       Enddo
       close(10)

       pi = acos(-1.d0)

       !First Moment xmom1
       !Differs from xmom1 = 1.0d0 in *p.f90
       xmom1 = pi * xqmc(1) 
       Allocate (Alpha_tot(N_alpha) )
       do nt = 1,N_alpha
          alpha_tot(nt) = alpha_st*(R**(nt-1))
       enddo

       close (50)
       Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER, Back_Trans_Aom, Beta, &
            &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm) 


     end Program Test
     

     ! Differs from Max_SAC_p.f90, in which
     ! XKER = exp(-tau*om) / ( 1.d0 + exp(-Beta*om) )

     ! IMPORTANT: here exp(-tau*om) + exp(-( beta - tau )*om
     ! accounts for the symmetry over beta/2 for two-particle correlation
     Real (Kind=8) function XKER(tau,om, beta)

       Implicit None
       real (Kind=8) :: tau, om, pi, beta

       pi = 3.1415927

!       XKER = (exp(-tau*om) + exp(-( beta - tau )*om ) ) / ( pi*(1.d0 + exp( - beta * om ) ) )

       !Switch to the same as Max_SAC_p.f90, which do not consider symmetry over beta/2
       XKER = exp(-tau*om) / ( pi*(1.d0 + exp( - beta * om ) ) )

     end function XKER

     ! Differs from Max_SAC_p.f90, in which
     ! Back_trans_Aom = Aom
     ! See Notes.pdf, which says 
     ! If you use this wrapper the output in the files Aom_ps_* corresponds to chi''/(1-exp(-beta*omega))
     ! which seems not correct, From Eq.(5) and below, should just be chi''
     Real (Kind=8) function Back_trans_Aom(Aom, om, beta)

       Implicit None
       real (Kind=8) ::  Aom, om, beta

       Back_trans_Aom =  Aom* ( exp(beta*om/2) - exp(-beta*om/2) )/&
            &                 ( exp(beta*om/2) + exp(-beta*om/2) )

     end function BACK_TRANS_AOM

     
