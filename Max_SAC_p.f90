     Program Test

       Use MaxEnt_stoch_mod
       Use Files_mod

       Implicit Real (KIND=8) (A-G,O-Z)
       Implicit Integer (H-N)

       Real (Kind=8), Dimension(:), allocatable :: XQMC, XTAU, Alpha_tot, xom, A, &
            & XDATA,FDATA,ERROR, ARES, XQMC_ST
       Real (Kind=8), Dimension(:,:), allocatable :: XCOV, Xn, Cov_st
       Real (Kind=8), External :: XKER, Back_trans_Aom, F_Fit
       Character (64) :: File_root, File1

       open(unit=30,file='paramSAC',status='old',action='read', iostat=io_error)
       if (io_error.eq.0) then
          read(30,*) Ngamma, OM_st, OM_en, Ndis, NBins, NSweeps, Nwarm
          read(30,*) N_alpha, alpha_st, R
       else
          write(6,*) 'No file paramSAC! '
          stop
       endif
       close(30)
       open (unit=10,File="g_dat", status="unknown")
       read(10,*) ntau
       close(10)
       Allocate ( XCOV(NTAU,NTAU), XQMC(NTAU), XQMC_ST(NTAU), XTAU(NTAU) )


       L_cov = 0
       Xcov = 0.d0
       open (unit=10,File="g_dat", status="unknown")
       read(10,*) n
       do nt = 1,ntau
          read(10,*) Xtau(nt),Xqmc(nt), err
          xcov(nt,nt) = err*err
       enddo
       if (L_cov.eq.1) then
          do nt = 1,ntau
             do nt1 = 1,ntau
                read(10,*) xcov(nt,nt1)
             enddo
          enddo
       endif
       close(10)

       xqmc_st = xqmc
       dtau = Xtau(2) - Xtau(1)
       Beta = ntau*Dtau

       xmom1 = 1.0d0
       Allocate (Alpha_tot(N_alpha) )
       do nt = 1,N_alpha
          alpha_tot(nt) = alpha_st*(R**(nt-1))
       enddo
       Open(unit=50,File='Info',Status="unknown")
       write(50,*) 'First Moment, Beta  ', Xmom1, Beta
       write(50,*) '# Number of nodes   ', Isize
       close(50)
       

       Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER, Back_Trans_Aom, Beta, &
            & Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm)


       N_alpha_1 = N_alpha - 10
       File_root = "Aom_ps"
       File1 = File_i(File_root,N_alpha_1)
       Open(Unit=66,File=file1,status="unknown")
       Allocate (xom(Ndis), A(Ndis))
       do nw = 1,Ndis
          read(66,*) xom(nw), A(nw), x, x1, x2
       enddo
       close(66)
       Dom = xom(2) - xom(1)
       
       Open (Unit=70,File="data_out", status="unknown")
       do nt = 1,Ntau
          X = 0.d0
          do nw = 1, Ndis
             X = X + Xker(Xtau(nt),xom(nw), beta)*A(nw)
          enddo
          X = X *dom
          write(70,2005) xtau(nt), xqmc_st(nt), sqrt(xcov(nt,nt)), X
       enddo
2005   format(F14.7,2x,F14.7,2x,F14.7,2x,F14.7)
       Close(70)
     end Program Test
     
     
     
     Real (Kind=8) function XKER(tau,om, beta)
       
       Implicit None
       real (Kind=8) :: tau, om, pi, beta
       
       !pi = 3.1415927
       XKER = exp(-tau*om) / ( 1.d0 + exp(-Beta*om) ) ! /pi
       
     end function XKER
     
     Real (Kind=8) function Back_trans_Aom(Aom, om, beta)
       
       Implicit None
       real (Kind=8) :: Aom, om, beta
       
       Back_trans_Aom = Aom
       
     end function BACK_TRANS_AOM
