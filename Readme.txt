
To compile:

make

But to make sure that it is well compiled, I usually remove the .o
files first then do make.

The make will produce two executables, Max_SAC_ph.out and
Max_SAC_p.out. 

The ph one is the ``particle-hole'' and is to be used
for density-density and spin-spin correlation functions. 

If the correlation function is spin-spin, the first line in the input
file, must be called g_dat, has three entries: beta, number of data
lines to follow, and then 0. After that the data follow in three
columns: tau, value, error. tau should start at 0 and go to
beta-dtau. 

If the correlation function is <n(tau)n(tau')>, then instead of
the 0 on the first line we enter 1. This tells the code that it must
shift the data by <n(beta/2)>^2. 

The p code is for the single particle Green fct. The input file is
g_dat, the first line only has the number of time slices to
follow. Then the data follow in the same way as above.

The output files are the Aom_ps_#. I usually look at a few before I
decide which one to keep.

Below is Fakher's email with a few comments by me.

*****************************************************************

Hi George,

here are the files with analysis. In principle, you can now use things
as a black box!  Well I hope so.  Just tell me if it doesn't work.

I have slightly changed the wrapper Max_SAC_ph.f90 in the directory
FAKHER.  Note that "_ph" stands for particle-hole so that you should
use this wrapper for the analysis of the den-den or spin-spin
correlations.

[GGB comment: I checked with Fakher. This Max_SAC_ph.f90 is for both
Fermions and Bosons. No changes needed to use it for one or the
other.]

I have done two test tuns (those are rather short runs) which you will
find in the directories denN6L80 and spinN6L80.

The QMC data is in the file g_dat. Here I have added a first line with
three numbers beta, # of tau points you are providing in the range 0
to beta, offset.

The code fits the data  g(tau) - offset

The offset is equal to unity for den den-den correlations and zero for
the spin-spin correlations.  As far as I understand you have measured
< n (tau) n (0) > and not < n(tau) n(0) > - < n(tau)> < n(0)>
(correct?)  and this is why I had to insert the offset.

[I checked with Fakher. I have a density-density file where in fact I
calculate <n(tau)n(0)>-<n(tau)><n(0)>, then I should just use
offset=0].

The parameter file is the usual one and is now called paramSAC_PH.

You can get smoother results  by enhancing the number of bins.

Hope that this helps, Greetings, Fakher.


************THE HOW TO FILE FROM THE NEWMAXENT DIR**************

I think it is relevant here too.

****************************************************************
****************************************************************
****************************************************************

The following is for the new boson and fermion codes:

FROM FAKHER:
I have added a wrapper for  fermions: MaxEnt_p.f90 and renamed 
the old wrapper for bosonic excitations to MaxEnt_ph.f90. 
(ph stands for particle-hole excitations which are bosonic).
The input file is the same as before. Since the G(tau) is not 
symmetric around beta/2 any more you have to provide the whole
G(tau) for tau from [ 0, beta[.

****************************************************************
****************************************************************
****************************************************************

Instructions on the use of Fakher Assaad's stochastic MaxEnt code.

To compile: make

It will compile both the Fermionic and bosonic guys. 

Uses the ifort compiler, so the right path should be in the makefile.

THE BOSONIC EXECUTABLE is called Max_SAC_ph.out and should be run in
the same directory as the input file. Just type Max_SAC_ph in the same
directory where g_cov and paramSAC are sitting.

THE FERMIONIC EXECUTABLE is called Max_SAC_p.out and should be run in
the same directory as the input file. Just type Max_SAC_ph in the same
directory where g_dat and paramSAC are sitting.


THE BOSONIC INPUT FILE is called g_cov. The first entry in the file is
the number of data lines to follow. Then right after, you have the
data in the form

     tau   corr-fct   error

and since these correlation functions in tau are periodic, you only
give the program half the period. So, if your beta=10, the data go
from tau=0 to tau=5. Clearly, corr-fct is the time correlation of
interest. 

THE FERMIONIC INPUT FILE is called g_dat. The first entry in the file
is the number of data lines to follow. Then right after, you have the
data in the form

     tau   corr-fct   error

and since these correlation functions in tau are NOT periodic, you
must give the program the whole corr function in the interval
[0,beta[. So, if your beta=10, the data go from tau=0 to
tau=10-dtau. NOTE: THE DATA MUST BE POSITIVE EVEN IF THE CORRELATION
FUNCTION IS NEGATIVE IN THE FERMION CASE.

IMPORTANT: FOR BOTH THE FERMIONS AND BOSONS, tau IS NOT THE INTEGER
TIME SLICE, IT IS THE IMAGINARY TIME.

BOSONIC PARAMETER FILE is called paramSAC. I do not understand most of
the parameters in it. But, if you want to change the range of omega
you want, you have two parameters OM_st, OM_en correspeonding to
starting and ending values of omega. I guess we always start at
omega=0 but you can easily change the ending omega. A typical file
looks like this:

250,  0.0,   20.0,  1000,  200,    100,  20
28      1.0     1.2
0.1
Ngamma, OM_st, OM_en, NDis, NBins, NSweeps, NWarm
N_alpha, alpha_st, R
rel_err


FERMIONIC PARAMETER FILE is called paramSAC. Fakher recommends the
values:

250,  -10.0  10.0,  1000,  50,    100,  20
20      0.1     1.2
0.1
Ngamma, OM_st, OM_en, NDis, NBins, NSweeps, NWarm
N_alpha, alpha_st, R
rel_err 

Note that OM_st=-10 now not 0 as for the bosons since the results is
not necessarily symmetric in gamma.

Also Fakher writes: ``With the old parameter N_alpha = 28 alpha_st=1.0
and R= 1.2 I believe that we were over fitting the data.  In
Wuerzburg, we use the following N_alpha= 20 alpha_st= 0.1 and R= 1.2.
The bottom line is that we were averaging over spectral functions whit
chisq << # of data points. With the new choice we allow for large
value of chisq.''



ENERGIES FILE: energies. This is essentially the chisq as a function
of 1/T.


OUTPUT FILES are called Aom_ps_NUMBER. Spectral functions. Each
simulation at a given temperature yields a spectral function with a
given chisq, or energy. If the temperature is too low you are
overfitting and of the temperature is to high you are
underfitting. The files Aom_ps_$ average of the spectral functions
from T=$ to the lowest temperature.  In the present form of the
stochastic MaxEnt, there is no clear cut way of deciding below which
temperature you should start averaging.  However, if you take a look
at the energy values, then you will see that they first decrease
rapidly before saturating.  K. Beach proposes to start averaging below
this crossover scale. 



FORT.44: Not important.





