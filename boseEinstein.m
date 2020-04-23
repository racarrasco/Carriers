function[boseResults] = boseEinstein(x, a, b, theta)

%BoseEinstein according to Vina cardona
%boseResults = a - b.*(1 + 2./(exp(theta./x)-1));


boltzmann = 1.38064852e-23;
e = 1.60217663e-19;

%0.0862 meV/K
kmeV = boltzmann/e * 1000;

%Giving direct physical significance to Bose Einstein fit parameters
%Formula found in P. T. Webster, N. A. Riordan, S. Liu, E. H. Steenbergen,
%R. A. Synowicki, Y. -H. Zhang, and S. R. Johnson
%J. Appl. Phys. 118, 245706 (2015). Measurement of InAsSb bandgap energy
%and InAs/InAsSb band edge positions using spectroscopic ellipsometry and
%photluminescence spectroscopy
boseResults = a - b*kmeV*theta./(exp(theta./x)-1);