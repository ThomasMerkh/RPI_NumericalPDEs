%% advection examples

LF  = 1; % Lax-Friedrichs
LW  = 2; % Lax-Wendroff
FOU = 3; % first-order upwind

GP  = 1; % Gaussian pulse
D   = 2; % discontinuity


%% first we look at large CFL and smooth solution
% Lax-Friedrichs
c =   1.; 
tf =  .5;
CFL = .9;
advection1D( 401,c,CFL,tf,LF,GP )
pause

advection1D( 401,c,CFL,tf,LW,GP )
pause

advection1D( 401,c,CFL,tf,FOU,GP )
pause

%% large CFL and discontinuous solution
% Lax-Friedrichs
advection1D( 401,c,CFL,tf,LF,D )
pause

advection1D( 401,c,CFL,tf,LW,D )
pause

advection1D( 401,c,CFL,tf,FOU,D )