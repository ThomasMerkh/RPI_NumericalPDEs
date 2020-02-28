nu = 1.0;
tf = 1.0;

%%%%%
%% TZ with Crank Nicolson (dt~dx)
%%%%%
Nr = 5;
Ns = 40;
Heat2DAnnulus( Nr,Ns,nu,tf,1,1,.5 );


Nr = 10;
Ns = 80;
Heat2DAnnulus( Nr,Ns,nu,tf,1,1,.5 );

Nr = 20;
Ns = 160;
Heat2DAnnulus( Nr,Ns,nu,tf,1,1,.5 );


%%%%%
%% boundary driven CN (dt~dx)
%%%%%
Nr = 20;
Ns = 160;
Heat2DAnnulus( Nr,Ns,nu,tf,3,1,.5 );

%%%%%
%% boundary driven CN (dt~dx^2)
%%%%%
Nr = 20;
Ns = 160;
Heat2DAnnulus( Nr,Ns,nu,tf,3,1,.5+eps );

%%%%%
%% tophat CN (dt~dx^2)
% Fun Example 1
%%%%%
Nr = 20;
Ns = 160;
Heat2DAnnulus( Nr,Ns,nu,tf,2,1,.5+eps );

%%%%%
%% tophat CN (dt~dx)
%%%%%
Nr = 20;
Ns = 160;
Heat2DAnnulus( Nr,Ns,nu,tf,2,1,.5 );

%%%%%
%% tophat BE (dt~dx)
% Fun Example 2: Line 182 comment out 'fixed axis' for adjusting axis.
%%%%%
Nr = 20;
Ns = 160;
Heat2DAnnulus( Nr,Ns,nu,tf,2,1,0.0 );