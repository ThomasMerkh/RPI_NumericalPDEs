function [x,y,u,err] = twoWave( Nr,Ns,c,sigma0,tf,iOption,iTZ )
%%
% fun example : twoWave( 160,160,1,0.9,2,2,1 );
% Nr,Ns: number cells
% c: input parameter
% sigma0: estimated CFL number
% tf: final time
% iOption: selects case (1: TZ, 2: topHat)

% Note: r: radius
 % s: theta
 % t: time

 %% set some parameters
 ra = 1.; % inner radius
 rb = 2.; % outter radius
sa = 0.;
 sb = 0.5*pi;
 dr = (rb-ra)/(Nr);
 ds = (sb-sa)/(Ns);

 %% estimate time step
 circ = 0.5*pi*ra;
 h = min( dr,circ/Ns );
 dtGuess = sigma0*h/(2*c);
 nStep = ceil(tf/dtGuess);
 dt = tf/nStep;

 NRT = Nr+1+2; %% 2 ghost in radial
 NST = Ns+1; %% no ghost in theta
 %% set up polar grid
 r = linspace( ra-dr,rb+dr,NRT );
 s = linspace( sa,sb,NST );

 Ntot = NRT*NST;

 f = zeros(NRT,NST);
 g = zeros(NRT,NST);
 x = zeros(NRT,NST);
 y = zeros(NRT,NST);
 %% set ICs and physical grid
 t = 0;
 for k = 1:NST
 for j = 1:NRT
 x(j,k) = r(j)*cos(s(k));
 y(j,k) = r(j)*sin(s(k));
 f(j,k) = getICf( r(j), s(k), t, iOption, iTZ );
 g(j,k) = getICg( r(j), s(k), t, iOption, iTZ );
 end
 end

 f = setBCs( f,r,s,0,iOption,iTZ );
 utt = getdu( f,r,s,0,c,iOption,iTZ );

 unm1 = f;
 un = f+dt*g+dt^2/2*utt;
 un = setBCs( un,r,s,dt,iOption,iTZ );
 u = 0*un;


 mesh( x,y,un );
 shading interp
 xlabel( 'x' );
 ylabel( 'y' );
 pause

 %% time-stepping loop
 told = dt;
 for n = 1:nStep
 tnew = told+dt;

 utt = getdu( un,r,s,told,c,iOption,iTZ );

 u = 2.*un-unm1+dt^2*utt;
 u = setBCs( u,r,s,tnew,iOption,iTZ );
 unm1 = un;
un = u;
told = tnew;

if( iOption == 1 )
ue = 0*u;
for k = 1:NST
for j = 1:NRT
ue(j,k) = getU( r(j),s(k),tnew,iTZ );
end
end
surf( x,y,u-ue );
%surf( x,y,u );
xlabel( 'x' );
ylabel( 'y' );
shading interp
pause( 0.01 );
%pause
else
surf( x,y,u );
xlabel( 'x' );
ylabel( 'y' );
shading interp
view([0,90]);
 %pause
 pause( 0.01 );
 end

 told = tnew;
 end

 err = 0;
 if( iOption == 1 )
 ue = 0*u;
 for k = 1:NST
 for j = 1:NRT
 ue(j,k) = getU( r(j),s(k),tnew,iTZ );
 end
 end
 err = max(max(abs(u-ue)));
 fprintf( 'max error: %e\n', err );
 end

 indr = 1:NRT;
 inds = 2:NST-1;
 u = u(indr,inds);
 x = x(indr,inds);
 y = y(indr,inds);

 return
 end

 function up = getdu( u,r,s,t,c,iOption,iTZ )

 [NRT,NST] = size(u);
 up = zeros(NRT,NST);
 dr = r(2)-r(1);
 ds = s(2)-s(1);
 for k = 2:NST-1
 for j = 2:NRT-1
 rjp = 0.5*(r(j+1)+r(j));
 rjm = 0.5*(r(j)+r(j-1));
 up(j,k) = c^2*(...
 u(j,k)*(-1./(r(j)*dr^2)*(rjp+rjm)-2./(r(j)^2*ds^2))+...
 u(j-1,k)/(r(j)*dr^2)*rjm+...
 u(j+1,k)/(r(j)*dr^2)*rjp+...
 u(j,k-1)/(r(j)^2*ds^2)+...
 u(j,k+1)/(r(j)^2*ds^2));


 %urr = (u(j+1,k)-2.*u(j,k)+u(j-1,k))/(dr^2);
 %ur = (u(j+1,k)-u(j-1,k))/(2*dr);
 %uss = (u(j,k+1)-2.*u(j,k)+u(j,k-1))/(ds^2);
 %up(j,k) = c^2*(urr+1/r(j)*ur+1/r(j)^2*uss);

 %% add forcing if doing TZ
 if( iOption == 1 )
 up(j,k) = up(j,k)+getQ( r(j),s(k),t,c,iTZ );
 end
 end
 end

 return
 end

 function u = setBCs( u,r,s,t,iOption,iTZ )

 [NRT,NST] = size(u);
 dr = r(2)-r(1);
 ds = s(2)-s(1);

 %% Dirichlet BCs
 for j = 1:NRT
 u(j,1) = alphaL( r(j),s(1),t,iOption,iTZ );
 u(j,NST) = alphaR( r(j),s(NST),t,iOption,iTZ );
 end

 %% Neumann BCs
 for k = 2:NST-1
 u(1,k) = -2.*dr*betaL( r(2),s(k),t,iOption,iTZ )+u(3,k);
 u(NRT,k) = 2.*dr*betaR( r(NRT-1),s(k),t,iOption,iTZ )+u(NRT-2,k);
 end

 return
 end

 function z = getICf( r,s,t,iOption,iTZ )
 %% function to compute initial condition
 if( iOption == 1 ) %% TZ
 z = getU( r,s,t,iTZ );
 
 elseif( iOption == 2 ) %% topHat
    x = r*cos(s);
    y = r*sin(s);
    x0 = 1.;
    y0 = 1.;
    rad2 = (x-x0)^2+(y-y0)^2;
    if( rad2 < .25^2 )
      z = 1.e5;
    else
      z = 0.;
    end
 else %% problem from HW
 x = r*cos(s);
 y = r*sin(s);
 z = exp(-100*((r-1.5)^2+(s-pi/4.)^2));
 %z = exp(-100*((x-1.1)^2+(y-1.1)^2));
 end

 return
 end

 function z = getICg( r,s,t,iOption,iTZ )
 %% function to compute initial condition
 if( iOption == 1 ) %% TZ
 z = getU_t( r,s,t,iTZ );
 else %% problem from HW
 z = 0;
 end
 return
 end

 function a = alphaL( r,s, t, iOption,iTZ )
 %% function to BC at left
 if( iOption == 1 )
 a = getU( r,s,t,iTZ );
 else
 a = 0.;
 end
 return
 end

 function a = alphaR( r,s,t,iOption,iTZ )
 %% function to BC at right
 if( iOption == 1 )
 a = getU( r,s,t,iTZ );
 else
 a = 0;
 end
 return
 end

 function b = betaL( r,s,t,iOption,iTZ )
 %% function to BC at bottom
 if( iOption == 1 )
 b = getU_r( r,s,t,iTZ );
 else
 b = 0.;
 end
 return
 end

 function b = betaR( r,s,t,iOption,iTZ )
 %% function to BC at top
 if( iOption == 1 )
 b = getU_r( r,s,t,iTZ );
 else
 b = 0.;
 end
 return
 end

 function z = getU( r,s,t,iTZ )
 %% function to compute exact TZ solution
 if( iTZ == 1 )
 z = sin(3*r)*cos(2*s+1.)*cos(t);
 elseif( iTZ == 2 )
 z = s^2+s+r^2+r+t^2+t;
 else
 z = 1;
 end
 return
 end

 function z = getU_t( r,s,t,iTZ )
 %% function to compute exact TZ solution (t derivative)
 if( iTZ == 1 )
 z = -sin(3*r)*cos(2*s+1.)*sin(t);
 elseif( iTZ == 2 )
 z = 2*t+t;
 else
 z = 0.;
 end
 return
 end

 function z = getU_tt( r,s,t,iTZ )
 %% function to compute exact TZ solution (t derivative)
 if( iTZ == 1 )
 z = -sin(3*r)*cos(2*s+1.)*cos(t);
 elseif( iTZ == 2 )
 z = 2;
 else
 z = 0.;
 end
 return
 end

 function z = getU_r( r,s,t,iTZ )
 %% function to compute exact TZ solution (r derivative)
 if( iTZ == 1 )
 z = 3*cos(3*r)*cos(2*s+1.)*cos(t);
 elseif( iTZ == 2 )
 z = 2*r+1;
 else
 z = 0.;
 end
 return
 end

 function z = getU_rr( r,s,t,iTZ )
 %% function to compute exact TZ solution (rr derivative)
 if( iTZ == 1 )
 z = -9*sin(3*r)*cos(2*s+1.)*cos(t);
 elseif( iTZ == 2 )
 z = 2;
 else
 z = 0.;
 end
 return
 end

 function z = getU_s( r,s,t,iTZ )
 %% function to compute exact TZ solution (s derivative)
 if( iTZ == 1 )
 z = -2*sin(3*r)*sin(2*s+1.)*cos(t);
 elseif( iTZ == 2 )
z = 2*s+1;
 else
 z = 0.;
 end
 return
 end

 function z = getU_ss( r, s, t,iTZ )
 %% function to compute exact TZ solution (ss derivative)
 if( iTZ == 1 )
 z = -4*sin(3*r)*cos(2*s+1.)*cos(t);
 elseif( iTZ == 2 )
 z = 2;
 else
 z = 0.;
 end
 return
 end

 function Q = getQ( r,s,t,c,iTZ )
 %% function to compute TZ forcing
 ut = getU_t( r,s,t,iTZ );
 utt = getU_tt( r,s,t,iTZ );
 ur = getU_r( r,s,t,iTZ );
 urr = getU_rr( r,s,t,iTZ );
 us = getU_s( r,s,t,iTZ );
 uss = getU_ss( r,s,t,iTZ );

 Q = utt-c^2*(1./r*(ur+r*urr)+1./(r^2)*uss);
 return
 end
