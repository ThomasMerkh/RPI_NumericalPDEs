function [x,t,uST] = wave1D_4th( Nx,Nt,tf,c,iOption )
  xa = 0;
  xb = 1;
  
  dx = (xb-xa)/Nx;
  dt = tf/Nt;
  sigma = c*dt/dx;
  
  % set the number of ghost cells
  ng  = 2; 
  NXT = Nx+1+2*ng;
  
  x = linspace(xa-ng*dx,xb+ng*dx,NXT);
  t = linspace(0,tf,Nt+1);
  
  %% setup space time data array for later plotting
  uST = zeros( NXT,Nt+1 );
  
  %% setup u,un, unm1
  u    = zeros(NXT,1);
  un   = zeros(NXT,1);
  unm1 = zeros(NXT,1);
  
  %% set initial conditions over domain interior
  for j = ng+1:NXT-ng
    fm2 = f(x(j-2),iOption);
    fm1 = f(x(j-1),iOption);
    f0  = f(x(j)  ,iOption);
    fp1 = f(x(j+1),iOption);
    fp2 = f(x(j+2),iOption);
    
    gm2 = g(x(j-2),iOption);
    gm1 = g(x(j-1),iOption);
    g0  = g(x(j)  ,iOption);
    gp1 = g(x(j+1),iOption);
    gp2 = g(x(j+2),iOption);
    
    unm1(j) = f0;
    un(j)   = ...
       f0...
      +dt*g0...
      +c^2*dt^2/(2.*12.*dx^2)*(-fp2+16.*fp1-30.*f0+16.*fm1-fm2)...
      +c^2*dt^3/(6.*12.*dx^2)*(-gp2+16.*gp1-30.*g0+16.*gm1-gm2)...
      +c^4*dt^4/(24.*dx^4)   *(fp2 -4.*fp1 +6.*f0 -4.*fm1+fm2);
  end
  
  %% set BCs
  unm1 = setBCs( unm1,t(1),dx,c,Nx,ng,iOption );
  un   = setBCs( un,  t(2),dx,c,Nx,ng,iOption );
    
  %% add first 2 solns to plotting space-time array
  uST(:,1) = unm1;
  uST(:,2) = un;
  
  %% plot solution at t=dt
  figure(1)
  plot( x,un,'k-' );
  xlabel( 'x' );
  ylabel( 'u' );
  axis( [xa,xb,-1.1,1.1] );
  pause
 
  %% time stepping loop
  for n = 3:Nt+1
    %% loop over interior
    for j = ng+1:NXT-ng
      u(j) = 2.*un(j)-unm1(j)...
        +sigma^2*(un(j+1)-2.*un(j)+un(j-1))...
        +sigma^2/12.*(sigma^2-1.)*(un(j+2)-4.*un(j+1)+6.*un(j)-4.*un(j-1)+un(j-2));
    end
    
    %% set BCs
    u = setBCs( u,t(n),dx,c,Nx,ng,iOption );

    %% plot to show movie
    plot( x,u,'k-' );
    xlabel( 'x' );
    ylabel( 'u' );
    axis( [xa,xb,-1.1,1.1] );
    pause( 0.01 );
    %pause
    
    %% update old solutions and plotting space-time array
    unm1 = un;
    un   = u;
    uST(:,n) = u;
  end
  
  inds = ng+1:NXT-ng;
  x = x(inds);
  uST = uST(inds,:);
  
  return
end

function u = setBCs( u,t,dx,c,Nx,ng,iOption )
    %% compatibility BCs    
    %% set BCs at left
    ii = ng+1;
    u(ii) = alpha(t,iOption);
    u(ii-1) = ...
       dx^2/c^2*alpha_tt(t,iOption)...
      +dx^4/(12.*c^4)*alpha_tttt(t,iOption)...
      +2.*u(ii)...
      -u(ii+1);
    u(ii-2) = ...
      dx^4/c^4*alpha_tttt(t,iOption)...
      +4.*u(ii-1)...
      -6.*u(ii)...
      +4.*u(ii+1)...
      -u(ii+2);
    
    %% set BCs at right
    ii = Nx+1+ng;
    u(ii) = beta(t,iOption);
    u(ii+1) = ...
       dx^2/c^2*beta_tt(t,iOption)...
      +dx^4/(12.*c^4)*beta_tttt(t,iOption)...
      +2.*u(ii)...
      -u(ii-1);
    u(ii+2) = ...
      dx^4/c^4*beta_tttt(t,iOption)...
      +4.*u(ii+1)...
      -6.*u(ii)...
      +4.*u(ii-1)...
      -u(ii-2);
    
return
end

function z = f(x,iOption)
  if( iOption == 1 )
    z = exp(-200.*(x-0.3)^2);
  elseif( iOption == 2 )
    z = 0.;
  elseif( iOption == 3 )
    if( (x-0.3)^2 < .15^2 )
      z = 1.;
    else
      z = 0.;
    end
  else
    z = sin(5.*x);
  end
  return
end

function z = g(x,iOption)
  if( iOption == 1 )
    z = 0.;
  elseif( iOption == 2 )
    z = 0.;
  elseif( iOption == 3 )
    z = 0.;
  else
    z =-5.*.9*cos(5.*x);
  end
  return
end

function z = alpha(t,iOption)
  if( iOption == 1 )
    z = 0.;
  elseif( iOption == 2 )
    z = 0.5*(sin(10*t));
  elseif( iOption == 3 )
    z = 0.;
  else
    z = sin(-5.*.9*t);
  end
  return
end

function z = alpha_tt(t,iOption)
  if( iOption == 1 )
    z = 0.;
  elseif( iOption == 2 )
    z = -0.5*10^2*(sin(10*t));
  elseif( iOption == 3 )
    z = 0.;
  else
    z = -(-5.*.9)^2*sin(-5.*.9*t);
  end
  return
end

function z = alpha_tttt(t,iOption)
  if( iOption == 1 )
    z = 0.;
  elseif( iOption == 2 )
    z = 0.5*10^4*(sin(10*t));
  elseif( iOption == 3 )
    z = 0.;
  else
    z = (-5.*.9)^4*sin(-5.*.9*t);
  end
  return
end

function z = beta(t,iOption)
  if( iOption == 1 )
    z = 0.;
  elseif( iOption == 2 )
    z = 0.;
  elseif( iOption == 3 )
    z = 0.;
  else
    z = sin(5.*(1.-.9*t));
  end
  return
end

function z = beta_tt(t,iOption)
  if( iOption == 1 )
    z = 0.;
  elseif( iOption == 2 )
    z = 0.;
  elseif( iOption == 3 )
    z = 0.;
  else
    z = -(-5.*.9)^2*sin(5.*(1.-.9*t));
  end
  return
end

function z = beta_tttt(t,iOption)
  if( iOption == 1 )
    z = 0.;
  elseif( iOption == 2 )
    z = 0.;
  elseif( iOption == 3 )
    z = 0.;
  else
    z = (-5.*.9)^4*sin(5.*(1.-.9*t));
  end
  return
end
    