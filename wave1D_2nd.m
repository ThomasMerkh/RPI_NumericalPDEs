function [x,t,uST] = wave1D_2nd( Nx,Nt,tf,c,iOption )
  xa = 0;
  xb = 1;
  
  dx = (xb-xa)/Nx;
  dt = tf/Nt;
  sigma = c*dt/dx;
  
  x = linspace(xa,xb,Nx+1);
  t = linspace(0,tf,Nt+1);
  
  %% setup space time data array for later plotting
  uST = zeros( Nx+1,Nt+1);
  
  %% setup u,un, unm1
  u    = zeros(Nx+1,1);
  un   = zeros(Nx+1,1);
  unm1 = zeros(Nx+1,1);
  
  %% set initial conditions
  for j = 2:Nx
    unm1(j) = f(x(j),iOption);
    un(j)   = f(x(j),iOption)...
      +dt*g(x(j),iOption)...
      +0.5*sigma^2*(f(x(j+1),iOption)-2.*f(x(j),iOption)+f(x(j-1),iOption));
  end
  
  %% set BCs
  unm1(1)    = alpha(t(1),iOption );
  unm1(Nx+1) = beta(t(1),iOption );
  un(1)      = alpha(t(2),iOption );
  un(Nx+1)   = beta(t(2),iOption );
  
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
    for j = 2:Nx
      u(j) = 2.*un(j)-unm1(j)+sigma^2*(un(j+1)-2.*un(j)+un(j-1));
    end
    
    %% set BCs
    u(1)    = alpha(t(n),iOption);
    u(Nx+1) = beta(t(n),iOption);
    
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
    