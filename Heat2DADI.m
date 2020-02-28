function [x,y,u,ue] = Heat2DADI( Nx,Ny,nStep,tf,iOption,cOption )
  %%
  % Nx,Ny: number cells
  % nstep: number of time steps
  % nu: input parameter
  % tf: final time
  % iOption: selects TZ function
  % cOption: selects coefficient variation
   
  %% set some parameters
  xa    = 0;  % left of domain
  xb    = 1;  % right of domain
  ya    = 0;  % bottom of domain
  yb    = 1;  % rom of domain
  dx    = (xb-xa)/(Nx);
  dy    = (yb-ya)/(Ny);
  dt    = tf/nStep;
  
  %% number of grid cells in x and y (note we have no ghost cells here)
  NXT   = Nx+1;
  NYT   = Ny+1;
  %% set up spatial grid
  x     = linspace( xa,xb,NXT );
  y     = linspace( ya,yb,NYT );
  
  rx    = dt/dx^2;
  ry    = dt/dy^2;
  
  u = zeros(NXT,NYT);
  %% set ICs
  t = 0;
  for k = 1:NYT
    for j = 1:NXT
      u(j,k) = getEX( x(j), y(k), t, iOption );
    end
  end
  
  D = zeros(NXT,NYT);
  %% set variable coefficients
  t = 0;
  for k = 1:NYT
    for j = 1:NXT
      D(j,k) = getD( x(j), y(k), cOption );
    end
  end
  
  %% allocate discretization matrices and RHSs
  Ax = zeros(NXT,NXT);
  qx = zeros(NXT,1);
  Ay = zeros(NYT,NYT);
  qy = zeros(NYT,1);
  
  %% time-stepping loop
  told = 0;
  for n = 1:nStep
    tnew  = told+dt;
    thalf = 0.5*(told+tnew);
        
    %% 1st half step (implicit in x)
    uold = u;

    %% loop over k
    for k = 1:NYT
      %% set exact solution at y=0, y=1
      if( k == 1 | k == NYT )
        for j = 1:NXT
          u(j,k) = getEX( x(j), y(k), thalf, iOption );
        end
      else
      
        %% otherwise do discretization
        Ax(1,1) = 1.;
        qx(1)   = getEX( x(1), y(k), thalf, iOption );
        Ax(NXT,NXT) = 1.;
        qx(NXT)     = getEX( x(NXT), y(k), thalf, iOption );
        for j = 2:NXT-1
          %% metric averages
          Dmx = 0.5*(D(j-1,k)+D(j,k));
          Dpx = 0.5*(D(j,k)+D(j+1,k));
          Dmy = 0.5*(D(j,k-1)+D(j,k));
          Dpy = 0.5*(D(j,k)+D(j,k+1));
  
          %% matrix entries
          Ax(j,j-1) = -0.5*rx*Dmx;
          Ax(j,j)   = 1.+0.5*rx*(Dmx+Dpx);
          Ax(j,j+1) = -0.5*rx*Dpx;
        
          %% RHS entry
          uyp = uold(j,k+1)-uold(j,k);
          uym = uold(j,k)-uold(j,k-1);
        
          qx(j) = uold(j,k)+0.5*ry*(Dpy*uyp-Dmy*uym)+...
            0.5*dt*getQ( x(j),y(k),told+0.25*dt,iOption,cOption );
        end
        u(:,k) = Ax\qx;  
      end
    end
    
    %% 2nd half step (implicit in y)
    uold = u;

    %% loop over j
    for j = 1:NXT
      %% set exact solution at x=0, x=1
      if( j == 1 | j == NXT )
        for k = 1:NYT
          u(j,k) = getEX( x(j), y(k), tnew, iOption );
        end
      else
      
        %% otherwise do discretization
        Ay(1,1) = 1.;
        qy(1)   = getEX( x(j), y(1), tnew, iOption );
        Ay(NYT,NYT) = 1.;
        qy(NYT)     = getEX( x(j), y(NYT), tnew, iOption );
        for k = 2:NYT-1
          %% metric averages
          Dmx = 0.5*(D(j-1,k)+D(j,k));
          Dpx = 0.5*(D(j,k)+D(j+1,k));
          Dmy = 0.5*(D(j,k-1)+D(j,k));
          Dpy = 0.5*(D(j,k)+D(j,k+1));
  
          %% matrix entries
          Ay(k,k-1) = -0.5*ry*Dmy;
          Ay(k,k)   = 1.+0.5*ry*(Dmy+Dpy);
          Ay(k,k+1) = -0.5*ry*Dpy;
        
          %% RHS entry
          uxp = uold(j+1,k)-uold(j,k);
          uxm = uold(j,k)-uold(j-1,k);
        
          qy(k) = uold(j,k)+0.5*rx*(Dpx*uxp-Dmx*uxm)+...
            0.5*dt*getQ( x(j),y(k),told+0.75*dt,iOption,cOption );
        end
        u(j,:) = (Ay\qy)';  
      end
    end
      
    told = tnew;
  end
  
  %% compute exact solution
  ue = zeros(NXT,NYT);
  for k = 1:NYT
    for j = 1:NXT
      ue(j,k) = getEX( x(j), y(k), tnew, iOption );
    end
  end
  
  return
end

function Q = getQ( x, y, t, iOption,cOption )
  %% function to compute TZ forcing
  ut  = getEX_t(  x,y,t,iOption );
  uxx = getEX_xx( x,y,t,iOption );
  uyy = getEX_yy( x,y,t,iOption );
  ux  = getEX_x(  x,y,t,iOption );
  uy  = getEX_y(  x,y,t,iOption );
  
  Dx  = getD_x(   x,y,cOption );
  Dy  = getD_y(   x,y,cOption );
  D   = getD(     x,y,cOption );
  
  Q = ut-(Dx*ux+D*uxx+Dy*uy+D*uyy);
  return
end

function z = getD( x, y, cOption )
  %% function to set variable coefficients
  if( cOption == 1 )
    z = 1.;
  elseif( cOption == 2 )
    z = 1.+x+y;
  else
    k = 2.*pi;
    z = 1.+0.1*cos(k*x)*sin(k*y);
  end
  return
end

function z = getD_x( x, y, cOption )
  %% function to set variable coefficient derivatie
  if( cOption == 1 )
    z = 0.;
  elseif( cOption == 2 )
    z = 1.;
  else
    k = 2.*pi;
    z = -0.1*k*sin(k*x)*sin(k*y);
  end
  return
end

function z = getD_y( x, y, cOption )
  %% function to set variable coefficient derivative
  if( cOption == 1 )
    z = 0.;
  elseif( cOption == 2 )
    z = 1.;
  else
    k = 2.*pi;
    z = 0.1*k*cos(k*x)*cos(k*y);
  end
  return
end

function z = getEX( x, y, t, iOption )
  %% function to compute exact TZ solution
  if( iOption == 1 )
    z = 1.+x+x^2+y+y^2+t+t^2;
  elseif( iOption == 2 )
    z = 1.+x+y+t;
  else
    z = sin(x)*log(1.+y)*exp(-t);
  end
  return
end

function z = getEX_t( x, y, t, iOption )
  %% function to compute exact TZ solution (t derivative)
  if( iOption == 1 )
    z = 1.+2.*t;
  elseif( iOption == 2 )
    z = 1.;
  else
    z = -sin(x)*log(1.+y)*exp(-t);
  end
  return
end

function z = getEX_x( x, y, t, iOption )
  %% function to compute exact TZ solution (x derivative)
  if( iOption == 1 )
    z = 1.+2.*x;
  elseif( iOption == 2 )
    z = 1.;
  else
    z = cos(x)*log(1.+y)*exp(-t);
  end
  return
end

function z = getEX_xx( x, y, t, iOption )
  %% function to compute exact TZ solution (xx derivative)
  if( iOption == 1 )
    z = 2.;
  elseif( iOption == 2 )
    z = 0.;
  else
    z = -sin(x)*log(1.+y)*exp(-t);
  end
  return
end

function z = getEX_y( x, y, t, iOption )
  %% function to compute exact TZ solution (y derivative)
  if( iOption == 1 )
    z = 1.+2.*y;
  elseif( iOption == 2 )
    z = 1.;
  else
    z = sin(x)/(1.+y)*exp(-t);
  end
  return
end

function z = getEX_yy( x, y, t, iOption )
  %% function to compute exact TZ solution (yy derivative)
  if( iOption == 1 )
    z = 2.;
  elseif( iOption == 2 )
    z = 0.;
  else
    z = -sin(x)/((1.+y)^2)*exp(-t);
  end
  return
end
