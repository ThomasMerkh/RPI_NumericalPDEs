function [x,y,u,ue] = Heat2DCrankNicolsonGhost( Nx,Ny,nStep,nu,tf,iOption )
  %%
  % Nx,Ny: number cells
  % nstep: number of time steps
  % nu: input parameter
  % tf: final time
  % iOption: selects TZ function
  
  %% sample runs, see 
 
  %% set some parameters
  xa    = 0;  % left of domain
  xb    = 1;  % right of domain
  ya    = 0;  % bottom of domain
  yb    = 1;  % rom of domain
  dx    = (xb-xa)/(Nx);
  dy    = (yb-ya)/(Ny);
  dt    = tf/nStep;
  
  NXT   = Nx+1+2;
  NYT   = Ny+1+2;
  %% set up spatial grid
  x     = linspace( xa-dx,xb+dx,NXT );
  y     = linspace( ya-dy,yb+dy,NYT );
  
  Ntot  = NXT*NYT;
  rx    = nu*dt/dx^2;
  ry    = nu*dt/dy^2;
  
  u = zeros(NXT,NYT);
  %% set ICs
  t = 0;
  for k = 1:NYT
    for j = 1:NXT
      u(j,k) = getEX( x(j), y(k), t, iOption );
    end
  end
  
  %% allocate discretization matrix and RHS
  A = zeros(Ntot,Ntot);
  q = zeros(Ntot,1);
  
  %% setup discretization matrix
  % BCs at y = 0 
  k = 1; 
  for j = 1:NXT
    ind = getInd( j,k,NXT,NYT );
    A(ind,getInd( j,k,  NXT,NYT )) = 0.5;
    A(ind,getInd( j,k+2,NXT,NYT )) = 0.5;
  end
  
  % loop over lines of constant y
  for k = 2:Ny+2
    % BCs at x=0
    j = 1;
    ind = getInd( j,k,NXT,NYT );
    A(ind,getInd( j,k,  NXT,NYT )) = -1./(2.*dx);
    A(ind,getInd( j+2,k,NXT,NYT )) = 1./(2.*dx);
    
    % interior
    for j = 2:Nx+2
      ind = getInd( j,k,NXT,NYT );
      A(ind,getInd( j,  k,NXT,NYT )) = 1.+rx+ry;
      A(ind,getInd( j-1,k,NXT,NYT )) = -rx/2.;
      A(ind,getInd( j+1,k,NXT,NYT )) = -rx/2.;
      A(ind,getInd( j,k-1,NXT,NYT )) = -ry/2.;
      A(ind,getInd( j,k+1,NXT,NYT )) = -ry/2.;
    end
    
    % BCs at x=1
    j = NXT;
    ind = getInd( j,k,NXT,NYT );
    A(ind,getInd( j,k,  NXT,NYT )) = 0.5;
    A(ind,getInd( j-2,k,NXT,NYT )) = 0.5;
  end
  
  % BCs at y = 1
  k = NYT; 
  for j = 1:NXT
    ind = getInd( j,k,NXT,NYT );
    A(ind,getInd( j,k,  NXT,NYT )) = 1./(2.*dy);;
    A(ind,getInd( j,k-2,NXT,NYT )) = -1./(2.*dy);;
  end
  
  %% factor matrix
  [L,U,P] = lu(A);
  
  %% time-stepping loop
  told = 0;
  for n = 1:nStep
    tnew  = told+dt;
    thalf = 0.5*(told+tnew);

    %% create RHS vector
    % BCs at y=0
    k = 1;
    for j = 1:NXT
      ind = getInd( j,k,NXT,NYT );
      q(ind,1) = getb0( x(j), tnew, iOption );
    end
    
    for k = 2:Ny+2
      % BCs at x=0
      j = 1;
      ind = getInd( j,k,NXT,NYT );
      q(ind,1) = geta0( y(k), tnew, iOption );
      
      % interior
      for j = 2:Nx+2
        ind = getInd( j,k,NXT,NYT );
        q(ind,1) = u(j,k)+...
          0.5*rx*(u(j+1,k)-2.*u(j,k)+u(j-1,k))+...
          0.5*ry*(u(j,k+1)-2.*u(j,k)+u(j,k-1))+...
          +dt*getQ( x(j),y(k),thalf,nu,iOption );
      end
      
      % BCs at x=1
      j = NXT;
      ind = getInd( j,k,NXT,NYT );
      q(ind,1) = geta1( y(k), tnew, iOption );
    end
      
    % BCs at y=1
    k = NYT;
    for j = 1:NXT
      ind = getInd( j,k,NXT,NYT );
      q(ind,1) = getb1( x(j), tnew, iOption );
    end
    
    %% now do solve 
    p = U\(L\(P*q));
    
    % translate new solution from vector to array
    for j = 1:NXT
      for k = 1:NYT
        u(j,k) = p(getInd( j,k,NXT,NYT ),1);
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
  
  indx = 2:Nx+2;
  indy = 2:Ny+2;
  u = u(indx,indy);
  ue = ue(indx,indy);
  x = x(indx);
  y = y(indy);
  
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

function a0 = geta0( y, t, iOption )
  %% function to compute a0
  a0 = getEX_x( 0.,y,t,iOption );
  return
end

function a1 = geta1( y, t, iOption )
  %% function to compute a1
  a1 = getEX( 1.,y,t,iOption );
  return
end

function b0 = getb0( x, t, iOption )
  %% function to compute b0
  b0 = getEX( x,0.,t,iOption );
  return
end

function b1 = getb1( x, t, iOption )
  %% function to compute b1
  b1 = getEX_y( x,1.,t,iOption );
  return
end

function Q = getQ( x, y, t, nu, iOption )
  %% function to compute TZ forcing
  ut = getEX_t( x,y,t,iOption );
  uxx = getEX_xx( x,y,t,iOption );
  uyy = getEX_yy( x,y,t,iOption );
  
  Q = ut-nu*(uxx+uyy);
  return
end

function ind = getInd( j,k,Nxtot,Nytot )
  %% indexing function
  ind = (k-1)*Nxtot+j;
  %ind = (j-1)*Nytot+k;
  return
end
  