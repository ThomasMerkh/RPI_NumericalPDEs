function [x,y,u,ue] = Heat2DCrankNicolsonDD( Nx,Ny,nStep,nu,tf,iOption )
  %%
  % Nx,Ny: number cells
  % nstep: number of time steps
  % nu: input parameter
  % tf: final time
  % iOption: selects TZ function
  
  %% sample runs, see runHeat2DCNDD.m
 
  %% set some parameters
  xa    = 0;  % left of domain
  xb    = 1;  % right of domain
  ya    = 0;  % bottom of domain
  yb    = 1;  % rom of domain
  dx    = (xb-xa)/(Nx);
  dy    = (yb-ya)/(Ny);
  dt    = tf/nStep;
  
  %% set up spatial grid
  x     = linspace( xa,xb,Nx+1 );
  y     = linspace( ya,yb,Ny+1 );
  
  Ntot  = (Nx+1)*(Ny+1);
  rx    = nu*dt/dx^2;
  ry    = nu*dt/dy^2;
  
  u = zeros(Nx+1,Ny+1);
  %% set ICs
  t = 0;
  for k = 1:Ny+1
    for j = 1:Nx+1
      u(j,k) = getEX( x(j), y(k), t, iOption );
    end
  end
  
  %% allocate discretization matrix and RHS
  A = zeros(Ntot,Ntot);
  q = zeros(Ntot,1);
  
  %% setup discretization matrix
  % Basically, the first Nx+1 diag values of A are 1 corresponding to the
  % (y=0) boundary condition.  Same with the last Nx+1 diag values,
  % corrsponding to (y = 1) boundary condition.  Every apprx Nx+1 values in
  % the interior will be 1's corresponding to the boundaries at x=0 and
  % x=1, since we are "climbing" up the unit box and every so often we hit
  % the boundary at x=0 or x=1.  Coundary conditions correspond to a rows
  % with diag=1 and the rest equal to zero.  Interior points will have up
  % to 5 non-zero values, on a tridiag (left center right) and 2 offset
  % diag (up and down values around a point).
  
  % BCs at y = 0 
  k = 1; 
  for j = 1:Nx+1
    ind = getInd( j,k,Nx+1,Ny+1 );
    A(ind,ind) = 1.0;
  end
  
  % loop over lines of constant y
  for k = 2:Ny
    % BCs at x=0
    j = 1;
    ind = getInd( j,k,Nx+1,Ny+1 );
    A(ind,ind) = 1.0;
    
    % interior
    for j = 2:Nx
      ind = getInd( j,k,Nx+1,Ny+1 );
      A(ind,getInd( j,  k,Nx+1,Ny+1 )) = 1.+rx+ry;
      A(ind,getInd( j-1,k,Nx+1,Ny+1 )) = -rx/2.;
      A(ind,getInd( j+1,k,Nx+1,Ny+1 )) = -rx/2.;
      A(ind,getInd( j,k-1,Nx+1,Ny+1 )) = -ry/2.;
      A(ind,getInd( j,k+1,Nx+1,Ny+1 )) = -ry/2.;
    end
    
    % BCs at x=1
    j = Nx+1;
    ind = getInd( j,k,Nx+1,Ny+1 );
    A(ind,ind) = 1.0;
  end
  
  % BCs at y = 1
  k = Ny+1; 
  for j = 1:Nx+1
    ind = getInd( j,k,Nx+1,Ny+1 );
    A(ind,ind) = 1.0;
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
    for j = 1:Nx+1
      ind = getInd( j,k,Nx+1,Ny+1 );
      q(ind,1) = getEX( x(j), y(k), tnew, iOption );
    end
    
    for k = 2:Ny
      % BCs at x=0
      j = 1;
      ind = getInd( j,k,Nx+1,Ny+1 );
      q(ind,1) = getEX( x(j), y(k), tnew, iOption );
      
      % interior
      for j = 2:Nx
        ind = getInd( j,k,Nx+1,Ny+1 );
        q(ind,1) = u(j,k)+...
          0.5*rx*(u(j+1,k)-2.*u(j,k)+u(j-1,k))+...
          0.5*ry*(u(j,k+1)-2.*u(j,k)+u(j,k-1))+...
          +dt*getQ( x(j),y(k),thalf,nu,iOption );
      end
      
      % BCs at x=1
      j = Nx+1;
      ind = getInd( j,k,Nx+1,Ny+1 );
      q(ind,1) = getEX( x(j), y(k), tnew, iOption );
    end
      
    % BCs at y=1
    k = Ny+1;
    for j = 1:Nx+1
      ind = getInd( j,k,Nx+1,Ny+1 );
      q(ind,1) = getEX( x(j), y(k), tnew, iOption );
    end
    
    %% now do solve 
    p = U\(L\(P*q));
    
    % translate new solution from vector to array
    for j = 1:Nx+1
      for k = 1:Ny+1
        u(j,k) = p(getInd( j,k,Nx+1,Ny+1 ),1);
      end
    end
    
    told = tnew;
  end
  
  %% compute exat solution
  ue = zeros(Nx+1,Ny+1);
  for k = 1:Ny+1
    for j = 1:Nx+1
      ue(j,k) = getEX( x(j), y(k), tnew, iOption );
    end
  end
  
  return
end

function u = getEX( x, y, t, iOption )
  %% function to compute exact TZ solution
  if( iOption == 1 )
    u = 1.+x+x^2+y+y^2+t+t^2;
  else
    u = sin(x)*log(1.+y)*exp(-t);
  end
  return
end

function u = getQ( x, y, t, nu, iOption )
  %% function to compute TZ forcing
  if( iOption == 1 )
    u = (2.*t+1.)-4.*nu;
  else
    u = sin(x)*exp(-t)*(nu/((1.+y)^2)+log(1.+y)*(nu-1.));
  end
  return
end

function ind = getInd( j,k,Nxtot,Nytot )
  %% indexing function
  ind = (k-1)*Nxtot+j;
  %ind = (j-1)*Nytot+k;
  return
end