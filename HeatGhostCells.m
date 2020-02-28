function []=HeatGhostCells( M,nu,r0,tf )
  %%
  % M: number of interior cells
  % nu: input parameter
  % r0: nu*dt/dx^2 (guess)
  % tf: final time
 
  %% set some parameters
  N     = M+2; % total number of grid points
  xL    = -1;  % left of domain
  xR    = 1;   % right of domain
  dx    = (xR-xL)/(M-1);
  
  %% set up spatial grid
  x     = linspace( xL-dx,xR+dx,N );
  
  %% calculate the number of time steps
  %%  including some safety in case dt doesn't fit an integral # of steps
  t     = 0.;  % current time
  dt    = r0*dx^2/nu;
  nStep = ceil((tf-t)/dt);
  dt    = (tf-t)/nStep;
  r     = nu*dt/dx^2;
  

  %% allocate space
  u     = zeros(N,1);
  uold  = zeros(N,1);

  %% set initial condition (Heaviside function)
  for j = 1:N
    if( x(j) < 0 )
      uold(j) = 0.0;
    else
      uold(j) = 1.0;
    end
  end

  %% time stepping loop
  for n = 0:nStep


    %% space loop (do not include ghost cells)
    for j = 2:N-1
      %% FD formula
      u(j) = uold(j)+r*(uold(j+1)-2.*uold(j)+uold(j-1));
    end
    
    t = t+dt;
    
    %% boundary conditions (at new time)
    uxL  = pi ^ (-0.1e1 / 0.2e1) * exp(-xL ^ 2 / nu / t / 0.4e1) * (nu * t) ^ (-0.1e1 / 0.2e1) / 0.2e1;
    u(1) = u(3)-2.0*dx*uxL;
    
    utR  = -pi ^ (-0.1e1 / 0.2e1) * exp(-xR ^ 2 / nu / t / 0.4e1) * xR * (nu * t) ^ (-0.3e1 / 0.2e1) * nu / 0.4e1;
    u(N) = 2.*u(N-1)-u(N-2)+dx^2/nu*utR;
    
    %% update old soln
    uold = u;
    
    %% make some plots
    uex = 0.5*(1.0+erf(x/(sqrt(4*nu*t))));
    plot( x(2:N-1),u(2:N-1),'bx', x(2:N-1),uex(2:N-1),'k-' );
    %plot( x,u,'bx', x,uex,'k-' );
    xlabel( 'x' );
    ylabel( 'solution' );
    pause(0.01);
    if( n == 0 )
      pause
    end
      
  end

  