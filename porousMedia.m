function []=porousMedia( N,L,gam,nPlots,tf )
  % N: number cells
  % L: length of domain
  % gam: exponent (as in (u^gam u_x)_x )
  % tf: final time
  
  %% sample runs
  % porousMedia( 200,2,.5,5,.1 )
  % porousMedia( 200,2,.9,5,.1 )
  % porousMedia( 200,2,1,5,.1 )
  % porousMedia( 200,2,1.5,5,.1 )
  % porousMedia( 200,2,1.9,5,.1 )
  % porousMedia( 200,2,2.0,5,.1 )
  %

  %% set some parameters
  xa    = 0;  % left of domain
  xb    = L;  % right of domain
  dx    = (xb-xa)/(N);
  nStep = N*10;
  dt    = tf/nStep;
  r = (gam+1)*dt/(2.^gam*dx^2);
  
  %% set some parameters for the Newton iteration
  maxIt = 10;
  tol   = 1e-10;
  
  %% set tolerance to measure zero concentration front
  uMin   = 1e-2;
  xFront = [];
  tFront = [];
  
  %% set up the grid, ICs, and matrix
  Ntot = N+2;
  x = linspace( xa-0.5*dx,xb+0.5*dx,Ntot );
  u = zeros(Ntot,1);
  A = zeros(Ntot,Ntot);
  f = zeros(Ntot,1);
  for j = 1:length(x)
    A = getNormalization(gam);
    u(j) = A*(5./((cosh(5*x(j)))^2)); 
  end
  %% set BCs on initial conditions
  u(1)    = u(2);
  u(Ntot) = u(Ntot-1);
  
  %% plot initial conditions
  fs = 16;
  lineWidth = 2;
  ms = 16;    
  pInd = 2:Ntot-1;
  figure(1);
  set(gca,'FontSize',fs);
  plot( x(pInd),u(pInd),'b-' );
  xlabel( 'x' );
  ylabel( 'u' );
  figure(2);
  set(gca,'FontSize',fs);
  plot( [0],[0] );
  xlabel( 'time' );
  ylabel( 'front location' );
  pause
  
  
  %%%%%
  %% time-stepping loop
  %%%%%
  told = 0;
  for n = 1:nStep
    tnew  = told+dt;
    
    uold = u;
    %% Newton iteration
    for it = 1:maxIt
      %%
      %% setup Jacobian matrix
      %%
      % Neumann BCs
      A(1,1) = 1.;
      A(1,2) = -1.;
      f(1) = 0.;
      A(Ntot,Ntot)   = 1.;
      A(Ntot,Ntot-1) = -1.;
      f(Ntot) = 0.;
      
      % interior equations
      for j = 2:Ntot-1
        A(j,j+1) = -0.5*r*(...
          gam*(u(j+1)+u(j))^(gam-1.)*(u(j+1)-u(j))...
          +(u(j+1)+u(j))^gam);
        A(j,j)   = 1.0-0.5*r*(...
          gam*(u(j+1)+u(j))^(gam-1.)*(u(j+1)-u(j))...
          -(u(j+1)+u(j))^gam...
          -gam*(u(j)+u(j-1))^(gam-1.)*(u(j)-u(j-1))...
          -(u(j)+u(j-1))^gam);
        A(j,j-1) = -0.5*r*(...
          -gam*(u(j)+u(j-1))^(gam-1.)*(u(j)-u(j-1))...
          +(u(j)+u(j-1))^gam);
        
        f(j) = u(j)-0.5*r*(...
            (u(j+1)+u(j))^gam*(u(j+1)-u(j))-...
            (u(j)+u(j-1))^gam*(u(j)-u(j-1)))...
          -uold(j)-0.5*r*(...
            (uold(j+1)+uold(j))^gam*(uold(j+1)-uold(j))-...
            (uold(j)+uold(j-1))^gam*(uold(j)-uold(j-1)));
      end
      
      %% compute solution increment and stop Newton if within tol
      du = A\f;
      u = u-du;
      res = max(abs(du));
      if( res < tol )
        break;
      else
        if( it == maxIt )
          fprintf( 'Iteration did not converge, t=%e, n=%i\n',told,n );
          return
        end
      end
    end
    
    %%%%%
    %% The remainder of this code is essentially plotting 
    %%%%%
    
    %% estimate location of zero u
    for j = 2:Ntot-1
      if u(j) < uMin
        m = (u(j)-u(j-1))/dx;
        xFront = [xFront,x(j-1)+(uMin-u(j-1))/m];
        tFront = [tFront,tnew];
        break
      end
    end
    
    if( mod(n,ceil(nStep/nPlots)) == 0 )      
      alpha = 1/(gam+2);
      xp = x(pInd);
      t = tnew;
      figure(1);
      pInd = 2:Ntot-1;
      hold on
      % plot numerical solution
      plot( xp,u(pInd),'b-', 'lineWidth',lineWidth,'MarkerSize',ms );
      % plot Barenblatt solution
      xf = sqrt(2)*sqrt(alpha*gam*t^(2*alpha)*(gam+1))/(alpha*gam);
      xp = linspace(0,xf-10*eps,N );
      plot( xp,t^(-alpha)*(1-alpha*gam/(2*(gam+1))*xp.^2./(t^(2*alpha))).^(1/gam),'k', 'lineWidth',lineWidth,'MarkerSize',ms );
      hold off
      xlabel( 'x' );
      ylabel( 'u' );
      
      figure(2);
      % plot estimated front location
      plot( tFront,xFront,'rx', 'lineWidth',lineWidth,'MarkerSize',ms );
      hold on
      % plot front location from Barenblatt
      plot( tFront, sqrt(2)*sqrt(alpha*gam*tFront.^(2*alpha)*(gam+1))/(alpha*gam), 'k', 'lineWidth',lineWidth,'MarkerSize',ms );
      hold off
      xlabel( 'time' );
      ylabel( 'front location' );
      pause(0.01);
      pause
    end
 
    told = tnew;
        
  end
  
  return
end

function A = getNormalization( g )
  if( abs(g-2) < 10*eps )
    A = (1/4)*pi*sqrt(2)*sqrt(6);
  else
    A = -pi ^ (-0.1e1 / 0.2e1) * sqrt(0.2e1) * ((-g / (g + 2) / (g + 1)) ^ (-0.1e1 / 0.2e1)) * gamma(-((g + 2) / g) / 0.2e1) / g * gamma((1 / g * (g + 1))) * (sqrt((1 / (g + 2) * g * (g + 1))) * ((-g / (g + 2) / (g + 1)) ^ (((g + 2) / g) / 0.2e1)) * ((g + 2) / g * sqrt((1 / (g + 2) * g * (g + 1)))) ^ (2 / g) * g + g * sin(pi / g) + 0.2e1 * sqrt((1 / (g + 2) * g * (g + 1))) * ((-g / (g + 2) / (g + 1)) ^ (((g + 2) / g) / 0.2e1)) * ((g + 2) / g * sqrt((1 / (g + 2) * g * (g + 1)))) ^ (2 / g)) / 0.2e1;
    A = real(A);
  end
end
  