function [x,y,u] = Heat2DAnnulus( Nr,Ns,nu,tf,iOption,iTZ,gamma )
  %%
  % Nr,Ns: number cells
  % nStep: number of time steps
  % nu: input parameter
  % tf: final time
  % iOption: selects case (1: TZ, 2: topHat)
  % gamma: selects time stepper (1: FE, 0: BE, .5: CN)
     
  % Note: r: radius
  %       s: theta
  %       t: time
    
  %% set some parameters
  ra    = .5;  % inner radius
  rb    = 1.;  % outter radius
  sa    = 0.;  
  sb    = 2.*pi; 
  dr    = (rb-ra)/(Nr);
  ds    = (sb-sa)/(Ns);
  
  %% set some constants for the time stepper
  gl    = 1.0-gamma;
  gr    = gamma;
  
  %% estimate time step
  circ = pi*2.*ra;
  h = min( dr,circ/Ns );
  if( gamma > 0.5 )
    dtGuess = h^2/(4.*nu);
  else
    dtGuess = h/nu;
  end
  nStep = ceil(tf/dtGuess);
  dt = tf/nStep;
  
  NRT   = Nr+1;   %% no ghost in radial
  NST   = Ns+1+2; %% 2 ghost for periodic BCs in theta
  %% set up polar grid
  r     = linspace( ra,rb,NRT );
  s     = linspace( sa-ds,sb+ds,NST );
  
  Ntot  = NRT*NST;

  u = zeros(NRT,NST);
  x = zeros(NRT,NST);
  y = zeros(NRT,NST);
  %% set ICs and physical grid
  t = 0;
  for k = 1:NST
    for j = 1:NRT
      x(j,k) = r(j)*cos(s(k));
      y(j,k) = r(j)*sin(s(k));
      u(j,k) = getIC( r(j), s(k), t, iOption, iTZ );
    end
  end
  
  surf( x,y,u );
  AXXX = axis;
  shading interp
  xlabel( 'x' );
  ylabel( 'y' );
  pause
  
  %% allocate discretization matrix and RHS
  A = zeros(Ntot,Ntot);
  q = zeros(Ntot,1);
  
  %% setup discretization matrix
  % periodic BCs
  for j = 1:NRT
    eqn = getInd( j,1,NRT,NST );
    A(eqn,getInd( j,1,NRT,NST )) = 1.;
    A(eqn,getInd( j,NST-2,NRT,NST )) = -1.;
    
    eqn = getInd( j,NST,NRT,NST );
    A(eqn,getInd( j,NST,NRT,NST )) = 1.;
    A(eqn,getInd( j,3,NRT,NST )) = -1.;
  end
  
  % loop over lines of constant s
  for k = 2:NST-1
    % BCs at r=a
    j = 1;
    eqn = getInd( j,k,NRT,NST );
    A(eqn,getInd( j,k,NRT,NST )) = 1.;
    
    % interior
    for j = 2:NRT-1
      eqn = getInd( j,k,NRT,NST );
      rjp = 0.5*(r(j+1)+r(j));
      rjm = 0.5*(r(j)+r(j-1));
      A(eqn,getInd( j,  k,NRT,NST )) = 1.+...
        gl*dt*nu/(r(j)*dr^2)*(rjp+rjm)+...
        gl*2.*dt*nu/(r(j)^2*ds^2);
      A(eqn,getInd( j-1,k,NRT,NST )) = -gl*dt*nu*rjm/(r(j)*dr^2);
      A(eqn,getInd( j+1,k,NRT,NST )) = -gl*dt*nu*rjp/(r(j)*dr^2);
      A(eqn,getInd( j,k-1,NRT,NST )) = -gl*dt*nu/(r(j)^2*ds^2);
      A(eqn,getInd( j,k+1,NRT,NST )) = -gl*dt*nu/(r(j)^2*ds^2);
    end
    
    % BCs at r=b
    j = NRT;
    eqn = getInd( j,k,NRT,NST );
    A(eqn,getInd( j,k,NRT,NST )) = 1.;
  end
  
  
  %% factor matrix
  [L,U,P] = lu(A);
  
  %% time-stepping loop
  told = 0;
  for n = 1:nStep
    tnew  = told+dt;
    thalf = 0.5*(told+tnew);

    %% create RHS vector
    % periodic BCs
    for j = 1:NRT
      eqn = getInd( j,1,NRT,NST );
      q(eqn,1) = 0.;
      
      eqn = getInd( j,NST,NRT,NST );
      q(eqn,1) = 0.;
    end
    
    for k = 2:NST-1
      % BCs at r=a
      j = 1;
      eqn = getInd( j,k,NRT,NST );
      q(eqn,1) = alpha( r(j), s(k), tnew, iOption,iTZ );
      
      % interior
      for j = 2:NRT-1
        eqn = getInd( j,k,NRT,NST );
        rjp = 0.5*(r(j+1)+r(j));
        rjm = 0.5*(r(j)+r(j-1));
        q(eqn,1) = ...
          u(j,k)*(1.-gr*dt*nu/(r(j)*dr^2)*(rjp+rjm)-gr*2.*dt*nu/(r(j)^2*ds^2))+...
          u(j-1,k)*gr*dt*nu/(r(j)*dr^2)*rjm+...
          u(j+1,k)*gr*dt*nu/(r(j)*dr^2)*rjp+...
          u(j,k-1)*gr*dt*nu/(r(j)^2*ds^2)+...
          u(j,k+1)*gr*dt*nu/(r(j)^2*ds^2);
                
        %% add forcing if doing TZ
        if( iOption == 1 )
          q(eqn,1) = q(eqn,1)+dt*getQ( r(j),s(k),thalf,nu,iTZ );
        end
      end
      
      % BCs at r=b
      j = NRT;
      eqn = getInd( j,k,NRT,NST );
      q(eqn,1) = beta( r(j), s(k), tnew, iOption,iTZ );
    end
        
    %% now do solve 
    p = U\(L\(P*q));
 
    % translate new solution from vector to array
    for k = 1:NST
      for j = 1:NRT
        u(j,k) = p(getInd( j,k,NRT,NST ),1);
      end
    end
    
    if( iOption == 1 )
      ue = 0*u;
      for k = 1:NST
        for j = 1:NRT
          ue(j,k) = getU( r(j),s(k),tnew,iTZ );
        end
      end
      surf( x,y,u-ue );
      xlabel( 'x' );
      ylabel( 'y' );
      shading interp
      pause( 0.01 );
    else
      surf( x,y,u );
      axis([AXXX]);  % Comment out for adjusting axis.
      xlabel( 'x' );
      ylabel( 'y' );
      %shading interp
      %view([0,90]);
      pause
      pause( 0.01 );
    end
    
    told = tnew;
  end
    
  
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

function z = getIC( r,s,t,iOption,iTZ )
  %% function to compute initial condition
  if( iOption == 1 ) %% TZ
    z = getU( r,s,t,iTZ );
  elseif( iOption == 2 ) %% topHat
    x = r*cos(s);
    y = r*sin(s);
    x0 = 0.;
    y0 = .75;
    rad2 = (x-x0)^2+(y-y0)^2;
    if( rad2 < .15^2 )
      z = 1.e8;
    else
      z = 0.;
    end
  else %% wall hot spot
    z = 0.0;
  end
  return
end


function a = alpha( r,s, t, iOption,iTZ )
  %% function to compute alpha (solution on inner radius)
  if( iOption == 1 )
    a = getU( r,s,t,iTZ );
  elseif( iOption == 2 )
    a = 0.;
  else
    if( r*cos(s) > 0 )
      a = 1.0;
    else
      a = 0.0;
    end
  end
  return
end

function b = beta( r, s, t, iOption,iTZ )
  %% function to compute beta (solution on outter radius)
  if( iOption == 1 )
    b = getU( r,s,t,iTZ );
  elseif( iOption == 2 )
    b = 0.;
  else
    if( r*sin(s) < 0 )
      b = -1.0;
    else
      b = 0.0;
    end
  end
  return
end

function z = getU( r, s, t,iTZ )
  %% function to compute exact TZ solution
  if( iTZ == 1 )
    z = sin(r)*cos(s)*cos(t);
  else
    z = 1;
  end
  return
end

function z = getU_t( r,s,t,iTZ )
  %% function to compute exact TZ solution (t derivative)
  if( iTZ == 1 )
    z = -sin(r)*cos(s)*sin(t);
  else
    z = 0.;
  end
  return
end

function z = getU_r( r,s,t,iTZ )
  %% function to compute exact TZ solution (r derivative)
  if( iTZ == 1 )
    z = cos(r)*cos(s)*cos(t);
  else
    z = 0.;
  end
  return
end

function z = getU_rr( r,s,t,iTZ )
  %% function to compute exact TZ solution (rr derivative)
  if( iTZ == 1 )
    z = -sin(r)*cos(s)*cos(t);
  else
    z = 0.;
  end
  return
end

function z = getU_s( r,s,t,iTZ )
  %% function to compute exact TZ solution (s derivative)
  if( iTZ == 1 )
    z = -sin(r)*sin(s)*cos(t);
  else
    z = 0.;
  end
  return
end

function z = getU_ss( r, s, t,iTZ )
  %% function to compute exact TZ solution (ss derivative)
  if( iTZ == 1 )
    z = -sin(r)*cos(s)*cos(t);
  else
    z = 0.;
  end
  return
end

function Q = getQ( r,s,t,nu,iTZ )
  %% function to compute TZ forcing
  ut  = getU_t(  r,s,t,iTZ );
  ur  = getU_r(  r,s,t,iTZ );
  urr = getU_rr( r,s,t,iTZ );
  us  = getU_s(  r,s,t,iTZ );
  uss = getU_ss( r,s,t,iTZ );
  
  Q = ut-nu*(1./r*(ur+r*urr)+1./(r^2)*uss);
  return
end

function ind = getInd( j,k,Nrtot,Nstot )
  %% indexing function
  ind = (k-1)*Nrtot+j;
  %ind = (j-1)*Nstot+k;
  return
end
  