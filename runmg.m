function [v0,res] = runmg( N,M,nSmooth,nLevel,nCycle )
  dx = 1/N;
  dy = 1/M;

  x = zeros(N+1,M+1);
  y = zeros(N+1,M+1);
  for j = 0:N
    for k = 0:M
      x(j+1,k+1) = j*dx;
      y(j+1,k+1) = k*dy;
    end
  end

  uex = sin(2.*pi*x).*sin(2.*pi*y);
  f = -8.*pi^2*uex;
  
  %uex = sinh(2.*x).*sin(y);
  %f = 3.*uex;
  
  %% initial guess
  v0  = uex;
  v0(2:N,2:M) = zeros(N-1,M-1);
  res = norm( residual(v0,f) );
   for p = 1:nCycle
     v1 = mg( v0,f,nSmooth,nLevel );
     v0 = v1;
     res = [res,norm(residual(v0,f))];
   end

  return
end

function v = mg( v0,f,nSmooth,nLevel )

  [N,M] = size(v0);
  vs = MGJacobi( v0,f,nSmooth );
  r  = residual( vs,f );
  rc = restrict( r );
  if( nLevel == 2 )
    ec = solveSystem( rc );
  else
    ec0 = zeros(size(rc));
    ec = mg( ec0,rc,nSmooth,nLevel-1 );
  end
  e = interpolate( ec,N,M );
  v = MGJacobi( vs+e,f,nSmooth );
  
  return
end

function vSmooth = MGJacobi( v,f,nSmooth )

  vSmooth = v;
  [N,M] = size( vSmooth );
  dx = 1./(N-1);
  dy = 1./(M-1);
  
  omega = 4/5;
  
  for k = 1:nSmooth
    r = residual( vSmooth,f );
    vSmooth = vSmooth-omega*r/(2./dx^2+2./dy^2);
  end
  
  return
end

function r = residual( v,f )
  [N,M] = size( v );
  dx = 1./(N-1);
  dy = 1./(M-1);
  
  r = zeros( N,M );
  for j = 2:N-1
    for k = 2:M-1
      vxx = (v(j+1,k)-2.0*v(j,k)+v(j-1,k))/dx^2;
      vyy = (v(j,k+1)-2.0*v(j,k)+v(j,k-1))/dy^2;
      r(j,k) = f(j,k)-vxx-vyy;
    end
  end
  
  return
end

function rc = restrict( rf )
  %% restrict to 2x coarser grid 
  [Nf,Mf] = size( rf );
  xf = linspace(0,1,Nf);
  yf = linspace(0,1,Mf);
    
  Nc = floor((Nf+1)/2);
  Mc = floor((Mf+1)/2);
  if( min(Nc,Mc) <= 2 )
    fprintf( 'coarse grid too small: too many multigrid levels' );
    stop
  end
  xc = linspace(0,1,Nc);
  yc = linspace(0,1,Mc);
  
  rc = (interp2( xf,yf',rf',xc,yc' ))';
  
  return
end

function v = solveSystem( f )
  %% solve the discrete system for u_xx+u_yy = f on unit square with 
  %%  zero dirichlet BCs
  
  [N,M] = size( f );
  dx = 1./(N-1);
  dy = 1./(M-1);
  cx = 1./dx^2;
  cy = 1./dy^2;
  d = -2.*(cx+cy);
  
  A = zeros( (N-2)*(M-2),(N-2)*(M-2) );
  b = zeros( (N-2)*(M-2),1 );
  
  %% form A and b
  for j = 1:N-2
    for k = 1:M-2
      ieq = getInd(j,k,N-2,M-2);
      A(ieq,ieq) = d;
      b(ieq,1)   = f(j+1,k+1);
      if( j > 1 )
        A(ieq,getInd(j-1,k,N-2,M-2)) = cx;
      end
      if( j < N-2 )
        A(ieq,getInd(j+1,k,N-2,M-2)) = cx;
      end
      
      if( k > 1 )
        A(ieq,getInd(j,k-1,N-2,M-2)) = cy;
      end
      if( k < M-2 )
        A(ieq,getInd(j,k+1,N-2,M-2)) = cy;
      end
    end
  end
  
  %% solve
  u = A\b;
  
  %% unpack
  v = zeros( N,M );
  for j = 1:N-2
    for k = 1:M-2
      v(j+1,k+1) = u(getInd(j,k,N-2,M-2));
    end
  end
  
  return
end

function ind = getInd( j,k,N,M )
  %% indexing function
  ind = (k-1)*N+j;
  %ind = (j-1)*M+k;
  return
end

function ef = interpolate( ec,Nf,Mf )
  % course to fine interpolation (using matlab interpolation)
  [Nc,Mc] = size( ec );
  xc = linspace(0,1,Nc);
  yc = linspace(0,1,Mc);
  
  xf = linspace(0,1,Nf);
  yf = linspace(0,1,Mf);

  ef = (interp2( xc,yc',ec',xf,yf' ))';
       
  return
end
