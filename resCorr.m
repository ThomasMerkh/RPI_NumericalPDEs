% This program runs Jacobi, GaussSeidel, and SOR.
itCountJ = [];
itCountGS = [];
itCountSOR = [];
%runs = [5,10,20,40,80];
runs = 5:5:40;
nStep = 100000;
tol = 1e-10;

for N0 = runs
  N  = N0;
  M  = N0;
  
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

  uex = sinh(2.*x).*sin(y);
  f = 3.*uex;
  
  %% compute exact solution using Jacobi with lots of ieterations
  %%   ... and exact solution as starting guess
  u0  = uex;
  ref = uex;
  [v,res,refDiff] = Jacobi( u0,f,dx,dy,N,M,nStep,ref,tol );
  
  %% now compute the solution using Jacobi and a zero initial guess
  u0 = uex;
  u0(2:N,2:M) = zeros(N-1,M-1);
  ref = v;
  [vJ,resJ,refDiffJ] = Jacobi( u0,f,dx,dy,N,M,nStep,ref,tol );
  
  
  %% now compute the solution using Gauss Seidel and a zero initial guess
  u0 = uex;
  u0(2:N,2:M) = zeros(N-1,M-1);
  ref = v;
  [vGS,resGS,refDiffGS] = GaussSeidel( u0,f,dx,dy,N,M,nStep,ref,tol );
  
  %% now compute the solution using SOR and a zero initial guess
  u0 = uex;
  u0(2:N,2:M) = zeros(N-1,M-1);
  ref = v;
  [vSOR,resSOR,refDiffSOR] = SOR( u0,f,dx,dy,N,M,nStep,ref,tol );
  
  %% plots of residual
  figure(1)
  semilogy( 1:length(resJ),resJ,'bx' );
  hold on
  semilogy( 1:length(resGS),resGS,'rx' );
  semilogy( 1:length(resSOR),resSOR,'mx' );
  hold off
  legend( 'Jacobi', 'Gauss Seidel', 'SOR', 'Location', 'NorthEast' );
  xlabel( 'iteration' )
  ylabel( 'residual' )
  
  %% plots of error
  figure(2)
  semilogy( 1:length(refDiffJ), refDiffJ, 'bx' );
  hold on
  semilogy( 1:length(refDiffGS), refDiffGS, 'rx' );
  semilogy( 1:length(refDiffSOR), refDiffSOR, 'mx' );
  hold off
  legend( 'Jacobi', 'Gauss Seidel', 'SOR', 'Location', 'NorthEast' );
  xlabel( 'iteration' )
  ylabel( 'error' )
  
  itCountJ   = [itCountJ,length(resJ)];
  itCountGS  = [itCountGS,length(resGS)];
  itCountSOR = [itCountSOR,length(resSOR)];
  
  pause
  end

figure(3)
plot( runs,itCountJ,'bx-' );
hold on
plot( runs,itCountGS,'rx-' );
plot( runs,itCountSOR,'mx-' );
hold off
legend( 'Jacobi', 'Gauss Seidel', 'SOR', 'Location', 'NorthWest' );
xlabel( 'N' );
ylabel( 'iterations' );
