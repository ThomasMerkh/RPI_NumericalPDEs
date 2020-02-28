%% set some parameters
res   = 4;
r     = .4;
N     = 20*res;
x     = linspace( 0,1,N);
dx    = 1/(N-1);
nStep = 50*res^2;

%% calculate nu*dt 
nudt = r*dx^2;

%% allocate space
u   = zeros(nStep,N);
uex = zeros(nStep,N);

%% set initial condition
for j = 1:N
  u(1,j)   = x(j)+sin(pi*x(j))+sin(3.*pi*x(j));
  uex(1,j) = u(1,j);
end

%% time stepping loop
for n = 1:nStep-1

  %% boundary conditions
  u(n+1,1) = 0;
  u(n+1,N) = 1;

  %% space loop
  for j = 2:N-1
    %% FD formula
    u(n+1,j) = u(n,j)+r*(u(n,j+1)-2.*u(n,j)+u(n,j-1));
  end

  
  %% exact soln
  for j = 1:N
    nut = (n+1)*nudt;
    uex(n+1,j) = x(j)+exp(-nut*pi^2)*sin(pi*x(j))+exp(-nut*9*pi^2)*sin(3*pi*x(j));
  end
  
end


%%% do some output
fs = 16;

figure
set(gca,'FontSize',fs);
surf( x,[1:nStep],u );
view([-60,15]);
shading interp
xlabel( 'x' )
ylabel( 'time step' );
title( sprintf( 'N=%i, r=%e, approximation',N,r ) );
plotName = sprintf('v_N%irp%i',N,round(10*r));
fprintf('Saving file=[%s]\n',plotName);
print('-depsc2',plotName);

figure
set(gca,'FontSize',fs);
surf( x,[1:nStep],u-uex );
view([-125,50]);
shading interp
xlabel( 'x' )
ylabel( 'time step' );
title( sprintf( 'N=%i, r=%e, error',N,r ) );
plotName = sprintf('e_N%irp%i',N,round(10*r));
fprintf('Saving file=[%s]\n',plotName);
print('-depsc2',plotName);

fprintf( 'maximum error: %e\n', max(max(abs(u-uex))) );