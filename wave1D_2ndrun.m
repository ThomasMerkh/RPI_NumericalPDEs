%%%%%
%% case 1a
%%%%%
[x,t,u] = wave1D_2nd( 100,300,3,.95,1 );
surf( x,t,u' )
shading interp
xlabel( 'x' );
ylabel( 't' );

%%%%%
%% case 1b
%%%%%
[x,t,u] = wave1D_2nd( 100,300,3,1.001,1 );
surf( x,t,u' )
shading interp
xlabel( 'x' );
ylabel( 't' );

%%%%%
%% case 2
%%%%%
[x,t,u] = wave1D_2nd( 100,1000,10,.95,2 );
surf( x,t,u' )
shading interp
xlabel( 'x' );
ylabel( 't' );

%%%%%
%% case 3
%%%%%
[x,t,u] = wave1D_2nd( 100,300,3,.95,3 );
surf( x,t,u' )
shading interp
xlabel( 'x' );
ylabel( 't' );

%%%%%
%% convergence study
%%%%%
tf = 1.;
c = .9;
N = 20;
[x,t,u] = wave1D_2nd( N,N,tf,c,4 );
fprintf( 'N=%i, e=%e\n', N, max(abs(u(:,N+1)-sin(5*(x'-c*tf)))));

N = 40;
[x,t,u] = wave1D_2nd( N,N,tf,c,4 );
fprintf( 'N=%i, e=%e\n', N, max(abs(u(:,N+1)-sin(5*(x'-c*tf)))));

N = 80;
[x,t,u] = wave1D_2nd( N,N,tf,c,4 );
fprintf( 'N=%i, e=%e\n', N, max(abs(u(:,N+1)-sin(5*(x'-c*tf)))));

N = 160;
[x,t,u] = wave1D_2nd( N,N,tf,c,4 );
fprintf( 'N=%i, e=%e\n', N, max(abs(u(:,N+1)-sin(5*(x'-c*tf)))));

N = 320;
[x,t,u] = wave1D_2nd( N,N,tf,c,4 );
fprintf( 'N=%i, e=%e\n', N, max(abs(u(:,N+1)-sin(5*(x'-c*tf)))));
