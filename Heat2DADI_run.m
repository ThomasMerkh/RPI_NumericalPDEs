%%%%%
%% quadratic TZ, constant D
%%%%%
tic(); [x,y,u,ue] = Heat2DADI( 10,12,10,2,1,1 ); toc()
fprintf( 'max error: %e\n', max(max(abs(u-ue))) );
surf( x,y,(u-ue)' );
xlabel( 'x' );
ylabel( 'y' );

%%%%%
%% quadratic TZ,  linear D
%%%%%
tic(); [x,y,u,ue] = Heat2DADI( 10,12,10,2,1,2 ); toc()
fprintf( 'max error: %e\n', max(max(abs(u-ue))) );
surf( x,y,(u-ue)' );
xlabel( 'x' );
ylabel( 'y' );

%%%%%
%% quadratic TZ,  trig D
%%%%%
tic(); [x,y,u,ue] = Heat2DADI( 10,12,10,2,1,3 ); toc()
fprintf( 'max error: %e\n', max(max(abs(u-ue))) );
surf( x,y,(u-ue)' );
xlabel( 'x' );
ylabel( 'y' );
uiwait;

tic(); [x,y,u,ue] = Heat2DADI( 20,24,20,2,1,3 ); toc()
fprintf( 'max error: %e\n', max(max(abs(u-ue))) );
surf( x,y,(u-ue)' );
xlabel( 'x' );
ylabel( 'y' );
uiwait;

surf( x,y,u' );

%%%%%
%% trig TZ,  trig D
%%%%%
tic(); [x,y,u,ue] = Heat2DADI( 10,12,10,2,3,3 ); toc()
fprintf( 'max error: %e\n', max(max(abs(u-ue))) );
surf( x,y,(u-ue)' );
xlabel( 'x' );
ylabel( 'y' );
uiwait;

tic(); [x,y,u,ue] = Heat2DADI( 20,20,20,2,3,3 ); toc()
fprintf( 'max error: %e\n', max(max(abs(u-ue))) );
surf( x,y,(u-ue)' );
xlabel( 'x' );
ylabel( 'y' );
uiwait;

tic(); [x,y,u,ue] = Heat2DADI( 40,48,40,2,3,1 ); toc()
fprintf( 'max error: %e\n', max(max(abs(u-ue))) );
surf( x,y,(u-ue)' );
xlabel( 'x' );
ylabel( 'y' );
uiwait;
%%
tic(); [x,y,u,ue] = Heat2DADI( 80,96,80,2,3,3 ); toc()
fprintf( 'max error: %e\n', max(max(abs(u-ue))) );
surf( x,y,(u-ue)' );
xlabel( 'x' );
ylabel( 'y' );
uiwait;

surf( x,y,u' );

