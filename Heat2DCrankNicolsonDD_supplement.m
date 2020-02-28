%% quadratic TZ
tic(); [x,y,u,ue] = Heat2DCrankNicolsonDD( 10,12,10,2,1,1 ); toc()
max(max(abs(u-ue)))


%% nonpoly TZ
tic(); [x,y,u,ue] = Heat2DCrankNicolsonDD( 10,12,10,2,1,2 ); toc()
fprintf( 'max error: %e\n', max(max(abs(u-ue))) );
surf( x,y,(u-ue)' );
xlabel( 'x' );
ylabel( 'y' );
uiwait;

tic(); [x,y,u,ue] = Heat2DCrankNicolsonDD( 20,24,20,2,1,2 ); toc()
fprintf( 'max error: %e\n', max(max(abs(u-ue))) );
surf( x,y,(u-ue)' );
xlabel( 'x' );
ylabel( 'y' );
uiwait;

tic(); [x,y,u,ue] = Heat2DCrankNicolsonDD( 40,48,40,2,1,2 ); toc()
fprintf( 'max error: %e\n', max(max(abs(u-ue))) );
surf( x,y,(u-ue)'*1e7 );
xlabel( 'x' );
ylabel( 'y' );
uiwait;

tic(); [x,y,u,ue] = Heat2DCrankNicolsonDD( 80,96,80,2,1,2 ); toc()
fprintf( 'max error: %e\n', max(max(abs(u-ue))) );
surf( x,y,(u-ue)'*1e7 );
xlabel( 'x' );
ylabel( 'y' );
