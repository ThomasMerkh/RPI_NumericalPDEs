%% quadratic TZ
tic(); [x,y,u,ue] = Heat2DCrankNicolsonGhost( 10,12,10,2,1,1 ); toc()
fprintf( 'max error: %e\n', max(max(abs(u-ue))) );
surf( x,y,(u-ue)' );
xlabel( 'x' );
ylabel( 'y' );


%% linear TZ
tic(); [x,y,u,ue] = Heat2DCrankNicolsonGhost( 10,12,10,2,1,2 ); toc()
fprintf( 'max error: %e\n', max(max(abs(u-ue))) );


%% nonpoly TZ
tic(); [x,y,u,ue] = Heat2DCrankNicolsonGhost( 10,12,10,2,1,3 ); toc()
fprintf( 'max error: %e\n', max(max(abs(u-ue))) );
surf( x,y,(u-ue)' );
xlabel( 'x' );
ylabel( 'y' );
tic(); [x,y,u,ue] = Heat2DCrankNicolsonGhost( 20,24,20,2,1,3 ); toc()
fprintf( 'max error: %e\n', max(max(abs(u-ue))) );
surf( x,y,(u-ue)' );
xlabel( 'x' );
ylabel( 'y' );

tic(); [x,y,u,ue] = Heat2DCrankNicolsonGhost( 40,48,40,2,1,3 ); toc()
fprintf( 'max error: %e\n', max(max(abs(u-ue))) );
surf( x,y,(u-ue)'*1e7 );
xlabel( 'x' );
ylabel( 'y' );

tic(); [x,y,u,ue] = Heat2DCrankNicolsonGhost( 80,96,80,2,1,3 ); toc()
fprintf( 'max error: %e\n', max(max(abs(u-ue))) );
surf( x,y,(u-ue)'*1e7 );
xlabel( 'x' );
ylabel( 'y' );
