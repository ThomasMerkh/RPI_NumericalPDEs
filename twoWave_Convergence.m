%%%%%
 %% convergence study of two wave
 %%%%%

 tf = 1.;
 c = 0.9;
 sigma = 0.9;
 iOption = 1;
 iTZ = 1;
 N0 = 20;

 m = 4;
 err = zeros(m,1);
 h = zeros(m,1);
 for j = 0:m-1
N = N0*2^j;
 [x,t,u,err(j+1)] = twoWave( N,N,c,sigma,tf,iOption,iTZ );
 h(j+1) = x(2)-x(1);
 end

 figure
 fs = 16;
 lineWidth = 2;
 ms = 16;
 set(gca,’FontSize’,fs);
 loglog( h,err,’rx’,’lineWidth’,lineWidth,’MarkerSize’,ms );
 hold on
 loglog( h,1e0*h.^2,’k-’,’lineWidth’,lineWidth,’MarkerSize’,ms );
 xlabel( ’h’ );
 ylabel( ’max error’ );
 legend( ’numerics’, ’h^2 ref’,’Location’,’NorthWest’ );

 title( sprintf( ’convergence study’ ) );
 plotName = sprintf(’images/twoConv.eps’);
 fprintf(’Saving file=[%s]\n’,plotName);
 print(’-depsc2’,plotName);
