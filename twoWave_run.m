%% create results
 Nr = 160;
 Ns = 160;
 c = 1;
 iOption = 2; %% zero initial condition
 iTZ = 2; % doesn't in fact matter since not doing TZ
 sigma = .9;

 tf = 1.0000e-15;
   
 [x,y,u_0,~] = twoWave(Nr,Ns,c,sigma,tf,iOption,iTZ);
 tf = .1;
 [x,y,u_1,~] = twoWave( Nr,Ns,c,sigma,tf,iOption,iTZ );
 tf = .5;
 [x,y,u_5,~] = twoWave( Nr,Ns,c,sigma,tf,iOption,iTZ );
 tf = 1.5;
 [x,y,u_15,~] = twoWave( Nr,Ns,c,sigma,tf,iOption,iTZ );

 figure
 fs = 16;
 lineWidth = 2;
 ms = 16;
 set(gca,'FontSize',fs);
surf( x,y,u_0 );
 shading interp
 xlabel( 'x' );
 ylabel( 'y' );
 zlabel( 'z' );
 title( sprintf( 't=0' ) );
 plotName = sprintf('images/two_t0.eps');
 fprintf('Saving file=[%s]\n',plotName);
 print('-depsc2',plotName);

 figure
 fs = 16;
 lineWidth = 2;
 ms = 16;
 set(gca,'FontSize',fs);
 surf( x,y,u_1 );
 shading interp
 xlabel( 'x' );
 ylabel( 'y' );
 zlabel( 'z' );
 title( sprintf( 't=0.1' ) );
 plotName = sprintf('images/two_t1.eps');
 fprintf('Saving file=[%s]\n',plotName);
 print('-depsc2',plotName);

 figure
 fs = 16;
 lineWidth = 2;
 ms = 16;
 set(gca,'FontSize',fs);
 surf( x,y,u_5 );
 shading interp
 xlabel( 'x' );
 ylabel( 'y' );
 zlabel( 'z' );
 title( sprintf( 't=0.5' ) );
 plotName = sprintf('images/two_t5.eps');
 fprintf('Saving file=[%s]\n',plotName);
print('-depsc2',plotName);

 figure
 fs = 16;
 lineWidth = 2;
 ms = 16;
 set(gca,'FontSize',fs);
 surf( x,y,u_15 );
 shading interp
 xlabel('x');
 ylabel('y' );
 zlabel( 'z' );
 title( sprintf( 't=1.5' ) );
 plotName = sprintf('images/two_t15.eps');
 fprintf('Saving file=[%s]\n',plotName);
 print('-depsc2',plotName);


 figure
 fs = 16;
 lineWidth = 2;
 ms = 16;
 set(gca,'FontSize',fs);
 surf( x,y,u_5 );
 shading interp
 xlabel( 'x' );
 ylabel( 'y' );
 zlabel( 'z' );
 title( sprintf( 't=0.5' ) );
 plotName = sprintf('images/two_t5.eps');
 fprintf('Saving file=[%s]\n',plotName);
 print('-depsc2',plotName);

 figure
 fs = 16;
 lineWidth = 2;
 ms = 16;
 set(gca,'FontSize',fs);
 xb = x(1,:);
 yb = y(1,:);
 theta = atan(yb./xb);
 plot( theta,u_0(1,:),'r-','lineWidth',lineWidth,'MarkerSize',ms );
 hold on
 plot( theta,u_1(1,:),'g-','lineWidth',lineWidth,'MarkerSize',ms );
 plot( theta,u_5(1,:),'b-','lineWidth',lineWidth,'MarkerSize',ms );
 plot( theta,u_15(1,:),'k-','lineWidth',lineWidth,'MarkerSize',ms );
 hold off;
 xlabel( '\theta' );
 ylabel( 'u');
 legend( 't=0', 't=.1', 't=.5', 't=1.5' );
 plotName = sprintf('images/twoSlice.eps');
 fprintf('Saving file=[%s]\n',plotName);
 print('-depsc2',plotName);
