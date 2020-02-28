
t2 = [];
t3 = [];
t4 = [];

nSmooth = 3;
nCycle = 10;
N0 = [40:10:100];

for N = N0
  M = N;
  tic();
  [~,res2] = runmg( N,M,nSmooth,2,nCycle );
  t2 = [t2,toc()];

  tic();
  [~,res3] = runmg( N,M,nSmooth,3,nCycle );
  t3 = [t3,toc()];

  tic();
  [~,res4] = runmg( N,M,nSmooth,4,nCycle );
  t4 = [t4,toc()];



  figure()
  semilogy( res2,'bx' );
  hold on
  semilogy( res3,'rs' );
  semilogy( res4,'md' )
  hold off
  legend( '2-level', '3-level', '4-level', 'Location','NorthEast' );
  xlabel( 'v-cycle' );
  ylabel( 'residual' );
  title( sprintf('N=%i', N) );
end

figure
plot( N0,t2,'bx' );
hold on
plot( N0,t3,'rs' );
plot( N0,t4,'md' );
hold off
legend( '2-level', '3-level', '4-level', 'Location','NorthEast' );
xlabel( 'N' );
ylabel( 'time' );

