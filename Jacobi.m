function [v,resHist,refDiff]=Jacobi( u,f,dx,dy,N,M,nStep,ref,tol )
  v   = u;
  resHist = [];
  refDiff = [];
  for n = 0:nStep-1
    vold = v;
    denom = 2./(dx^2)+2./(dy^2);
    for j = 2:N
      for k = 2:M
        uxx = (vold(j+1,k)-2.0*vold(j,k)+vold(j-1,k))/(dx^2);
        uyy = (vold(j,k+1)-2.0*vold(j,k)+vold(j,k-1))/(dy^2);
        res(j,k) = f(j,k)-uxx-uyy;
        v(j,k) = vold(j,k)-res(j,k)/denom;
      end
    end
    figure( 1 ); surf( res ); shading interp; pause;
    resLoc       = max(max(abs(res)));
    refDiffLoc   = max(max(abs(v-ref)));
    resHist = [resHist,resLoc];
    refDiff = [refDiff,refDiffLoc];
    if( resLoc < tol )
      break
    end
  end
  return
end
    