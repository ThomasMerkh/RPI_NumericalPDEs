function [v,resHist,refDiff]=SOR( u,f,dx,dy,N,M,nStep,ref,tol )
% successive over relaxation

  denom = 2./(dx^2)+2./(dy^2);
  
  %% optimal omega
  sigma = 2.*(cos(pi/N)/dx^2+cos(pi/M)/dy^2)/denom;
  omega = 2./(1.+sqrt(1.-sigma^2))
  
  %% guess at omega
  %omega = 1.9;
  
  v   = u;
  resHist = [];
  refDiff = [];
  for n = 0:nStep-1
    for j = 2:N
      for k = 2:M
        uxx = (v(j+1,k)-2.0*v(j,k)+v(j-1,k))/(dx^2);
        uyy = (v(j,k+1)-2.0*v(j,k)+v(j,k-1))/(dy^2);
        res(j,k) = f(j,k)-uxx-uyy;
        v(j,k) = v(j,k)-omega*res(j,k)/denom;
      end
    end
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