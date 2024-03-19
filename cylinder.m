lx     = 400;
ly     = 100;
obst_x = lx/5+1;
obst_y = ly/2+3;
obst_r = ly/10+1;
uMax   = 0.1;
Re     = 100;
nu     = uMax * 2.*obst_r / Re;
omega  = 1. / (3*nu+1./2.);
maxT   = 400000;
tPlot  = 50;
t  = [4/9  1/9 1/9 1/9 1/9  1/36 1/36 1/36 1/36];
cx = [  0    1   0  -1   0     1   -1   -1    1];
cy = [  0    0   1   0  -1     1    1   -1   -1];
opp = [ 1    4   5   2   3     8    9    6    7];
col = 2:(ly-1);
[y,x] = meshgrid(1:ly,1:lx);
obst = (x-obst_x).^2 + (y-obst_y).^2 <= obst_r.^2;
obst(:,[1,ly]) = 1;
bbRegion = find(obst);
L = ly-2; y_phys = y-1.5;
u = 4 * uMax / (L*L) * (y_phys.*L-y_phys.*y_phys);
v = zeros(lx,ly);
rho = 1;
for i=1:9
  cu = 3*(cx(i)*u+cy(i)*v);
  f0(i,:,:) = rho .* t(i) .* ( 1 + cu + 1/2*(cu.*cu) - 3/2*(u.^2+v.^2) );
end

for cycle = 1:maxT
  rho = sum(f0);
  u  = reshape ( (cx * reshape(f0,9,lx*ly)), 1,lx,ly) ./rho;
  v  = reshape ( (cy * reshape(f0,9,lx*ly)), 1,lx,ly) ./rho;
			 # MACROSCOPIC (DIRICHLET) BOUNDARY CONDITIONS
			 # 1let: Poiseuille profile
  y_phys = col-1.5;
  u(:,1,col) = 4 * uMax / (L*L) * (y_phys.*L-y_phys.*y_phys);
  v(:,1,col) = 0;
  rho(:,1,col) = 1 ./ (1-u(:,1,col)) .* ( sum(f0([1,3,5],1,col)) + 2*sum(f0([4,7,8],1,col)) );

				# Lxlet: Constant pressure
  rho(:,lx,col) = 1;
  u(:,lx,col) = -1 + 1 ./ (rho(:,lx,col)) .* ( sum(f0([1,3,5],lx,col)) + 2*sum(f0([2,6,9],lx,col)) );
  v(:,lx,col)  = 0;

		  # MICROSCOPIC BOUNDARY CONDITIONS: 1LET (Zou/He BC)
  f0(2,1,col) = f0(4,1,col) + 2/3*rho(:,1,col).*u(:,1,col);
  f0(6,1,col) = f0(8,1,col) + 1/2*(f0(5,1,col)-f0(3,1,col)) ...
		  + 1/2*rho(:,1,col).*v(:,1,col) ...
		  + 1/6*rho(:,1,col).*u(:,1,col);
  f0(9,1,col) = f0(7,1,col) + 1/2*(f0(3,1,col)-f0(5,1,col)) ...
		  - 1/2*rho(:,1,col).*v(:,1,col) ...
		  + 1/6*rho(:,1,col).*u(:,1,col);

		 # MICROSCOPIC BOUNDARY CONDITIONS: LXLET (Zou/He BC)
  f0(4,lx,col) = f0(2,lx,col) - 2/3*rho(:,lx,col).*u(:,lx,col);
  f0(8,lx,col) = f0(6,lx,col) + 1/2*(f0(3,lx,col)-f0(5,lx,col)) ...
		   - 1/2*rho(:,lx,col).*v(:,lx,col) ...
		   - 1/6*rho(:,lx,col).*u(:,lx,col);
  f0(7,lx,col) = f0(9,lx,col) + 1/2*(f0(5,lx,col)-f0(3,lx,col)) ...
		   + 1/2*rho(:,lx,col).*v(:,lx,col) ...
		   - 1/6*rho(:,lx,col).*u(:,lx,col);
  for i=1:9
    cu = 3*(cx(i)*u+cy(i)*v);
    fEq(i,:,:)  = rho .* t(i) .* ...
		  ( 1 + cu + 1/2*(cu.*cu)  - 3/2*(u.^2+v.^2) );
    f1(i,:,:) = f0(i,:,:) - omega .* (f0(i,:,:)-fEq(i,:,:));
  end
				# OBSTACLE (BOUNCE-BACK)
  for i=1:9
    f1(i,bbRegion) = f0(opp(i),bbRegion);
  end
  for i=1:9
    f0(i,:,:) = circshift(f1(i,:,:), [0,cx(i),cy(i)]);
  end
  if (mod(cycle,tPlot)==0)
    path = sprintf("cyl.%09d.raw", cycle);
    fid = fopen(path, "w");
    for field = {u, v, rho, rho}
      fwrite(fid, field{1}, "double");
    end
    fclose(fid);
    disp([var(u(:)), var(v(:)), var(rho(:))]);
  end
end
