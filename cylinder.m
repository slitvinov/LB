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
in  = 1;
out = lx;

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
			 # Inlet: Poiseuille profile
  y_phys = col-1.5;
  u(:,in,col) = 4 * uMax / (L*L) * (y_phys.*L-y_phys.*y_phys);
  v(:,in,col) = 0;
  rho(:,in,col) = 1 ./ (1-u(:,in,col)) .* ( sum(f0([1,3,5],in,col)) + 2*sum(f0([4,7,8],in,col)) );

				# Outlet: Constant pressure
  rho(:,out,col) = 1;
  u(:,out,col) = -1 + 1 ./ (rho(:,out,col)) .* ( sum(f0([1,3,5],out,col)) + 2*sum(f0([2,6,9],out,col)) );
  v(:,out,col)  = 0;

		  # MICROSCOPIC BOUNDARY CONDITIONS: INLET (Zou/He BC)
  f0(2,in,col) = f0(4,in,col) + 2/3*rho(:,in,col).*u(:,in,col);
  f0(6,in,col) = f0(8,in,col) + 1/2*(f0(5,in,col)-f0(3,in,col)) ...
		  + 1/2*rho(:,in,col).*v(:,in,col) ...
		  + 1/6*rho(:,in,col).*u(:,in,col);
  f0(9,in,col) = f0(7,in,col) + 1/2*(f0(3,in,col)-f0(5,in,col)) ...
		  - 1/2*rho(:,in,col).*v(:,in,col) ...
		  + 1/6*rho(:,in,col).*u(:,in,col);

		 # MICROSCOPIC BOUNDARY CONDITIONS: OUTLET (Zou/He BC)
  f0(4,out,col) = f0(2,out,col) - 2/3*rho(:,out,col).*u(:,out,col);
  f0(8,out,col) = f0(6,out,col) + 1/2*(f0(3,out,col)-f0(5,out,col)) ...
		   - 1/2*rho(:,out,col).*v(:,out,col) ...
		   - 1/6*rho(:,out,col).*u(:,out,col);
  f0(7,out,col) = f0(9,out,col) + 1/2*(f0(5,out,col)-f0(3,out,col)) ...
		   + 1/2*rho(:,out,col).*v(:,out,col) ...
		   - 1/6*rho(:,out,col).*u(:,out,col);
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
