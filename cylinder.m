nx     = 400;
ny     = 100;
obst_x = nx/5+1;
obst_y = ny/2+3;
obst_r = ny/10+1;
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
col = 2:(ny-1);
[y,x] = meshgrid(1:ny,1:nx);
obst = (x-obst_x).^2 + (y-obst_y).^2 <= obst_r.^2;
obst(:,[1,ny]) = 1;
bb = find(obst);
L = ny-2; y_phys = y-1.5;
u = 4 * uMax / (L*L) * (y_phys.*L-y_phys.*y_phys);
v = zeros(nx,ny);
rho = 1;
for i=1:9
  cu = 3*(cx(i)*u+cy(i)*v);
  f(i,:,:) = rho .* t(i) .* ( 1 + cu + 1/2*(cu.*cu) - 3/2*(u.^2+v.^2) );
end

for cycle = 0 : maxT - 1
  rho = sum(f);
  u  = reshape ( (cx * reshape(f,9,nx*ny)), 1,nx,ny) ./rho;
  v  = reshape ( (cy * reshape(f,9,nx*ny)), 1,nx,ny) ./rho;
			 # MACROSCOPIC (DIRICHLET) BOUNDARY CONDITIONS
			 # inlet: Poiseuille profile
  y_phys = col-1.5;
  u(:,1,col) = 4 * uMax / (L*L) * (y_phys.*L-y_phys.*y_phys);
  v(:,1,col) = 0;
  rho(:,1,col) = 1 ./ (1-u(:,1,col)) .* ( sum(f([1,3,5],1,col)) + 2*sum(f([4,7,8],1,col)) );

				# outlet: Constant pressure
  rho(:,nx,col) = 1;
  u(:,nx,col) = -1 + 1 ./ (rho(:,nx,col)) .* ( sum(f([1,3,5],nx,col)) + 2*sum(f([2,6,9],nx,col)) );
  v(:,nx,col)  = 0;

		  # MICROSCOPIC BOUNDARY CONDITIONS: inlet (Zou/He BC)
  f(2,1,col) = f(4,1,col) + 2/3*rho(:,1,col).*u(:,1,col);
  f(6,1,col) = f(8,1,col) + 1/2*(f(5,1,col)-f(3,1,col)) ...
		  + 1/2*rho(:,1,col).*v(:,1,col) ...
		  + 1/6*rho(:,1,col).*u(:,1,col);
  f(9,1,col) = f(7,1,col) + 1/2*(f(3,1,col)-f(5,1,col)) ...
		  - 1/2*rho(:,1,col).*v(:,1,col) ...
		  + 1/6*rho(:,1,col).*u(:,1,col);

		 # MICROSCOPIC BOUNDARY CONDITIONS: outlet (Zou/He BC)
  f(4,nx,col) = f(2,nx,col) - 2/3*rho(:,nx,col).*u(:,nx,col);
  f(8,nx,col) = f(6,nx,col) + 1/2*(f(3,nx,col)-f(5,nx,col)) ...
		   - 1/2*rho(:,nx,col).*v(:,nx,col) ...
		   - 1/6*rho(:,nx,col).*u(:,nx,col);
  f(7,nx,col) = f(9,nx,col) + 1/2*(f(5,nx,col)-f(3,nx,col)) ...
		   + 1/2*rho(:,nx,col).*v(:,nx,col) ...
		   - 1/6*rho(:,nx,col).*u(:,nx,col);
  for i=1:9
    cu = 3*(cx(i)*u+cy(i)*v);
    feq(i,:,:)  = rho .* t(i) .* ...
		  ( 1 + cu + 1/2*(cu.*cu)  - 3/2*(u.^2+v.^2) );
    f1(i,:,:) = f(i,:,:) - omega .* (f(i,:,:)-feq(i,:,:));
  end
				# OBSTACLE (BOUNCE-BACK)
  for i=1:9
    f1(i,bb) = f(opp(i),bb);
  end

  for i=1:9
    f(i,:,:) = circshift(f1(i,:,:), [0,cx(i),cy(i)]);
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
