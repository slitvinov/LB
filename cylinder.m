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

				# D2Q9 LATTICE CONSTANTS
t  = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
cx = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1];
cy = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1];
opp = [ 1,   4,  5,  2,  3,    8,   9,   6,   7];
col = [2:(ly-1)];
[y,x] = meshgrid(1:ly,1:lx);
obst = ...                   # Location of cylinder
(x-obst_x).^2 + (y-obst_y).^2 <= obst_r.^2;
obst(:,[1,ly]) = 1;    # Location of top/bottom boundary
bbRegion = find(obst); # Boolean mask for bounce-back cells

		# INITIAL CONDITION: Poiseuille profile at equilibrium
L = ly-2; y_phys = y-1.5;
ux = 4 * uMax / (L*L) * (y_phys.*L-y_phys.*y_phys);
uy = zeros(lx,ly);
rho = 1;
for i=1:9
  cu = 3*(cx(i)*ux+cy(i)*uy);
  fIn(i,:,:) = rho .* t(i) .* ( 1 + cu + 1/2*(cu.*cu) - 3/2*(ux.^2+uy.^2) );
end

for cycle = 1:maxT
  rho = sum(fIn);
  ux  = reshape ( (cx * reshape(fIn,9,lx*ly)), 1,lx,ly) ./rho;
  uy  = reshape ( (cy * reshape(fIn,9,lx*ly)), 1,lx,ly) ./rho;

			 # MACROSCOPIC (DIRICHLET) BOUNDARY CONDITIONS
			 # Inlet: Poiseuille profile
  y_phys = col-1.5;
  ux(:,1,col) = 4 * uMax / (L*L) * (y_phys.*L-y_phys.*y_phys);
  uy(:,1,col) = 0;
  rho(:,1,col) = 1 ./ (1-ux(:,1,col)) .* ( sum(fIn([1,3,5],1,col)) + 2*sum(fIn([4,7,8],1,col)) );

				# Outlet: Constant pressure
  rho(:,lx,col) = 1;
  ux(:,lx,col) = -1 + 1 ./ (rho(:,lx,col)) .* ( sum(fIn([1,3,5],lx,col)) + 2*sum(fIn([2,6,9],lx,col)) );
  uy(:,lx,col)  = 0;

		  # MICROSCOPIC BOUNDARY CONDITIONS: INLET (Zou/He BC)
  fIn(2,1,col) = fIn(4,1,col) + 2/3*rho(:,1,col).*ux(:,1,col);
  fIn(6,1,col) = fIn(8,1,col) + 1/2*(fIn(5,1,col)-fIn(3,1,col)) ...
		  + 1/2*rho(:,1,col).*uy(:,1,col) ...
		  + 1/6*rho(:,1,col).*ux(:,1,col);
  fIn(9,1,col) = fIn(7,1,col) + 1/2*(fIn(3,1,col)-fIn(5,1,col)) ...
		  - 1/2*rho(:,1,col).*uy(:,1,col) ...
		  + 1/6*rho(:,1,col).*ux(:,1,col);

		 # MICROSCOPIC BOUNDARY CONDITIONS: OUTLET (Zou/He BC)
  fIn(4,lx,col) = fIn(2,lx,col) - 2/3*rho(:,lx,col).*ux(:,lx,col);
  fIn(8,lx,col) = fIn(6,lx,col) + 1/2*(fIn(3,lx,col)-fIn(5,lx,col)) ...
		   - 1/2*rho(:,lx,col).*uy(:,lx,col) ...
		   - 1/6*rho(:,lx,col).*ux(:,lx,col);
  fIn(7,lx,col) = fIn(9,lx,col) + 1/2*(fIn(5,lx,col)-fIn(3,lx,col)) ...
		   + 1/2*rho(:,lx,col).*uy(:,lx,col) ...
		   - 1/6*rho(:,lx,col).*ux(:,lx,col);

				# COLLISION STEP
  for i=1:9
    cu = 3*(cx(i)*ux+cy(i)*uy);
    fEq(i,:,:)  = rho .* t(i) .* ...
		  ( 1 + cu + 1/2*(cu.*cu)  - 3/2*(ux.^2+uy.^2) );
    fOut(i,:,:) = fIn(i,:,:) - omega .* (fIn(i,:,:)-fEq(i,:,:));
  end

				# OBSTACLE (BOUNCE-BACK)
  for i=1:9
    fOut(i,bbRegion) = fIn(opp(i),bbRegion);
  end

				# STREAMING STEP
  for i=1:9
    fIn(i,:,:) = circshift(fOut(i,:,:), [0,cx(i),cy(i)]);
  end
  if (mod(cycle,tPlot)==0)
    path = sprintf("cyl.%09d.raw", cycle);
    fid = fopen(path, "w");
    for field = {ux, uy, rho, rho}
      fwrite(fid, field{1}, "double");      
    end
    fclose(fid);
  end
end
