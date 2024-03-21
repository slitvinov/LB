c     D2Q9 lattice, BGK version
c     0= rest particles, 1-4, nearest-neigh(nn), 5-8(nnn)
      implicit double precision(a-h, o-z)
      include 'bgk2.par'
      call input
      call inithydro
      call equil
      call initpop
      iconf=0
      do 10 istep = 1,nsteps
c     periodic boundary conditions
css   call pbc
c     mixed boundary conditions c(Poiseuille flow)
         call mbc
         call move
         call hydrovar
         call equil
         call colli
         if(iforce) then
            call force()
         endif
c     Obstacle ?
         if(iobst.eq.1) call obst
c     0 - dim diagnostic
         if(mod(istep,ndiag) .eq. 0) then
            call diag0D(istep)
         endif
c     movie
         if(mod(istep, 100) .eq. 0) then
            call movie(istep)
         endif
c     2d configs
         if(mod(istep, nout) .eq. 0) then
            call config(istep, iconf)
            iconf=iconf+1
         endif
 10   continue
      end

      subroutine Input
      implicit double precision(a-h, o-z)
      include 'bgk2.par'
      print *, 'Number of steps'
      read(5 ,*) nsteps
      print *, 'Number of steps between printing profile'
      read(5, *) nout
      print *, 'Number of steps between performing diagnostics'
      read(5 ,*) ndiag
      print *, 'Relaxation frequency omega'
      read(5, *) omega
      print *, 'Applied force(T or F)) ?'
      read(5 ,*) iforce
      print *, 'Initial density and velocity for the Poiseuille force'
      read(5, *) rho0, u0, v0
      print *, 'Final velocity for the Poise force'
      read(5 ,*) uf
      print *, 'Linear obstacle(T or F)?'
      read(5 ,*) iobst
      print *, 'Obstacle height?'
      read(5 ,*) nobst
      print *, 'Obstacle id', iobst
      print *, 'Length of the obstacle(multple of 2)', nobst
      print *, 'File for output : 5 chars'
      read(5, '(A)') fileout
      open(10, file=fileout      //   '.uy')
      open(11, file=fileout      //   '.vy')
      open(12, file=fileout      //   '.rhoy')
      open(14, file=fileout      //   '.uvx')
      open(16, file=fileout      //   '.pop')
      open(50, file=fileout      //   '.probe')
      open(35, file = 'porous')

      print * ,'*****************************************'
      print *, 'Lattice BGK model, 2D with 9 velocities'
      print * ,'*****************************************'
      print *, 'Number of cells :' , nx, '*', ny
      print *, 'Nsteps : ', nsteps
      print *, 'Relaxation frequency :', omega
      print *, 'Initial velocity for this Poiseuille force :', u0
      print *, 'Initial density :', rho0
      print *, 'Applied force :', iforce
      if(iobst.eq.1) then
         print *, 'Linear Obstacle with length :', nobst
      endif
      write(6, '(A)') 'Output file :', fileout
c     lattice weights
      w0 = 4.0d0/9.0d0
      w1 = 1.0/9.0d0
      w2 = 1.0/36.0d0
c     sound - speed and related constants
      cs2 = 1.0d0 / 3.0d0
      cs22 = 2.0d0 * cs2
      cssq = 2.0d0 / 9.0d0
c     viscosity and nominal Reynolds
      visc =(1.0d0 / omega - 0.5d0) * cs2
      rey = u0*ny / visc
      print *, 'Viscosity and nominal Reynolds:', visc, rey
      if(visc .lt. 0) stop 'OMEGA OUT of(0, 2) interval!!'

c     Applied force(based on Stokes problem)
      fpois = 8.0d0 * visc * uf / dfloat(ny) / dfloat(ny)
c     # of biased populations
      fpois = rho0 * fpois / 6.
      print *, 'Intensity of the applied force', fpois
      end

      subroutine Inithydro
      implicit double precision(a-h, o-z)
      include 'bgk2.par'
      write(6, *) 'u0', u0
      do j = 0, ny+1
         do i = 0, nx+1
            rho(i,j) = rho0
            u(i,j)   = u0
            v(i,j)   = v0
         enddo
      enddo
      end

      subroutine Initpop
      implicit double precision(a- h, o-z)
      include 'bgk2.par'
      iseed = 15391
c     random amplitude
      ramp = 0.01
      do j = 0, ny+1
         do i = 0, nx+1
            do ip = 0, npop - 1
               rr = 2.* rand(iseed) - 1.
               f(ip ,i,j) =(1. + ramp * rr)* rho0 / float(npop)
               iseed = iseed + 1
            end do
         end do
      end do
      end

      subroutine Move
      implicit double precision(a-h, o-z)
      include 'bgk2.par'
      do j = ny, 1, -1
         do i = 1, nx
            f(2, i, j) = f(2, i, j - 1)
            f(6, i, j) = f(6, i + 1, j -1)
         enddo
      enddo
      do j = ny,1,-1
         do i = nx, 1, -1
            f(1,i,j) = f(1, i - 1, j)
            f(5,i,j) = f(5, i - 1, j - 1)
         enddo
      enddo
      do j = 1, ny
         do i = nx, 1, -1
            f(4, i, j) = f(4, i, j + 1)
            f(8, i, j) = f(8, i - 1, j +1)
         enddo
      enddo
      do j = 1, ny
         do i = 1, nx
            f(3, i, j) = f(3, i + 1, j)
            f(7, i, j) = f(7, i + 1, j + 1)
         enddo
      enddo
      end

      subroutine Hydrovar
      implicit double precision(a-h, o-z)
      include 'bgk2.par'
c     hydro variables
      do j = 1, ny
         do i = 1, nx
            rho(i,j)= f(1, i, j) + f(2, i, j) + f(3, i, j)
     $           + f(4, i, j) + f(5, i, j) + f(6, i, j)
     $           + f(7, i, j) + f(8, i, j) + f(0, i, j)
            rhoi = 1. / rho(i, j)
            u(i, j) =(f(1, i, j) - f(3, i, j) + f(5, i, j) -
     $           f(6, i, j) - f(7, i, j) + f(8, i, j)) * rhoi
            v(i, j) =(f(5, i, j) + f(2, i, j)+ f(6, i, j)
     $           - f(7, i, j) - f(4, i, j) - f(8, i, j)) * rhoi
         enddo
      enddo

      end

      subroutine Equil
      implicit double precision(a-h, o-z)
      include 'bgk2.par'
c     equils are written explicitly to avoid multplications by zero
      do j = 0, ny+1
         do i = 0, nx+1
            ul = u(i,j) / cs2
            vl = v(i,j) / cs2
            uv = ul * vl
            usq = u(i,j)* u(i,j)
            vsq = v(i,j)* v(i,j)
            sumsq =(usq+vsq) / cs22
            sumsq2 = sumsq *(1.0d0 - cs2) / cs2
            u2 = usq / cssq
            v2 = vsq / cssq
            feq(0, i, j) = w0 *(1.0d0 - sumsq)

            feq(1, i, j) = w1 *(1.0d0 - sumsq + u2 + ul)
            feq(2, i, j) = w1 *(1.0d0 - sumsq + v2 + vl)
            feq(3, i, j) = w1 *(1.0d0 - sumsq + u2 - ul)
            feq(4, i, j) = w1 *(1.0d0 - sumsq + v2 - vl)

            feq(5, i, j) = w2 *(1.0d0 + sumsq2 + ul + vl + uv)
            feq(6, i, j) = w2 *(1.0d0 + sumsq2 - ul + vl - uv)
            feq(7, i, j) = w2 *(1.0d0 + sumsq2 - ul - vl + uv)
            feq(8, i, j) = w2 *(1.0d0 + sumsq2 + ul - vl - uv)
         enddo
      enddo
      end

      subroutine Colli
      implicit double precision(a-h, o-z)
      include 'bgk2.par'
      do k = 0, npop - 1
         do j = 1, ny
            do i = 1, nx
               f(k, i, j) = f(k, i, j) *(1.0d0 - omega) +
     $              omega* feq(k, i, j)
            enddo
         enddo
      enddo
      end

      subroutine Force()
      implicit double precision(a-h, o-z)
      include 'bgk2.par'
      do j = 1, ny
         do i = 1, nx
            f(1, i, j) = f(1, i, j) + fpois
            f(5, i, j) = f(5, i, j) + fpois
            f(8, i, j) = f(8, i, j) + fpois

            f(3, i, j) = f(3, i, j) - fpois
            f(6, i, j) = f(6, i, j) - fpois
            f(7, i, j) = f(7, i, j) - fpois
         enddo
      enddo
      end

      subroutine Pbc
      implicit double precision(a-h, o-z)
      include 'bgk2.par'
c     EAST
      do j = 1, ny
         f(1, 0, j) = f(1, nx, j)
         f(5, 0, j) = f(5, nx, j)
         f(8, 0, j) = f(8, nx, j)
      enddo
c     WEST
      do j = 1, ny
         f(3, nx +1, j) = f(3, 1, j)
         f(6, nx +1, j) = f(6, 1, j)
         f(7, nx +1, j) = f(7, 1, j)
      enddo
c     NORTH
      do i = 1, nx
         f(2, i, 0) = f(2, i, ny)
         f(5, i, 0) = f(5, i, ny)
         f(6, i, 0) = f(6, i, ny)
      enddo
c     SOUTH
      do i = 1, nx
         f(4, i, ny+1) = f(4, i, 1)
         f(7, i, ny+1) = f(7, i, 1)
         f(8, i, ny+1) = f(8, i, 1)
      enddo
      end
      subroutine Mbc
      implicit double precision(a-h, o-z)
      include 'bgk2.par'
c     WEST inlet
      do j = 1, ny
         f(1, 0, j) = f(1, nx, j)
         f(5, 0, j) = f(5, nx, j)
         f(8, 0, j) = f(8, nx, j)
      enddo
c     EAST outlet
      do j = 1, ny
         f(3, nx +1, j) = f(3, 1, j)
         f(6, nx +1, j) = f(6, 1, j)
         f(7, nx +1, j) = f(7, 1, j)
      enddo
c     NORTH solid
      do i = 1, nx
         f(4, i, ny+1) = f(2, i, ny)
         f(8, i, ny+1) = f(6, i +1, ny)
         f(7, i, ny+1) = f(5, i - 1, ny)
      enddo
c     SOUTH solid
      do i = 1, nx
         f(2, i ,0) = f(4, i ,1)
         f(6, i, 0) = f(8, i - 1 ,1)
         f(5, i, 0) = f(7, i +1 ,1)
      enddo
c     corners bounce-back
      f(8, 0, ny + 1) = f(6, 1, ny)
      f(5, 0, 0) = f(7, 1, 1)
      f(7, nx + 1, ny + 1) = f(5, nx, ny)
      f(6, nx + 1, 0) = f(8, nx, 1)
      end

      subroutine Obst
      implicit double precision(a-h, o-z)
      include 'bgk2.par'
      i         = nx / 4
      jbot = ny / 2 - nobst / 2
      jtop = ny / 2 + nobst / 2 + 1
      do j = ny / 2 - nobst / 2, ny / 2 + nobst / 2 + 1
         f(1 ,i,j) = f(3, i + 1, j)
         f(5 ,i,j) = f(7, i + 1, j +1)
         f(8 ,i,j) = f(6, i + 1, j - 1)
         f(3 ,i,j) = f(1, i - 1, j)
         f(7 ,i,j) = f(5, i - 1, j)
         f(6 ,i,j) = f(8, i - 1, j)
      enddo
c     top
      f(2, i, jtop)= f(4, i, jtop + 1)
      f(6, i, jtop)= f(8, i - 1, jtop + 1)
c     bot
      f(4, i, jbot)= f(2, i, jbot - 1)
      f(7, i, jbot)= f(5, i - 1, jbot - 1)
      end

      subroutine Movie(it)
      implicit double precision(a-h, o-z)
      include 'bgk2.par'
c----------------------------------------------------------
      character raw*(2 + 8 + 4)
      character xdmf*(2 + 8 + 6)
      parameter(nf = 4)
      character fields(nf)*4
      data fields /'u', 'v', 'rho', 'vort'/
      intrinsic trim

      write(raw, '(A, I8.8, A)') 'p.', it, '.raw'
      open(77, file = raw, status = 'REPLACE', form = 'UNFORMATTED',
     $     err = 101, access = 'DIRECT',
     $     recl = nf * 8 * nx * ny)
      write(77, rec = 1) u(1:nx, 1:ny), v(1:nx, 1:ny),
     $     rho(1:nx, 1:ny),
     $     ((u(i, j + 1) - u(i, j - 1) - v(i + 1, j) + v(i - 1, j),
     $     i = 1, nx), j = 1, ny)
      close(77)

      write(xdmf, '(A, I8.8, A)') 'p.', it, '.xdmf2'
      open (77, file=xdmf, status='REPLACE')
      write(77, '(A)') '<Xdmf>'
      write(77, '(A)') '  <Domain>'
      write(77, '(A)') '    <Grid>'
      write(77, '(A)') '      <Topology'
      write(77, '(A)') '	  TopologyType="2DCoRectMesh"'
      write(77, '(    ''          Dimensions="'', I8, I8, ''"/>'')')
     $     ny + 1, nx + 1
      write(77, '(A)') '      <Geometry'
      write(77, '(A)') '	  GeometryType="ORIGIN_DXDY">'
      write(77, '(A)') '	<DataItem'
      write(77, '(A)') '	    Dimensions="2">'
      write(77, '(A)') '	  0'
      write(77, '(A)') '	  0'
      write(77, '(A)') '	</DataItem>'
      write(77, '(A)') '	<DataItem'
      write(77, '(A)') '	    Dimensions="2">'
      write(77, '(A)') '	  1'
      write(77, '(A)') '	  1'
      write(77, '(A)') '	</DataItem>'
      write(77, '(A)') '      </Geometry>'
      do 100 i = 1, nf
         write(77, '(A)') '      <Attribute'
         write(77, '(A)') '	  Center="Cell"'
         write(77, '(    ''       Name="'', A, ''">'')') trim(fields(i))
         write(77, '(A)') '	<DataItem'
         write(77, '(A)') '            Format="Binary"'
         write(77, '(A)') '            Precision="8"'
         write(77, '(    ''            Seek="'', I8, ''"'')')
     $        (i - 1) * 8 * ny * nx
         write(77, '(    ''            Dimensions="'', I8, I8, ''">'')')
     $        ny, nx
         write(77, '(    ''          '', A)') raw
         write(77, '(A)') '	</DataItem>'
         write(77, '(A)') '      </Attribute>	'
 100  continue
      write(77, '(A)') '    </Grid>'
      write(77, '(A)') '  </Domain>'
      write(77, '(A)') '</Xdmf>'
      close(77)
 101  write (*, '(''bgk2: error: fail to write output'')')
      stop 1
      end

      subroutine Diag0D(istep)
      implicit double precision(a-h, o-z)
      include 'bgk2.par'
      densit= 0.0d0
      do k = 0, npop - 1
         do j = 1, ny
            do i = 1, nx
               densit = densit + f(k ,i,j)
            enddo
         enddo
      enddo
      densit= densit/ dfloat( nx*ny) / dfloat( npop)
      umoy = 0.0d0
      vmoy = 0.0d0

      do j = 1, ny
         do i = 1, nx
            umoy = umoy + u(i,j)
            vmoy = vmoy + v(i,j)
         enddo
      enddo

      umoy = umoy / dfloat( nx*ny)
      vmoy = vmoy / dfloat( nx*ny)

      print *, 'diagnostic 0D : istep density umoy and vmoy ',
     $     istep, densit, umoy, vmoy
      end

      subroutine Config(istep, iconf)
      implicit double precision(a-h, o-z)
      include 'bgk2.par'
      iout = 60 + iconf
      do j = 2, ny - 1
         do i = 2, nx - 1
            vor = u(i, j + 1) - u(i, j - 1) -
     $           v(i + 1, j) + v(i - 1, j)
            write(iout, *) i,j, u(i,j), v(i,j), vor
         enddo
         write( iout, '( bn) ')
      enddo
      write(6, *) 'configuration at time and file >>', istep, iout
      end
