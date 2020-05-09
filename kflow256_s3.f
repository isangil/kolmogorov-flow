c -----------------------------------------------------------------
c       This is the spectral code for kolmogorov flow
c       using the message-passing in Origin 2000 machines
c       isg, oct 26 98 REV aug 99
c       nx = ny = nz = 2**n   ( n = 4,5,...)
c       nx = x-diretion mesh size in each node
c       ny = y-diretion mesh size in each node
c       nz = z-direction mesh size in each node
c       nproc = number of process (nodes)
c       lx=nx,ly=ny,lz=nz/nproc
c ----------------------------------------------------------------- 

       program kolflow

       include 'include.h'

       complex,allocatable,dimension(:,:,:)::vx,vy,vz
       complex,allocatable,dimension(:,:,:)::wx,wy,wz
       complex,allocatable,dimension(:,:,:)::ox,oy,oz
       real,allocatable,dimension(:,:,:)::kx,ky,kz,k2
       real,allocatable,dimension(:,:,:)::tmp,tmp1,k2_e
       integer,allocatable,dimension(:)::iseed 
       complex, allocatable, dimension(:)::coeffy 
       real, allocatable, dimension(:)::coeffxz
       complex ff
       real ffr

c  setup MPP environment

       call mpi_init(ierror)
       call mpi_comm_rank(mpi_comm_world,id,ierror)
       call mpi_comm_size(mpi_comm_world,nproc,ierror)
       nallgrp = mpi_comm_world

c   read parameters from node id = 0

       ak0 = 4.75683

       if(id.eq.0)then
         open(90,file='parameter.d',status='unknown')
         write(*,*)'enter nx,ny,nz'
         read(90,*)nx,ny,nz
         write(*,*)'enter number of processors'
         read(90,*)nprocs
         write(*,*) 'enter seed (321)'
         read(90,*) iseed_m0
         write(*,*)'enter nstep'
         read(90,*) nstep
         write(*,*)'enter tout'
         read(90,*) itout
         write(*,*)'enter eout'
         read(90,*) ieout
         write(*,*)'enter dt'
         read(90,*) dt
         write(*,*)'enter rnu'
         read(90,*) rnu
         write(*,*)'enter u0'
         read(90,*) u0
         write(*,*)'enter time'
         read(90,*) time
         write(*,*)'enter new'
         read(90,*) new
         write(*,*)'enter idp for rerun'
         read(90,*) idp
         write(*,*)'enter force'
         read(90,*) ffr
         write(*,*)'enter kolmogorov forcing mode'
         read(90,*) is
         write(*,*)nx,ny,nz,nprocs,iseed_m0,nstep,itout,ieout,dt,rnu,
     1     u0,time,new,idp,ffr,is
       endif 

       call mpi_bcast(nx,1,MPI_INTEGER,0,nallgrp,ierror)
       call mpi_bcast(ny,1,MPI_INTEGER,0,nallgrp,ierror)
       call mpi_bcast(nz,1,MPI_INTEGER,0,nallgrp,ierror)
       call mpi_bcast(nprocs,1,MPI_INTEGER,0,nallgrp,ierror)      
       call mpi_bcast(iseed_m0,1,MPI_INTEGER,0,nallgrp,ierror)
       call mpi_bcast(nstep,1,MPI_INTEGER,0,nallgrp,ierror)
       call mpi_bcast(itout,1,MPI_INTEGER,0,nallgrp,ierror)
       call mpi_bcast(ieout,1,MPI_INTEGER,0,nallgrp,ierror)
       call mpi_bcast(dt,1,MPI_REAL,0,nallgrp,ierror)
       call mpi_bcast(rnu,1,MPI_REAL,0,nallgrp,ierror)
       call mpi_bcast(u0,1,MPI_REAL,0,nallgrp,ierror)
       call mpi_bcast(time,1,MPI_REAL,0,nallgrp,ierror)
       call mpi_bcast(new,1,MPI_LOGICAL,0,nallgrp,ierror)
       call mpi_bcast(idp,1,MPI_INTEGER,0,nallgrp,ierror)
       call mpi_bcast(ffr,1,MPI_REAL,0,nallgrp,ierror)
       call mpi_bcast(is,1,MPI_INTEGER,0,nallgrp,ierror)

       lx = nx
       ly = ny
       lz = nz/nproc
       lx1=lx+1
       lx2= lx*2
       lly=ly/nproc
       pi=4.*atan(1.)
       scale = 1.0/(nx*ny*nz)
     
c allocate memory 
      
       allocate( vx(lx1,ly,lz) )
       allocate( vy(lx1,ly,lz) )
       allocate( vz(lx1,ly,lz) )
       allocate( wx(lx1,ly,lz) )
       allocate( wy(lx1,ly,lz) )
       allocate( wz(lx1,ly,lz) )
       allocate( ox(lx1,ly,lz) )
       allocate( oy(lx1,ly,lz) )
       allocate( oz(lx1,ly,lz) )
       allocate( kx(lx1,ly,lz) )
       allocate( ky(lx1,ly,lz) )
       allocate( kz(lx1,ly,lz) )
       allocate( k2(lx1,ly,lz) )
       allocate( k2_e(lx1,ly,lz) )
       allocate( tmp(lx1,ly,lz) )
       allocate( tmp1(lx1,ly,lz) )
       allocate( iseed(nproc) )
       allocate( coeffy( ly+15 ) )
       allocate( coeffxz( (lx2+15)+2*(nz+15) ) )

       if(id.eq.0)print*,' memory allocated'

       ff = cmplx(ffr,0.0)
       iii = 1
       jjj = 1
       kkk = 1

       eout = ieout * dt
       tout = itout * dt
       nek=0.9*0.5*sqrt(4.0*nx*nx+ny*ny+nz*nz)
        
       if(id.eq.0)then
         write(*,*)ieout,itout,nek
         open(20,file='spectrum.d',status='unknown') 
         open(22,file='transfer.d',status='unknown') 
         open(70,file='initial.sp',status='unknown')
       end if
        
c---Generate initial condition -------

       if(new) then

c---Generate random number seed for each process ---

         if(id.eq.0)then
           call srand(iseed_m0)
           do i=1,nproc
             iseed(i) = irand()
           end do
         end if 
         call mpi_bcast(iseed,nproc,MPI_INTEGER,0,nallgrp,ierror)
         irr = iseed(id+1)
         call srand(irr)

         call gaussian(vx,pi)
         call gaussian(vy,pi)
         call gaussian(vz,pi)

         call wavenumber(kx,ky,kz,k2,k2_e,dt,rnu,id)

         tmp = (kx*real(vx) + ky*real(vy) + kz*real(vz))/k2
         vx = cmplx(real(vx) - kx*tmp, aimag(vx))
         vy = cmplx(real(vy) - ky*tmp, aimag(vy))
         vz = cmplx(real(vz) - kz*tmp, aimag(vz))
         tmp = (kx*aimag(vx) + ky*aimag(vy) + kz*aimag(vz))/k2
         vx = cmplx(real(vx), aimag(vx) - kx*tmp)
         vy = cmplx(real(vy), aimag(vy) - ky*tmp)
         vz = cmplx(real(vz), aimag(vz) - kz*tmp)
      
         cc1 = u0 * sqrt(8.*sqrt(2./pi)/(3.*pi*ak0**5))
         tmp = cc1 * sqrt(k2) * exp (-k2/(ak0*ak0))

         vx = vx * tmp
         vy = vy * tmp
         vz = vz * tmp

         call symmetrize(vx,id)
         call symmetrize(vy,id)
         call symmetrize(vz,id)

          call dealiasing(vx,k2)
          call dealiasing(vy,k2)
          call dealiasing(vz,k2)
        
c     calculating initial spectrum

         tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
         tmp(1,:,:) = 0.5*tmp(1,:,:)
         do i=1,nek
           ek = 0.0
           e_t = 0.0
           ek=sum(tmp(1:lx,:,:),mask=(abs(sqrt(k2)-i).lt.0.5))
           call mpi_reduce(ek,e_t,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)
           if(id.eq.0)then
             write(70,*)i,e_t
           end if
         end do
         if(id.eq.0)close(70)
          if (id.eq.0)write(*,*)'initial spec done'
       else

c ------ Rreading from DISK ------ 

        call input(vx,vy,vz,idp,id,nallgrp)
cc         read(40+id)(((vx(i,j,k),i=1,lx1),j=1,ly),k=1,lz)
cc         read(40+id)(((vy(i,j,k),i=1,lx1),j=1,ly),k=1,lz)
cc         read(40+id)(((vz(i,j,k),i=1,lx1),j=1,ly),k=1,lz)
         call symmetrize(vx,id)
         call symmetrize(vy,id)
         call symmetrize(vz,id)

c        write(*,*)id,sum(vx*conjg(vx)),
c    1     sum(vy*conjg(vy)),sum(vz*conjg(vz))
c       if(id.eq.0)then
C for testing
c       write(*,*)(vx(i,5,2),i=1,lx1)
c       write(*,*)(vy(i,5,2),i=1,lx1)
c       write(*,*)(vz(i,5,2),i=1,lx1)
c       endif

         if(id.eq.0)then
           write(*,*)'after reading'
         end if
c
      endif

      call wavenumber(kx,ky,kz,k2,k2_e,dt,rnu,id)

      call cfft1di(ly,coeffy)
      call scfft2dui(lx2,nz,coeffxz)

      ox=vx
      oy=vy
      oz=vz

      dt_h=.5*dt
      istep=-1
      
      if (id.eq.0) write(*,*)'begin main loop'
1     istep = istep + 1

c   write out the total energy in K space

       time1 = mpi_wtime()
       tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
       tmp(1,:,:)=0.5*tmp(1,:,:)
       ek = sum(tmp)
       call mpi_reduce(ek,e_t,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)
       call mpi_bcast(e_t,1,MPI_REAL,0,nallgrp,ierror)
	if(id.eq.0)then
         write(*,*)istep,time,e_t
         write(76,*)time,e_t
       end if
       if (e_t.gt.100.)   then
         write(*,*)'id',id, 'et way up',e_t
          call MPI_FINALIZE(ierror) 
       endif
       if(mod(istep,itout).eq.0.and.istep.ne.0) then
         call output(vx,vy,vz,iii,id,nallgrp)
         iii = iii + 1
       endif

       wx = (0.,1.) * (ky*vz - kz*vy)
       wy = (0.,1.) * (kz*vx - kx*vz)
       wz = (0.,1.) * (kx*vy - ky*vx)

c     energy spectrum calculation 
     
      if(mod(istep,ieout).eq.0) then
        tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
        tmp(1,:,:)=0.5*tmp(1,:,:)
        if(id.eq.0)then   
          write(20,*)(kkk-1)
        end if
        do i=1,nek
          ek = 0.0
          e_t = 0.0
          ek=sum(tmp,mask=(abs(sqrt(k2)-i).lt.0.5))
          call mpi_reduce(ek,e_t,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)
          if(id.eq.0)then
            write(20,*)i,e_t
          end if
        end do
        kkk = kkk + 1
      end if
     
      kx = real(vx)
      ky = aimag(vx)
      kz = real(vy)
      k2 = aimag(vy)
      k2_e = real(vz)
      tmp1 = aimag(vz)
c
      isign = 1
c    
      call mpifft(vx,isign,coeffy,coeffxz,id,nallgrp)
      call mpifft(vy,isign,coeffy,coeffxz,id,nallgrp)
      call mpifft(vz,isign,coeffy,coeffxz,id,nallgrp)      

      call mpifft(wx,isign,coeffy,coeffxz,id,nallgrp)
      call mpifft(wy,isign,coeffy,coeffxz,id,nallgrp)    
      call mpifft(wz,isign,coeffy,coeffxz,id,nallgrp)      

c  ---Check energy and vorticity in physical space -----
 
c      t_ex = 0.
c      t_ey = 0.
c      t_ez = 0.
      
c      t_wx = 0.
c      t_wy = 0.
c      t_wz = 0.

c      t_ex=sum(vx(1:lx,:,:)*vx(1:lx,:,:))
c      t_ey=sum(vy(1:lx,:,:)*vy(1:lx,:,:))
c      t_ez=sum(vz(1:lx,:,:)*vz(1:lx,:,:))

c      t_wx=sum(wx(1:lx,:,:)*wx(1:lx,:,:))
c      t_wy=sum(wy(1:lx,:,:)*wy(1:lx,:,:))
c      t_wz=sum(wz(1:lx,:,:)*wz(1:lx,:,:))

c      call mpi_reduce(t_ex,to_ex,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)
c      call mpi_reduce(t_ey,to_ey,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)
c      call mpi_reduce(t_ez,to_ez,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)
 
c      call mp_reduce(t_wx,to_wx,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)
c      call mp_reduce(t_wy,to_wy,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)
c      call mp_reduce(t_wz,to_wz,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)

c      if(id.eq.0)then
c        to_ex = to_ex*scale
c        to_ey = to_ey*scale
c        to_ez = to_ez*scale
c
c        to_wx = to_wx*scale
c        to_wy = to_wy*scale
c        to_wz = to_wz*scale
c
c        to_energy = to_ex + to_ey + to_ez
c        to_omega  = to_wx + to_wy + to_wz
c
c        write(11,*)time,to_ex,to_ey,to_ez,0.5*to_energy
c        write(11,*)time,to_wx,to_wy,to_wz,0.5*to_omega
c      end if
       
c      average over x,z & t
c
c      do i=1,ly
c        avx(i) = avx(i) + sum(vx(:,i,:))
c        avv(i) = avv(i) + sum(vy(:,i,:))
c        avw(i) = avw(i) + sum(vz(:,i,:))
c      enddo
       tmp = real(vy)*real(wz) - real(vz)*real(wy)
       wz = cmplx(real(vz)*real(wx) - real(vx)*real(wz), aimag(wz))
       wx = cmplx(real(vx)*real(wy) - real(vy)*real(wx), aimag(wx))
       wy = cmplx(tmp, aimag(wy))
       tmp = aimag(vy)*aimag(wz) - aimag(vz)*aimag(wy)
       wz = cmplx(real(wz), aimag(vz)*aimag(wx) - aimag(vx)*aimag(wz))
       wx = cmplx(real(wx), aimag(vx)*aimag(wy) - aimag(vy)*aimag(wx))
       wy = cmplx(real(wy), tmp)

       vx = cmplx(kx,ky)
       vy = cmplx(kz,k2)
       vz = cmplx(k2_e,tmp1)

       call wavenumber(kx,ky,kz,k2,k2_e,dt,rnu,id)

       isign = -1

       call mpifft(wx,isign,coeffy,coeffxz,id,nallgrp)
       call mpifft(wy,isign,coeffy,coeffxz,id,nallgrp)
       call mpifft(wz,isign,coeffy,coeffxz,id,nallgrp)

       call symmetrize(wx,id)
       call symmetrize(wy,id)
       call symmetrize(wz,id)

       if (id.eq.0) then
         wy(1,is+1,1) = wy(1,is+1,1) + ff     ! kolmogorov force now
         wy(1,ly-is+1,1) = wy(1,ly-is+1,1) + ff  ! fx=cos(is*y)
       endif
       call dealiasing(vx,k2)
       call dealiasing(vy,k2)
       call dealiasing(vz,k2)

       tmp = (kx*real(wy) + ky*real(wz) + kz*real(wx))/k2
       wy = cmplx(real(wy) - kx*tmp, aimag(wy))
       wz = cmplx(real(wz) - ky*tmp, aimag(wz))
       wx = cmplx(real(wx) - kz*tmp, aimag(wx))
       tmp = (kx*aimag(wy) + ky*aimag(wz) + kz*aimag(wx))/k2
       wy = cmplx(real(wy), aimag(wy) - kx*tmp)
       wz = cmplx(real(wz), aimag(wz) - ky*tmp)
       wx = cmplx(real(wx), aimag(wx) - kz*tmp)

       if(mod(istep,ieout).eq.0) then
         tmp = real(vx)*real(wy) + real(vy)*real(wz) + 
     &        real(vz)*real(wx) +
     &        aimag(vx)*aimag(wy) + aimag(vy)*aimag(wz) + 
     &        aimag(vz)*aimag(wx)
         tmp(1,:,:) = 0.5*tmp(1,:,:)
         if(id.eq.0)then
           write(22,*)(jjj-1)
         end if
 
         do i=1,nek
           ek = 0.
           e_t = 0. 
           ek=sum(tmp(1:lx,:,:),mask=(abs(sqrt(k2)-i).lt.0.5))
           call mpi_reduce(ek,e_t,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)
           if(id.eq.0)then
             write(22,*) i,ek
           end if
         end do
         jjj = jjj + 1
       endif

      if(istep.eq.0) then
        vx = vx * (1.-rnu*dt_h*k2) + dt_h*wy
        vy = vy * (1.-rnu*dt_h*k2) + dt_h*wz
        vz = vz * (1.-rnu*dt_h*k2) + dt_h*wx

        time2=mpi_wtime()
        if(id.eq.0) then
          write(*,fmt='(5X,A,F20.4,A)')
     &   'elapsed time =',1.0D3*(time2-time1),' msec'
        end if
        call output(wx,wy,wz,0,id,nallgrp)
        time = time + dt_h
      elseif(istep.eq.1) then
        vx = ox + dt*wy - rnu*dt*k2*vx
        vy = oy + dt*wz - rnu*dt*k2*vy
        vz = oz + dt*wx - rnu*dt*k2*vz
     
        time2=mpi_wtime()
        if(id.eq.0) then
          write(*,fmt='(5X,A,F20.4,A)')
     &   'elapsed time =',1.0D3*(time2-time1),' msec'
        end if
        call input(ox,oy,oz,0,id,nallgrp)
        time = time + dt_h
      else
         vx = vx + dt_h*(3.*wy - k2_e*ox)
         vy = vy + dt_h*(3.*wz - k2_e*oy)
         vz = vz + dt_h*(3.*wx - k2_e*oz)
         vx = vx*k2_e
         vy = vy*k2_e
         vz = vz*k2_e
         ox = wy
         oy = wz
         oz = wx

         time = time + dt
         time2 = MPI_wtime()

        if(id.eq.0) then
          write(*,fmt='(5X,A,F20.4,A)')
     &   'elapsed time =',1.0D3*(time2-time1),' msec'
        end if

      endif

      if(istep.lt.nstep) goto 1
      write(*,*)id,'finished '

      if(id.eq.0)then
        close(11)
        close(20)
        close(22)
      end if

      stop

100   format(' step',i4,':  t = ',e14.6,' E = ',e14.6,' W = ',e14.6)
200   format(8e13.5)
102   format(2(1x,E16.7))

      call MPI_FINALIZE(ierror)
      stop
      end

c --------------------------------------------------------------
       subroutine wavenumber(kx,ky,kz,k2,k2_e,dt,rnu,id)

       include 'include.h'

       real,dimension(lx1,ly,lz)::kx,ky,kz,k2,k2_e

       do i=1,lx1
         kx(i,:,:)= i-1 
       end do

       do j=1,ly
         ky(:,j,:) = mod(j-1+ly/2,ly)-ly/2
       end do
   
       do k=1,lz
         k1=id*lz+k
         kz(:,:,k) = mod(k1-1+nz/2,nz)-nz/2
       end do
 
       k2=kx*kx+ky*ky+kz*kz
       k2_e = exp(-k2*dt*rnu)
 
       if(id.eq.0)then
         k2(1,1,1)=1.e-5
         k2_e(1,1,1) = 1.0
       end if

       return
       end

c --------------------------------------------------------------------
      subroutine gaussian(u,pi)

      include 'include.h'

      complex,dimension(lx1,ly,lz) :: u

      u = 0.0

      do i=1,lx1
       do j=1,ly
        do k=1,lz
          t1 = rand()
          t2 = rand()
          if(t1.le.1.e-10) t1 = 1.e-10
          if(t2.le.1.e-10) t2 = 1.e-10
          t2 = 2*pi*t2
          u(i,j,k)=sqrt(-2.0*log(t1))*cmplx(cos(t2),sin(t2))
        end do
       end do
      end do
  
      return
      end

c ---------------------------------------------------------------------
      subroutine symmetrize(c,id)

      include 'include.h'
 
      complex,dimension(lx1,ly,lz)::c
      complex,dimension(ly,lz)::cc

c   zero k2=ny/2 modes
      c(:,ly/2+1,:) = (0.0,0.0)
c   zero k3=nz/2 modes
      nid=(nz/2+1)/lz
      nip=nz/2+1-nid*lz
      if(id.eq.nid) c(:,:,nip) = (0.0,0.0)
c   zero k1=nx modes
      c(lx1,:,:)=(0.0,0.0)
c  zero k1=k2=k3=0 mode
      if(id.eq.0) c(1,1,1) = (0.0,0.0)

      cc=c(1,:,:)

c for k1=0 plane, zero all modes with k2<0

c     do j = ly/2+2,ly 
c       c(1,j,:) = (0.0,0.0)
c     end do

c for the line k1=0 and k2=0,  zero all modes with k3<0
c     do k = 1,lz
c       k1=k+id*lz
c       if(k1.gt.nz/2.and.k1.le.nz) c(1,1,k)=(0.0,0.0)
c     end do

      return
      end

c ----------------------------------------------------------------------
c   2/3 rule    
c ----------------------------------------------------------------------
      subroutine dealiasing(vx,k2)
         
      include 'include.h' 

      complex,dimension(lx1,ly,lz)::vx
      real,dimension(lx1,ly,lz)::k2
      double precision ktr
  
      ktr = nz/3.0 

      where(sqrt(k2).ge.ktr)vx = (0.0,0.0)
      
      return
      end

c ----------------------------------------------------------------------
      subroutine output(ux,uy,uz,idump,id,nallgrp)

      include "include.h"

      complex,dimension(lx1,ly,lz)::ux,uy,uz
      integer i1d,i2d,i3d,i4d,i5d

      character*60 fnm1

      i1d=int(id/100)
      i2d=int((id-i1d*100)/10)
      i3d=id-i1d*100-i2d*10

      if(idump.lt.10) then
        fnm1='/scratch.o3.1/eneko/kol256/vel'//char(idump+48)//'.'
     .        //char(i1d+48)//char(i2d+48)//char(i3d+48)
      else
        i4d=int(idump/10)
        i5d=mod(idump,10)
        fnm1='/scratch.o3.1/eneko/kol256/vel'//char(i4d+48)
     .        //char(i5d+48)//'.'//char(i1d+48)
     .        //char(i2d+48)//char(i3d+48)
      endif

      open(10,file=fnm1,status='unknown',form='unformatted')
      write(10)ux
      write(10)uy
      write(10)uz
      close(10)

      return
      end

c ---------------------------------------------------------------------
      subroutine input(ux,uy,uz,idp,id,nallgrp)

      include "include.h"

      complex,dimension(lx1,ly,lz)::ux,uy,uz
      integer i1d,i2d,i3d,i4d,i5d

      character*60 fnm1

      i1d=int(id/100)
      i2d=int((id-i1d*100)/10)
      i3d=id-i1d*100-i2d*10
      if(idp.lt.10) then
        fnm1='/scratch.o3.1/eneko/kol256/vel'//char(idp+48)//'.'
     .        //char(i1d+48)//char(i2d+48)//char(i3d+48)
      else
        i4d=int(idp/10)
        i5d=mod(idp,10)
        fnm1='/scratch.o3.1/eneko/kol256/vel'//char(i4d+48)
     .        //char(i5d+48)//'.'//char(i1d+48)
     .        //char(i2d+48)//char(i3d+48)
      endif

      open(10,file=fnm1,status='unknown',form='unformatted')
      read(10)ux
      read(10)uy
      read(10)uz
      close(10)

      return
      end

c -------------------------------------------------------------------
        subroutine transpose(ux,id,nallgrp)

        include 'include.h'

        complex,dimension(lx+1,ly,lz)::ux
        complex,dimension(lx+1,ly/nproc,lz):: tmp1,tmp2
        complex,dimension(lx+1,ly/nproc,lz*nproc)::tmp
        integer isize,nzm,nzp,status(MPI_STATUS_SIZE,2),req(2)

        isize=(lx+1)*ly*lz/nproc
        lly=ly/nproc
        do i=1,nproc-1
          nzp=mod(id+i,nproc)
          nzm=mod(id-i+nproc,nproc)
          js=nzp*lly
          do j=1,lly
            j1=js+j
            tmp1(:,j,:)=ux(:,j1,:)
          end do
          call mpi_isend(tmp1,isize,MPI_COMPLEX,nzp,i,nallgrp,req(1),
     .                   ierror)
          call mpi_irecv(tmp2,isize,MPI_COMPLEX,nzm,i,
     .               nallgrp,req(2),ierror)
          call MPi_Waitall(2,req,status,ierror)
          ks=nzm*lz
          do k=1,lz
            k1=ks+k
            tmp(:,:,k1)=tmp2(:,:,k)
          end do
        end do

         ks=id*lz
         js=id*lly
        do k=1,lz
          k1=ks+k
          do j=1,ly/nproc
            j1=js+j
            tmp(:,j,k1)=ux(:,j1,k)
          end do
        end do

        do k=1,lz
          do j=1,ly
             ux(:,j,k)=tmp(:,k,j)
           end do
        end do

        return
        end

c ------------------------------------------------------------------------
       subroutine mpifft(ux,isign,coeffz,coeffxy,id,nallgrp)

       include "include.h"

       complex,dimension(lx+1,ly,lz)::ux
       complex,dimension(lx+1,ly/nproc,lz):: tmp1,tmp2
       complex,dimension(lx+1,ly/nproc,lz*nproc)::tmp
       complex, dimension(nz+15)::coeffz
       real, dimension((lx2+15)+2*(ly+15))::coeffxy
       real alpha1,alpha2

       if(isign.eq.-1) then

C      Step 1, 2D fft, real-to-complex in x, complex-to-complex in y

        alpha2=1.0/lx2/ly
c       time1=mpi_wtime()
        do k=1,lz
          call scfft2du(isign,lx2,ly,ux(:,:,k),lx2+2,coeffxy)
        end do
        ux=ux*alpha2
c       time=mpi_wtime()-time1
c       if(id.eq.0)write(*,*)'time for 2-d fft',time

c      Step 2, Transpose between z and y

c       time2=mpi_wtime()
        call transpose(ux,id,nallgrp)
c       time=mpi_wtime()-time2
c       if(id.eq.0)write(*,*)'time for transpose',time

C      Step 3, 1d complex-to-complex in z

c       time3=mpi_wtime()
        do j=1,ly/nproc 
          do i=1,lx
            call cfft1d(isign,nz,ux(i,:,j),1,coeffz)
          end do
        end do
        ux=ux/nz
c       time=mpi_wtime()-time3
c       if(id.eq.0)write(*,*)'time for 1-d fft',time

       else if(isign.eq.1)then

C      Step 1, 1d complex-to-complex in z

c       time1=mpi_wtime()
        do j=1,ly/nproc
        do i=1,lx
            call cfft1d(isign,nz,ux(i,:,j),1,coeffz)
c           ux(i,:,j)=ux(i,:,j)/nz
          end do
        end do
c       time=mpi_wtime()-time1
c       if(id.eq.0)write(*,*)'time for 1-d fft',time

c      Step 2, Reverse transpose between z and y

c       time2=mpi_wtime()
        call transpose(ux,id,nallgrp)
c       time=mpi_wtime()-time2
c       if(id.eq.0)write(*,*)'time for transpose',time

C      Step 3, 2D fft, complex-to-complex in y, complex-to-real in x

        alpha2=1.0/lx2/ly
c       time3=mpi_wtime()
        do k=1,lz
          call csfft2du(isign,lx2,ly,ux(:,:,k),lx2+2,coeffxy)
c         ux(:,:,k)=ux(:,:,k)*alpha2
c         call sscal2d(lx,ly,alpha2,ux(:,:,k),lx+2)
        end do
c       time=mpi_wtime()-time3
c       if(id.eq.0)write(*,*)'time for 2-d fft',time

       end if

       return
       end


