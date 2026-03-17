      subroutine genslpsame(rlam, nnodes, tabslpsame)
!***********************************************************************
!
!     Generate table for single layer potential with
!     nnodes*nnodes Chebyshev discretization points per square
!
!     Parallelized with OpenMP: all 225 independent dcuhre calls are
!     flattened into a single parallel loop (9 neighbors x 25 target
!     points).
!
!***********************************************************************

      use omp_lib
      implicit none
      external fslpsame
      integer nnodes
      integer key, n, nf, ndim, mincls, maxcls, nw
      integer kx, ky, ii, jj, nn, npoint, ic, jc, ic2, jc2, n1, n2
      integer itype, iwork, nwork, irem
      parameter (ndim = 2, nw = 1000000)
      real *8 absreq, relreq, rlam
      real *8 pi
      real *8 x(nnodes)
      real *8 cx, cy
      real *8, allocatable :: umat(:,:)
      real *8, allocatable :: vmat(:,:)
      real *8 xnodes(100)
      real *8 weights(100)
      real *8 tabslpsame(nnodes*nnodes,nnodes*nnodes,9)
      real *8 a(ndim), b(ndim)

!     Thread-private work arrays
      real *8, allocatable :: wrkstr(:)
      real *8, allocatable :: finest(:), absest(:)
      real *8, allocatable :: chebtab(:,:), chebtab2(:,:)
      real *8 fargs(4)
      real *8 xtarg, ytarg
      integer neval, ifail

      allocate(umat(nnodes,nnodes))
      allocate(vmat(nnodes,nnodes))
      itype = 2
      call chebexps(itype,nnodes,xnodes,umat,vmat,weights)
      nf = nnodes*nnodes

      pi = 4*datan(1.0d0)
      do n=1,nnodes
         x(n) = dcos((2*(nnodes-n)+1)*pi/(2*nnodes))
      enddo

      do n = 1,ndim
         a(n) =-1.0d0
         b(n) = 1.0d0
      enddo
      mincls = 0
      maxcls = 300000
      key = 0
      absreq = 0
      relreq = 1d-14

      tabslpsame = 0.0d0

!     Total work items: 9 neighbors x nnodes^2 target points
      nwork = 9 * nnodes * nnodes

!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP   SHARED(tabslpsame, x, umat, nnodes, nf, rlam, nwork,         &
!$OMP          a, b, mincls, maxcls, key, absreq, relreq)             &
!$OMP   PRIVATE(wrkstr, finest, absest, chebtab, chebtab2,             &
!$OMP           fargs, xtarg, ytarg, neval, ifail,                     &
!$OMP           iwork, irem, nn, ii, jj, cx, cy,                       &
!$OMP           kx, ky, npoint, n, n1, n2,                             &
!$OMP           ic, jc, ic2, jc2)

      allocate(wrkstr(nw))
      allocate(finest(100), absest(100))
      allocate(chebtab(nnodes,nnodes), chebtab2(nnodes,nnodes))

!$OMP DO SCHEDULE(DYNAMIC)
      do iwork = 1, nwork
!        Decode flat index into (nn, npoint)
         irem = iwork - 1
         nn     = irem / (nnodes*nnodes) + 1
         npoint = mod(irem, nnodes*nnodes) + 1

!        nn=1..(ii=-1,jj=-1), nn=2..(ii=0,jj=-1), ..., nn=9..(ii=1,jj=1)
         jj = (nn-1)/3 - 1
         ii = mod(nn-1, 3) - 1
         cx = dble(ii)
         cy = dble(jj)

         ky = (npoint-1)/nnodes + 1
         kx = npoint - nnodes*(ky-1)
         xtarg = 2.0d0*cx + x(kx)
         ytarg = 2.0d0*cy + x(ky)
         fargs(1) = xtarg
         fargs(2) = ytarg
         fargs(3) = nnodes
         fargs(4) = rlam
         call dcuhre(ndim,nf,a,b,mincls,maxcls,fslpsame,fargs,        &
           absreq,relreq,key,nw,0,finest,absest,neval,ifail,          &
           wrkstr)
         n = 0
         do jc = 1,nnodes
         do ic = 1,nnodes
            n = n+1
            chebtab(ic,jc)=finest(n)
         enddo
         enddo
         do jc = 1,nnodes
         do ic = 1,nnodes
            chebtab2(ic,jc)=0.0d0
            do jc2 = 1,nnodes
            do ic2 = 1,nnodes
               chebtab2(ic,jc)=chebtab2(ic,jc)+                       &
                 chebtab(ic2,jc2)*umat(ic2,ic)*umat(jc2,jc)
            enddo
            enddo
         enddo
         enddo
         do n2 = 1,nnodes
         do n1 = 1,nnodes
           n = (n2-1)*nnodes+n1
           tabslpsame(npoint,n,nn) = chebtab2(n1,n2)
         enddo
         enddo
      enddo
!$OMP END DO

      deallocate(wrkstr, finest, absest, chebtab, chebtab2)

!$OMP END PARALLEL

      return
      end


      subroutine fslpsame(ndim, z, nfun, fargs, f)
      implicit real *8 (a-h,o-z)
      integer ndim, nfun
      real *8 z(ndim), f(nfun), rx, ry, pi
      real *8 rr,reps,dgreen
      real *8 tx(100), ty(100)
      real *8 xtarg, ytarg, fargs(4)
      pi = 4.0d0*datan(1.0d0)

      xtarg = fargs(1)
      ytarg = fargs(2)
      rlam = fargs(4)
      rx = z(1) - xtarg
      ry = z(2) - ytarg
      rr = dsqrt(rx*rx +ry*ry)
      reps = 1.0d-14
      norder = nint(fargs(3))

      tx(1) = 1.0d0
      tx(2) = z(1)
      ty(1) = 1.0d0
      ty(2) = z(2)
      do n = 3,norder
      tx(n) = 2.0d0*z(1)*tx(n-1) - tx(n-2)
      ty(n) = 2.0d0*z(2)*ty(n-1) - ty(n-2)
      enddo

      if (rr.gt.reps) then
      dgreen = (exp(-rlam*rr)/rr)/(4.0d0*pi)
         inext = 0
         do k = 1,norder
         do j = 1,norder
            inext = inext+1
            f(inext)= tx(j)*ty(K) * dgreen
         enddo
         enddo
      else
         inext = 0
         do k = 1,norder
         do j = 1,norder
            inext = inext+1
            f(inext)= 0.0d0
         enddo
         enddo
      endif
      end
