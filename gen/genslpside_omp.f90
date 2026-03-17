      subroutine genslpside(rlam, nnodes, tabslpside)
!***********************************************************************
!
! This subroutine generates quadratures for SLP due to polynomial basis
! functions in two dimensions over the unit square centered at the
! origin at K by K grids on touching squares on different sides of cube.
!
! Parallelized with OpenMP: all 1800 independent dcuhre calls are
! flattened into a single parallel loop (24 face pairs x 3 steps x
! 25 target points).
!
!***********************************************************************

      use omp_lib
      implicit none
      external fslpside
      integer nnodes
      integer key, n, nf, ndim, mincls, maxcls, nw
      integer ifacet, ifaces, istep
      integer nfast, nslow, npoint, n1, n2
      integer jc, ic, jc2, ic2, itype
      integer iwork, nwork, npairs, ipair, irem
      parameter (ndim = 2, nw = 1000000)
      real *8 absreq, relreq, rlam
      real *8 pi
      real *8 x(nnodes)
      real *8, allocatable :: umat(:,:)
      real *8, allocatable :: vmat(:,:)
      real *8 xnodes(100)
      real *8 weights(100)
      real *8 tabslpside(nnodes*nnodes,nnodes*nnodes,6,6,-1:1)
      real *8 a(ndim), b(ndim)

!     Valid face pair table: pairs(1,k)=ifaces, pairs(2,k)=ifacet
      integer pairs(2,24)

!     Thread-private work arrays
      real *8, allocatable :: wrkstr(:)
      real *8, allocatable :: finest(:), absest(:)
      real *8, allocatable :: chebtab(:,:), chebtab2(:,:)
      real *8 fargs(6)
      real *8 xtarg, ytarg, ztarg
      integer neval, ifail

      allocate(umat(nnodes,nnodes))
      allocate(vmat(nnodes,nnodes))
      nf = nnodes*nnodes

      itype = 2
      call chebexps(itype,nnodes,xnodes,umat,vmat,weights)
      pi = 4*datan(1.0D0)

!     Create grid pts in 1d
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
      relreq = 1.0d-14

      tabslpside = 0.0d0

!     Pre-enumerate the 24 valid (ifaces, ifacet) pairs.
!     Skip pairs where ifacet and ifaces are on opposite faces:
!       faces 1,2 are opposite; faces 3,4 are opposite; faces 5,6 are opposite.
      npairs = 0
      do ifaces = 1,6
         do ifacet = 1,6
            if ( (ifaces.le.2 .and. ifacet.le.2) .or.                  &
                 (ifaces.ge.3 .and. ifaces.le.4 .and.                  &
                  ifacet.ge.3 .and. ifacet.le.4) .or.                  &
                 (ifaces.ge.5 .and. ifacet.ge.5) ) cycle
            npairs = npairs + 1
            pairs(1,npairs) = ifaces
            pairs(2,npairs) = ifacet
         enddo
      enddo

!     Total work items: npairs * 3 steps * nnodes^2 target points
      nwork = npairs * 3 * nnodes * nnodes

!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP   SHARED(tabslpside, x, umat, nnodes, nf, rlam, pairs, npairs, &
!$OMP          nwork, a, b, mincls, maxcls, key, absreq, relreq)      &
!$OMP   PRIVATE(wrkstr, finest, absest, chebtab, chebtab2,             &
!$OMP           fargs, xtarg, ytarg, ztarg, neval, ifail,              &
!$OMP           iwork, ipair, irem, ifaces, ifacet, istep,             &
!$OMP           npoint, nslow, nfast, n, n1, n2,                       &
!$OMP           ic, jc, ic2, jc2)

      allocate(wrkstr(nw))
      allocate(finest(100), absest(100))
      allocate(chebtab(nnodes,nnodes), chebtab2(nnodes,nnodes))

!$OMP DO SCHEDULE(DYNAMIC)
      do iwork = 1, nwork
!        Decode flat index into (ipair, istep, npoint)
         irem = iwork - 1
         ipair  = irem / (3*nnodes*nnodes) + 1
         irem   = mod(irem, 3*nnodes*nnodes)
         istep  = irem / (nnodes*nnodes) - 1
         npoint = mod(irem, nnodes*nnodes) + 1

         ifaces = pairs(1,ipair)
         ifacet = pairs(2,ipair)
         nslow = (npoint-1)/nnodes + 1
         nfast = npoint - nnodes*(nslow-1)

         if (ifacet.eq.1) then
            xtarg = -1.0d0
            ytarg = x(nslow)
            ztarg = x(nfast)
            if (ifaces.eq.3) ztarg = x(nfast)+istep*2.0d0
            if (ifaces.eq.4) ztarg = x(nfast)+istep*2.0d0
            if (ifaces.eq.5) ytarg = x(nslow)+istep*2.0d0
            if (ifaces.eq.6) ytarg = x(nslow)+istep*2.0d0
         else if (ifacet.eq.2) then
            xtarg = 1.0d0
            ytarg = x(nfast)
            ztarg = x(nslow)
            if (ifaces.eq.3) ztarg = x(nslow)+istep*2.0d0
            if (ifaces.eq.4) ztarg = x(nslow)+istep*2.0d0
            if (ifaces.eq.5) ytarg = x(nfast)+istep*2.0d0
            if (ifaces.eq.6) ytarg = x(nfast)+istep*2.0d0
         else if (ifacet.eq.3) then
            ytarg = -1.0d0
            ztarg = x(nslow)
            xtarg = x(nfast)
            if (ifaces.eq.1) ztarg = x(nslow)+istep*2.0d0
            if (ifaces.eq.2) ztarg = x(nslow)+istep*2.0d0
            if (ifaces.eq.5) xtarg = x(nfast)+istep*2.0d0
            if (ifaces.eq.6) xtarg = x(nfast)+istep*2.0d0
         else if (ifacet.eq.4) then
            ytarg = 1.0d0
            ztarg = x(nfast)
            xtarg = x(nslow)
            if (ifaces.eq.1) ztarg = x(nfast)+istep*2.0d0
            if (ifaces.eq.2) ztarg = x(nfast)+istep*2.0d0
            if (ifaces.eq.5) xtarg = x(nslow)+istep*2.0d0
            if (ifaces.eq.6) xtarg = x(nslow)+istep*2.0d0
         else if (ifacet.eq.5) then
            ztarg = -1.0d0
            ytarg = x(nfast)
            xtarg = x(nslow)
            if (ifaces.eq.1) ytarg = x(nfast)+istep*2.0d0
            if (ifaces.eq.2) ytarg = x(nfast)+istep*2.0d0
            if (ifaces.eq.3) xtarg = x(nslow)+istep*2.0d0
            if (ifaces.eq.4) xtarg = x(nslow)+istep*2.0d0
         else if (ifacet.eq.6) then
            ztarg = 1.0d0
            ytarg = x(nslow)
            xtarg = x(nfast)
            if (ifaces.eq.1) ytarg = x(nslow)+istep*2.0d0
            if (ifaces.eq.2) ytarg = x(nslow)+istep*2.0d0
            if (ifaces.eq.3) xtarg = x(nfast)+istep*2.0d0
            if (ifaces.eq.4) xtarg = x(nfast)+istep*2.0d0
         endif
         fargs(1) = xtarg
         fargs(2) = ytarg
         fargs(3) = ztarg
         fargs(4) = nnodes
         fargs(5) = ifaces
         fargs(6) = rlam
         call dcuhre(ndim,nf,a,b,mincls,maxcls,fslpside,              &
           fargs,absreq,relreq,key,nw,0,finest,absest,neval,          &
           ifail, wrkstr)

!        Having integrated polynomials, store in array chebtab
         n = 0
         do jc = 1,nnodes
         do ic = 1,nnodes
            n = n+1
            chebtab(ic,jc)=finest(n)
         enddo
         enddo

!        Convert to quad weights for tensor product grid instead of
!        polynomial coefficients.
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
            tabslpside(npoint,n,ifacet,ifaces,istep) =                &
              chebtab2(n1,n2)
         enddo
         enddo
      enddo
!$OMP END DO

      deallocate(wrkstr, finest, absest, chebtab, chebtab2)

!$OMP END PARALLEL

      return
      end


      subroutine fslpside(ndim, z, nfun, fargs, f)
      implicit none
      integer ndim, nfun,nnodes,inext,j,k,ifaces
      real *8 z(ndim), f(nfun), rx, ry, rz, pi
      real *8 rr,reps,dgreen,rlam
      real *8 tx(100), ty(100)
      real *8 xtarg,ytarg,ztarg,fargs(6)
      pi = 4.0d0*datan(1.0d0)

      xtarg = fargs(1)
      ytarg = fargs(2)
      ztarg = fargs(3)
      rlam = fargs(6)
      nnodes = nint(fargs(4))
      ifaces = nint(fargs(5))

      if (ifaces.eq.1) then
         rx = xtarg + 1.0d0
         ry = ytarg - z(2)
         rz = ztarg - z(1)
      else if (ifaces.eq.2) then
         rx = xtarg - 1.0d0
         ry = ytarg - z(1)
         rz = ztarg - z(2)
      else if (ifaces.eq.3) then
         rx = xtarg - z(1)
         ry = ytarg + 1.0d0
         rz = ztarg - z(2)
      else if (ifaces.eq.4) then
         rx = xtarg -  z(2)
         ry = ytarg - 1.0d0
         rz = ztarg -z(1)
      else if (ifaces.eq.5) then
         rx = xtarg -z(2)
         ry = ytarg - z(1)
         rz = ztarg + 1.0d0
      else if (ifaces.eq.6) then
         rx = xtarg -z(1)
         ry = ytarg - z(2)
         rz = ztarg - 1.0d0
      endif
      rr = dsqrt(rx*rx+ry*ry+rz*rz)
      reps = 1.0d-16

      tx(1) = 1.0d0
      tx(2) = z(1)
      do j = 3,nnodes
         tx(j) = 2*z(1)*tx(j-1)-tx(j-2)
      enddo

      ty(1) = 1.0d0
      ty(2) = z(2)
      do j = 3,nnodes
         ty(j) = 2*z(2)*ty(j-1)-ty(j-2)
      enddo

      if (rr.gt.reps) then
         dgreen = (exp(-rlam*rr)/rr)/(4.0d0*pi)
         inext = 0
         do k = 1,nnodes
         do j = 1,nnodes
            inext = inext+1
            f(inext)= tx(j)*ty(k) * dgreen
         enddo
         enddo
      else
         inext = 0
         do k = 1,nnodes
         do j = 1,nnodes
            inext = inext+1
            f(inext)= 0.0d0
         enddo
         enddo
      endif
      return
      end
