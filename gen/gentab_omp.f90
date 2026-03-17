      subroutine gentab_omp(rlam, nnodes, nthreads,                     &
     &   tabslpsame, tabslpside, tabdlpside)
!***********************************************************************
!
!     Wrapper that sets the OpenMP thread count and calls the three
!     table generators: genslpsame, genslpside, gendlpside.
!
!     Input:
!       rlam     - Yukawa parameter (scaled by caller)
!       nnodes   - number of Chebyshev nodes per dimension
!       nthreads - number of OpenMP threads (0 = use current default)
!
!     Output:
!       tabslpsame(nf,nf,9)       - SLP same-face quadrature weights
!       tabslpside(nf,nf,6,6,3)   - SLP side quadrature weights
!       tabdlpside(nf,nf,6,6,3)   - DLP side quadrature weights
!       where nf = nnodes*nnodes
!
!***********************************************************************

      use omp_lib
      implicit none
      integer nnodes, nthreads
      real *8 rlam
      real *8 tabslpsame(*)
      real *8 tabslpside(*)
      real *8 tabdlpside(*)

      if (nthreads .gt. 0) call omp_set_num_threads(nthreads)

      call genslpsame(rlam, nnodes, tabslpsame)
      call genslpside(rlam, nnodes, tabslpside)
      call gendlpside(rlam, nnodes, tabdlpside)

      return
      end
