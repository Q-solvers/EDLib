subroutine darnoldi(nloc2,vout,eout,ncv, Hstate0, nev,ierr,info)
  implicit none
  double precision, dimension(1), intent(inout) :: vout
  double precision, intent(inout) :: Hstate0
  double precision, dimension(1), intent(inout) :: eout !dimension(:),
  integer, intent(in) :: nloc2, nev
  integer, intent(inout) :: ncv
  integer, intent(out) ::ierr,info
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!

  integer       rc, myid

  Double precision&
          v(nloc2,ncv), workl(ncv*(ncv+8)),&
          workd(3*nloc2), d(ncv,2),&
          ax(nloc2),resid(nloc2),&
          temp, Zero,xnorm

  logical      select(ncv)
  integer      iparam(11), ipntr(11)
  integer v_win
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
  character    bmat*1, which*2
  integer      ido,ncv1, lworkl,&
          j,  ishfts, maxitr, mode1, nconv,ierr1,info1, i,iwsize
  logical      rvec
  Double precision&
          tol, sigma, tmpnorm
!, vtime1, vtime2, vtime
!
!     %-----------%
!     |  buffer   |
!     %-----------%
!
  double precision rbuf(nloc2),lbuf(nloc2)
  integer      displs(1), displ, rnloc, recvcounts(1)
!
!     %------------%
!     | Parameters |
!     %------------%

double precision pdsaupdtime, time1,time2
!
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
  Double precision dnrm2
  external     dnrm2, daxpy
!
!     %--------------------%
!     | Intrinsic function |
!     %--------------------%
!
  intrinsic    abs

  do j=1,nev
    eout(j)=0.d0
    do i=1,nloc2
      vout(i + nloc2*(j-1))=0.d0
    end do
  end do

  bmat  = 'I'
  which = 'SA'
  lworkl = ncv*(ncv+8)
  tol = 1e-14
  info = 0
  ido = 0
  v=0.d0
  resid=0.d0
  iparam=0
  ipntr=0
  ishfts = 1
  maxitr = 1000
  mode1 = 1
  !
  iparam(1) = ishfts
  !
  iparam(3) = maxitr
  !
  iparam(7) = mode1

  !
  !     %------------------------------------------------%
  !     | M A I N   L O O P (Reverse communication loop) |
  !     %------------------------------------------------%
  !
   10   continue
  !
  !        %---------------------------------------------%
  !        | Repeatedly call the routine DSAUPD and take |
  !        | actions indicated by parameter IDO until    |
  !        | either convergence is indicated or maxitr   |
  !        | has been exceeded.                          |
  !        %---------------------------------------------%
  call dsaupd ( ido, bmat, nloc2, which, nev, tol, resid,&
                        ncv, v, nloc2, iparam, ipntr, workd, workl,&
                        lworkl, info )
  if (ido .eq. -1 .or. ido .eq. 1) then
  !
  !           %--------------------------------------%
  !           | Perform matrix vector multiplication |
  !           |              y <--- OP*x             |
  !           | The user should supply his/her own   |
  !           | matrix vector multiplication routine |
  !           | here that takes workd(ipntr(1)) as   |
  !           | the input, and return the result to  |
  !           | workd(ipntr(2)).                     |
  !           %--------------------------------------%
  !
              !call av (nloc2, v_win, workd(ipntr(2)))
  !
  !           %-----------------------------------------%
  !           | L O O P   B A C K to call DSAUPD again. |
  !           %-----------------------------------------%
  !
  go to 10
  !
  end if !ido.eq.1
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
  if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message. Check the |
!        | documentation in DSAUPD. |
!        %--------------------------%
!
    if(myid.eq.0) then
      print *, ' '
      print *, ' Error with _saupd, info = ', info
      print *, ' Check documentation in _saupd ', iparam(5)
      print *, ' '
    end if
!
  else
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may be also computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        |                                           |
!        | The routine DSEUPD now called to do this  |
!        | post processing (Other modes may require  |
!        | more complicated post processing than     |
!        | mode1.)                                   |
!        |                                           |
!        %-------------------------------------------%
!
    rvec = .true.
!
    call dseupd (rvec, 'All', select, d, v, nloc2, sigma,&
              bmat, nloc2, which, nev, tol, resid, ncv, v, nloc2,&
              iparam, ipntr, workd, workl, lworkl, ierr )
!
!         %----------------------------------------------%
!         | Eigenvalues are returned in the first column |
!         | of the two dimensional array D and the       |
!         | corresponding eigenvectors are returned in   |
!         | the first NCONV (=IPARAM(5)) columns of the  |
!         | two dimensional array V if requested.        |
!         | Otherwise, an orthogonal basis for the       |
!         | invariant subspace corresponding to the      |
!         | eigenvalues in D is returned in V.           |
!         %----------------------------------------------%
!
    if ( ierr .ne. 0) then
!
!            %------------------------------------%
!            | Error condition:                   |
!            | Check the documentation of DSEUPD. |
!            %------------------------------------%
!
      print *, ' '
      print *, ' Error with _seupd, info = ', ierr
      print *, ' Check the documentation of _seupd. '
      print *, ' '
    else
      nconv =  iparam(5)
      do 20 j=1, nconv
        do i=1,nloc2
    !                   lbuf(i) = v(i,j)
          vout(i+nloc2*(j-1))=v(i,j)
        enddo
 20          continue
      call dmout(6, nconv, 2, d, ncv, -6,&
                 'Ritz values and relative residuals')
    end if !ierr
  end if !info
endsubroutine