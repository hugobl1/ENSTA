c See qpalm_dpstrf.m for an explanation of this subroutine.

#include "fintrf.h"

C     Gateway routine

      subroutine mexfunction (nlhs,plhs,nrhs,prhs)

      implicit none

c MexFunction arguments

      mwsize nlhs, nrhs
      mwpointer plhs(*), prhs(*)

c Declarations of the MATLAB API functions used in this MEX-file

      mwpointer mxGetPr, mxCreateDoubleMatrix

c Local variables

      mwpointer A_pr,uplo_pr,tol_pr,L_pr,piv_pr,rk_pr,info_pr

      mwsize m, n, size
      mwsize mxGetM, mxGetN

      integer*4 mxIsNumeric,mxIsChar, piv(4000)
      integer info,rk
      real*8  A(4000*4000),tol,work(2*4000),piv2(4000),rk2, info2
      character uplo(1)

      mwPointer mxCreateNumericArray
      integer*4 mxClassIDFromClassName
      integer*4 classid
      integer*4 complexflag
      mwSize ndim
      mwSize dims(2)

c     TEST ECRITURE
      integer*4 mexPrintf
      integer*4 :: k,i
      real(8) :: x = 1.23d0
      character(len=80) :: line


c Check for proper number and type of arguments. 

      if (nrhs .ne. 3) then
         call mexErrMsgTxt('Three input required.')
      elseif (nlhs .gt. 4) then
         call mexErrMsgTxt('Too many output arguments.')
      endif
      if (mxIsNumeric(prhs(1)) .eq. 0) then
         call mexErrMsgTxt('Input 1 must be a number.')
      endif
      if (mxIsNumeric(prhs(3)) .eq. 0) then
         call mexErrMsgTxt('Input 3 must be a number.')
      endif
      if (mxIsChar(prhs(2)) .eq. 0) then
         call mexErrMsgTxt('Input 2 must be a number.')
      endif



c Get the size of the input array

      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      size = m*n

c Create matrix for the return argument

c      classid = mxClassIDFromClassName('int8')
c      complexflag = 0
c      ndim = 2
c      dims(1) = 1
c      dims(2) = n
c      plhs(2) = mxCreateNumericArray(ndim, dims, 
c     +          classid, 
c     +          complexflag)

c      dims(2) = 1
c      plhs(3) = mxCreateNumericArray(ndim, dims, 
c     +          classid, 
c     +          complexflag)


      plhs(1) = mxCreateDoubleMatrix(m,n,0)
      plhs(2) = mxCreateDoubleMatrix(1,n,0)
      plhs(3) = mxCreateDoubleMatrix(1,1,0)
      plhs(4) = mxCreateDoubleMatrix(1,1,0)

      A_pr = mxGetPr(prhs(1))
      uplo_pr = mxGetPr(prhs(2))
      tol_pr = mxGetPr(prhs(3))

      L_pr = mxGetPr(plhs(1))
      piv_pr = mxGetPr(plhs(2))
      rk_pr = mxGetPr(plhs(3))
      info_pr = mxGetPr(plhs(4))
            
      call mxCopyPtrToReal8(A_pr,A,size)
      call mxCopyPtrToCharacter(uplo_pr,uplo,1)
      call mxCopyPtrToReal8(tol_pr,tol,1)

c Call the computational subroutine

c      call dpstrf(y,A, m, n,uplo)
      call dpstrf2( uplo, n, A, piv, rk, tol, work, info )
c      do 30 i = 1,n
c          mexPrintf('\n TEST \m')
c      continue

c Load the data into y_pr, which is the output to MATLAB

      do i=1,n
        piv2(i) = piv(i)
      enddo
      rk2 = rk
      info2 = info


c      write(line,*) "\n\n info = ",info,"\n"
c      k = mexPrintf(line//achar(10))



      call mxCopyReal8ToPtr(A,L_pr,size)  
      call mxCopyReal8ToPtr(piv2,piv_pr,n) 
      call mxCopyReal8ToPtr(rk2,rk_pr,1) 
      call mxCopyReal8ToPtr(info2,info_pr,1) 

C    TEST ECRITURE

c      integer*4 :: k
c      character(len=80) :: line
c      real(8) :: x = 1.23d0
c      k = mexPrintf(line//achar(10))

c      mxDestroyArray(plhd(3))   
c      call mxCopyReal8ToPtr(z,z_pr,size)
      return
      end

c-----------------------------------------------------------------------
c The Fortran subroutine
c-----------------------------------------------------------------------

      subroutine dpstrf2 (uplo, n, A, piv, rk, tol, work, info)

      integer i,j
      character uplo
      real*8  A(n,n)
      
      if ((uplo .eq. 'L') .or. (uplo .eq. 'l')) then
        do i=1,n
           do j=(i+1),n
              A(i,j)= 0.0
           enddo
        enddo
      else if ((uplo .eq. 'U') .or. (uplo .eq. 'u')) then
        do j=1,n
           do i=(j+1),n
              A(i,j)= 0.0
           enddo
       enddo
      endif

      call dpstrf (uplo, n, A, n, piv, rk, tol, work, info)

      return
      end
