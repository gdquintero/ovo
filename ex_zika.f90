Program ex_zika
    use sort

    implicit none 
    
    integer :: allocerr,iter,iter_sub,max_iter,max_iter_sub,i,j,kflag,outliers
    real(kind=8) :: alpha,epsilon,delta,sigmin,fxk,fxtrial,gaux1,gaux2,ebt,ti,a,b,c
    real(kind=8), allocatable :: xtrial(:),faux(:),indices(:),nu_l(:),nu_u(:),opt_cond(:)
    integer, allocatable :: Idelta(:)
    logical :: box,outlier

    ! COMMON INTEGERS
    integer :: samples,q

    ! COMMON SCALARS
    real(kind=8) :: sigma

    ! COMMON ARRAYS
    real(kind=8),   pointer :: t(:),y(:),grad(:,:),xk(:)

    ! COMMON BLOCKS
    common /integerData/ samples,q
    common /scalarData/ sigma
    common /vectorData/ t,y,xk,grad

    ! LOCAL SCALARS
    logical :: checkder
    integer :: hnnzmax,inform,jcnnzmax,m,n,nvparam
    real(kind=8) :: cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt,f,nlpsupn,snorm

    ! LOCAL ARRAYS
    character(len=80) :: specfnm,outputfnm,vparam(10)
    logical :: coded(11)
    logical,        pointer :: equatn(:),linear(:)
    real(kind=8),   pointer :: l(:),lambda(:),u(:),x(:)

    ! Set parameters
    n = 4
    samples = 34
    outlier = .true.
    outliers = 4
    q = samples - outliers
    max_iter = 1000000
    max_iter_sub = 1000
    alpha = 0.5d0
    epsilon = 1.0d-6
    delta = 1.0d-3
    sigmin = 1.0d-1

    allocate(t(samples),y(samples),x(n),xk(n-1),xtrial(n-1),l(n),u(n),&
    faux(samples),indices(samples),Idelta(samples),nu_l(n-1),nu_u(n-1),opt_cond(n-1),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error in main program'
        stop
    end if
  
    call read_data(outlier)

    ! Coded subroutines

    coded(1:6)  = .true.  ! evalf, evalg, evalh, evalc, evaljac, evalhc
    coded(7:11) = .false. ! evalfc,evalgjac,evalgjacp,evalhl,evalhlp


    ! Upper bounds on the number of sparse-matrices non-null elements
    jcnnzmax = 10000
    hnnzmax  = 10000

    ! Checking derivatives?
    checkder = .false.

    ! Parameters setting
    epsfeas   = 1.0d-08
    epsopt    = 1.0d-08
  
    efstain   = sqrt( epsfeas )
    eostain   = epsopt ** 1.5d0
  
    efacc     = sqrt( epsfeas )
    eoacc     = sqrt( epsopt )

    outputfnm = ''
    specfnm   = ''

    nvparam   = 1
    vparam(1) = 'ITERATIONS-OUTPUT-DETAIL 0' 

    !==============================================================================
    ! MAIN ALGORITHM
    !==============================================================================
    iter = 0

    ! Initial solution
    ! xk(:) = 0.05
    xk(:) = (/0.051,0.33,0.003/)

    ! Box-constrained? 
    box = .true. 

    if (box .eqv. .false.) then
        l(1:n) = -1.0d+20
        u(1:n) = 1.0d+20
    else
        l(1:n-1) = 0.0d0
        l(n) = -1.0d+20 
        u(1:n-1) = 1.0d+20
        u(n) = 0.0d0
    endif

    indices(:) = (/(i, i = 1, samples)/)

    kflag = 2

    ! Scenarios
    do i = 1, samples
        faux(i) = fi(xk,i,n)
    end do

    ! Sorting
    call DSORT(faux,indices,samples,kflag)

    ! q-Order-Value function 
    fxk = faux(q)

    call mount_Idelta(faux,indices,delta,Idelta,m)

    do
        iter = iter + 1

        allocate(equatn(m),linear(m),lambda(m),grad(m,n-1),stat=allocerr)

        if ( allocerr .ne. 0 ) then
            write(*,*) 'Allocation error in main program'
            stop
        end if

        equatn(:) = .false.
        linear(:) = .false.
        lambda(:) = 0.0d0

        a = xk(1)
        b = xk(2)
        c = xk(3)

        do i = 1, m
            ti = t(Idelta(i))
            ebt = exp(-1.0d0 * b * ti)
            gaux1 = model(xk,Idelta(i),n) - y(Idelta(i))
            gaux2 = exp((a / b) * ti * ebt + (1.0d0 / b) * ((a / b) - c) * (ebt - 1.0d0) - c * ti)

            grad(i,1) = (-1.0d0 / b) * (ebt * (ti - (1.0d0 / b)) - (1.0d0 / b))

            grad(i,2) = (a / b) * (ti**2) * ebt + (a / b**2) * ti * ebt + (a / b**3) * (ebt - 1.0d0) + &
                        (ti / b) * ((a / b) - c) * ebt + (1.0d0 / b**2) * ((a / b) - c) * (ebt - 1.0d0)

            grad(i,3) = ti + (1.0d0 / b) * (ebt - 1.0d0)

            grad(i,:) = gaux1 * gaux2 * grad(i,:)
        end do

        sigma = sigmin

        iter_sub = 1
        x(:) = (/xk(:),0.0d0/)
        ! x = 0.0d0

        ! Minimizing using ALGENCAN
        do 
            call algencan(myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc,   &
                myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp,jcnnzmax,    &
                hnnzmax,epsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm,   &
                specfnm,nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded,    &
                checkder,f,cnorm,snorm,nlpsupn,inform)

            xtrial(1:n-1) = x(1:n-1)

            indices(:) = (/(i, i = 1, samples)/)

            ! Scenarios
            do i = 1, samples
                faux(i) = fi(xtrial,i,n)
            end do

            ! Sorting
            call DSORT(faux,indices,samples,kflag)
    
            fxtrial = faux(q)
    
            ! Test the sufficient descent condition
            if (fxtrial .le. (fxk - alpha * norm2(xtrial(1:n-1) - xk(1:n-1))**2)) exit
            if (iter_sub .ge. max_iter_sub) exit

            sigma = 2.0d0 * sigma
            iter_sub = iter_sub + 1
        end do ! End of internal iterations

        opt_cond(:) = 0.0d0
        nu_l(:) = 0.0d0
        nu_u(:) = 0.0d0

        do j = 1, n-1
            if (xtrial(j) .le. l(j)) then 
                do i = 1, m
                    nu_l(j) = nu_l(j) + lambda(i) * grad(i,j)
                end do
            else if (xtrial(j) .ge. u(j)) then
                do i = 1, m
                    nu_u(j) = nu_u(j) - lambda(i) * grad(i,j)
                end do
            end if
        end do

        do i = 1, m
            opt_cond(:) = opt_cond(:) + lambda(i) * grad(i,:)
        enddo

        opt_cond(:) = opt_cond(:) + nu_u(:) - nu_l(:)

        print*, iter, iter_sub, fxtrial, norm2(opt_cond), m

        if (norm2(opt_cond) .le. epsilon) exit
        if (iter .ge. max_iter) exit

        deallocate(lambda,equatn,linear,grad)
        
        xk(1:n-1) = xtrial(1:n-1)
        fxk = fxtrial

        call mount_Idelta(faux,indices,delta,Idelta,m)

    end do ! End of Main Algorithm

    print*, xk
    call export(xk,n)

    CONTAINS

    !==============================================================================
    ! EXPORT RESULT TO PLOT
    !==============================================================================
    subroutine export(xsol,n)
        implicit none

        integer,        intent(in) :: n
        real(kind=8),   intent(in) :: xsol(n-1)

        Open(Unit = 10, File = "output/xstarovo.txt", ACCESS = "SEQUENTIAL")

        write(10,*) xsol(1)
        write(10,*) xsol(2)
        write(10,*) xsol(3)

        close(10)

    end subroutine export

    !==============================================================================
    ! MOUNT THE SET OF INDICES I(x,delta)
    !==============================================================================
    subroutine mount_Idelta(f,indices,delta,Idelta,m)
        implicit none

        real(kind=8),   intent(in) :: delta,f(samples),indices(samples)
        integer,        intent(out) :: Idelta(samples),m
        integer :: i
        real(kind=8) :: fq

        Idelta(:) = 0
        fq = f(q)
        m = 0

        do i = q, samples
            if (abs(fq - f(i)) .le. delta) then
                m = m + 1
                Idelta(m) = int(indices(i))
            else
                exit
            end if
        end do

        do i = q-1, 1, -1
            if (abs(fq - f(i)) .le. delta) then
                m = m + 1
                Idelta(m) = int(indices(i))
            else
                exit
            end if
        end do

    end subroutine

    !==============================================================================
    ! QUADRATIC ERROR OF EACH SCENARIO
    !==============================================================================
    function fi(x,i,n) result (res)
        implicit none

        integer,        intent(in) :: n,i
        real(kind=8),   intent(in) :: x(n-1)
        real(kind=8) :: res
        
        res = model(x,i,n) - y(i)
        res = 0.5d0 * (res**2)

    end function fi

    !==============================================================================
    ! MODEL TO BE FITTED TO THE DATA
    !==============================================================================
    function model(x,i,n) result(res)
        implicit none 

        integer,        intent(in) :: n,i
        real(kind=8),   intent(in) :: x(n-1)
        real(kind=8) :: res,a,b,c,ti,ebt

        a = x(1)
        b = x(2)
        c = x(3)
        ti = t(i)
        ebt = exp(-1.0d0 * b * ti)

        res = (a / b) * ti * ebt + (1.0d0 / b) * ((a / b) - c) * (ebt - 1.0d0) - c * ti
        res = 1.0d0 - exp(res)

    end function model

    !==============================================================================
    ! READ THE DATA
    !==============================================================================
    subroutine read_data(outlier)
        implicit none

        logical, intent(in) :: outlier
        ! SCALARS
        integer :: i

        if (outlier) then
            Open(Unit = 10, File = "output/zika_outliers.txt", ACCESS = "SEQUENTIAL")
        else
            Open(Unit = 10, File = "output/zika.txt", ACCESS = "SEQUENTIAL")
        end if

        do i = 1, samples
            read(10,*) t(i), y(i)
        enddo

        close(10)
    end subroutine read_data

    !==============================================================================
    ! SUBROUTINES FOR ALGENCAN
    !==============================================================================

    !******************************************************************************
    ! OBJECTIVE FUNCTION
    !******************************************************************************
    subroutine myevalf(n,x,f,flag)
        implicit none

        ! SCALAR ARGUMENTS
        integer, intent(in) :: n
        integer, intent(out) :: flag
        real(kind=8), intent(out) :: f

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x(n)

        ! Compute objective function

        flag = 0

        f =  x(n)

    end subroutine myevalf

    !******************************************************************************
    ! GRADIENT OF THE OBJECTIVE FUNCTION
    !******************************************************************************
    subroutine myevalg(n,x,g,flag)
        implicit none

        ! SCALAR ARGUMENTS
        integer, intent(in) :: n
        integer, intent(out) :: flag

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: g(n)

        ! Compute gradient of the objective function

        flag = 0

        g(1:n-1) = 0.0d0
        g(n)     = 1.0d0

    end subroutine myevalg

    !******************************************************************************
    ! HESSIAN FOR THE OBJECTIVE FUNCTION
    !******************************************************************************
    subroutine myevalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)
        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(out) :: lmem
        integer, intent(in) :: lim,n
        integer, intent(out) :: flag,hnnz

        ! ARRAY ARGUMENTS
        integer, intent(out) :: hcol(lim),hrow(lim)
        real(kind=8), intent(in)  :: x(n)
        real(kind=8), intent(out) :: hval(lim)

        ! Compute (lower triangle of the) Hessian of the objective function
        flag = 0
        lmem = .false.
        hnnz = 0
    end subroutine myevalh

    !******************************************************************************
    ! CONSTRAINTS
    !******************************************************************************
    subroutine myevalc(n,x,ind,c,flag)
        implicit none

        ! SCALAR ARGUMENTS
        integer, intent(in) :: ind,n
        integer, intent(out) :: flag
        real(kind=8), intent(out) :: c

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x(n)

        ! Compute ind-th constraint
        flag = 0

        c = dot_product(x(1:n-1) - xk(1:n-1),grad(ind,1:n-1)) + (sigma * 0.5d0) * &
            (norm2(x(1:n-1) - xk(1:n-1))**2) - x(n)

    end subroutine myevalc

    !******************************************************************************
    ! JACOBIAN OF THE CONSTRAINTS
    !******************************************************************************
    subroutine myevaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(out) :: lmem
        integer, intent(in) :: ind,lim,n
        integer, intent(out) :: flag,jcnnz

        ! ARRAY ARGUMENTS
        integer, intent(out) :: jcvar(lim)
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: jcval(lim)

        integer :: i

        flag = 0
        lmem = .false.

        jcnnz = n

        if ( jcnnz .gt. lim ) then
            lmem = .true.
            return
        end if

        jcvar(1:n) = (/(i, i = 1, n)/)
        jcval(1:n) = (/(grad(ind,i) + sigma * (x(i) - xk(i)), i = 1, n-1), -1.0d0/)

    end subroutine myevaljac

    !******************************************************************************
    ! HESSIAN OF THE CONSTRAINTS
    !******************************************************************************
    subroutine myevalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(out) :: lmem
        integer, intent(in) :: ind,lim,n
        integer, intent(out) :: flag,hcnnz

        ! ARRAY ARGUMENTS
        integer, intent(out) :: hccol(lim),hcrow(lim)
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: hcval(lim)

        flag = 0
        lmem = .false.
    
        hcnnz = n - 1
    
        if ( hcnnz .gt. lim ) then
            lmem = .true.
            return
        end if
    
        hcrow(1:n-1) = (/(i, i = 1, n-1)/)
        hccol(1:n-1) = (/(i, i = 1, n-1)/)
        hcval(1:n-1) = sigma

    end subroutine myevalhc

    ! ******************************************************************
    ! ******************************************************************

    subroutine myevalfc(n,x,f,m,c,flag)

        implicit none

        ! SCALAR ARGUMENTS
        integer, intent(in) :: m,n
        integer, intent(out) :: flag
        real(kind=8), intent(out) :: f

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: c(m)

        flag = - 1

    end subroutine myevalfc

    ! ******************************************************************
    ! ******************************************************************

    subroutine myevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag)

        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(out) :: lmem
        integer, intent(in) :: lim,m,n
        integer, intent(out) :: flag,jcnnz

        ! ARRAY ARGUMENTS
        integer, intent(out) :: jcfun(lim),jcvar(lim)
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: g(n),jcval(lim)

        flag = - 1

    end subroutine myevalgjac

    ! ******************************************************************
    ! ******************************************************************

    subroutine myevalgjacp(n,x,g,m,p,q,work,gotj,flag)

        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(inout) :: gotj
        integer, intent(in) :: m,n
        integer, intent(out) :: flag
        character, intent(in) :: work

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(inout) :: p(m),q(n)
        real(kind=8), intent(out) :: g(n)

        flag = - 1

    end subroutine myevalgjacp

    ! ******************************************************************
    ! ******************************************************************

    subroutine myevalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,lim,lmem,flag)

        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(out) :: lmem
        integer, intent(in) :: lim,m,n
        integer, intent(out) :: flag,hlnnz
        real(kind=8), intent(in) :: sf

        ! ARRAY ARGUMENTS
        integer, intent(out) :: hlcol(lim),hlrow(lim)
        real(kind=8), intent(in) :: lambda(m),sc(m),x(n)
        real(kind=8), intent(out) :: hlval(lim)

        flag = - 1

    end subroutine myevalhl

    ! ******************************************************************
    ! ******************************************************************

    subroutine myevalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(inout) :: goth
        integer, intent(in) :: m,n
        integer, intent(out) :: flag
        real(kind=8), intent(in) :: sf

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: lambda(m),p(n),sc(m),x(n)
        real(kind=8), intent(out) :: hp(n)

        flag = - 1

    end subroutine myevalhlp
end Program ex_zika
