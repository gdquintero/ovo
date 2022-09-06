Program main
    use sort

    implicit none 
    
    integer :: allocerr,iter,int_iter,i,kflag
    integer, parameter :: maxIter = 1000,maxIntIter = 1000
    real(kind=8), parameter :: alpha = 0.5d0,epsilon = 1.0d-7,delta=0.1d0
    real(kind=8), allocatable :: xk(:),xtrial(:),faux(:),gaux(:),indices(:)
    integer, allocatable :: Idelta(:)
    real(kind=8) :: fxk,fxtrial,aux,opt_cond,sigma
    logical :: box

    ! COMMON SCALARS
    integer :: samples,q

    ! COMMON ARRAYS
    real(kind=8),   pointer :: t(:), y(:)

    ! COMMON BLOCKS
    common /integerData/ samples,q
    common /realVectorData/ t,y

    ! LOCAL SCALARS
    logical :: checkder
    integer :: hnnzmax,inform,jcnnzmax,m,n,nvparam
    real(kind=8) :: cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt,f,nlpsupn,snorm

    ! LOCAL ARRAYS
    character(len=80) :: specfnm,outputfnm,vparam(10)
    logical :: coded(11)
    logical,        pointer :: equatn(:),linear(:)
    real(kind=8),   pointer :: l(:),lambda(:),u(:),x(:)

    n = 4
    samples = 34
    q = 33

    allocate(t(samples),y(samples),x(n),xk(n-1),xtrial(n-1),l(n),u(n),&
    faux(samples),indices(samples),Idelta(samples),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error in main program'
        stop
    end if
  
    call read_data()

    ! Coded subroutines

    coded(1:6)  = .true.  ! evalf, evalg, evalh, evalc, evaljac, evalhc
    coded(7:11) = .false. ! evalfc,evalgjac,evalgjacp,evalhl,evalhlp


    ! Upper bounds on the number of sparse-matrices non-null elements
    jcnnzmax = 10000
    hnnzmax  = 10000

    ! Checking derivatives?
    checkder = .false.

    ! Parameters setting
    epsfeas   =   1.0d-08
    epsopt    =   1.0d-08

    efstain   =   1.0d+20
    eostain   = - 1.0d+20

    efacc     = - 1.0d+20
    eoacc     = - 1.0d+20

    outputfnm = ''
    specfnm   = ''

    nvparam   = 1
    vparam(1) = 'ITERATIONS-OUTPUT-DETAIL 0' 

    !==============================================================================
    ! MAIN ALGORITHM
    !==============================================================================
    iter = 0

    xk = 1.0d0

    box = .false.

    if (box .eqv. .false.) then
        l(1:n)   = -1.0d+20
        u(1:n-1) = 1.0d+20; l(n) = 0.0d0
    else
        l(1:n-1) = 0.0d0;   l(n) = -1.0d+20 
        u(1:n-1) = 1.0d+20; l(n) = 0.0d0
    endif

    indices(:) = (/(i,i=1,samples)/)

    kflag = 2

    do i = 1, samples
        faux(i) = fi(xk,i,n)
    end do

    call DSORT(faux,indices,samples,kflag)

    call mount_Idelta(faux,xk,n,delta,indices,Idelta,m)

    print*, m

    CONTAINS

    !==============================================================================
    ! EXPORT RESULT TO PLOT
    !==============================================================================
    subroutine export(x,n)
        implicit none

        integer,        intent(in) :: n
        real(kind=8),   intent(in) :: x(n-1)

        Open(Unit = 10, File = "output/xstarovo.txt", ACCESS = "SEQUENTIAL")

        write(10,*) x(1)
        write(10,*) x(2)
        write(10,*) x(3)

        close(10)

    end subroutine export

    !==============================================================================
    !
    !==============================================================================
    subroutine mount_Idelta(f,x,n,delta,indices,Idelta,m)
        implicit none

        integer,        intent(in) :: n
        real(kind=8),   intent(in) :: delta,x(n),f(samples),indices(samples)
        integer,        intent(out) :: Idelta(samples),m
        integer :: samples,q,i
        real(kind=8) :: fq

        common /integerData/ samples,q

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
    !
    !==============================================================================
    subroutine mount_grad(dFmin,x,n,indices)
        implicit none

        integer,        intent(in) :: n,indices(q)
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(out) :: dFmin(n)
        integer :: i,j
        real(kind=8) :: gx

        ! COMMON SCALARS
        integer :: samples,p,q

        ! COMMON ARRAYS
        real(kind=8),   pointer :: t(:), y(:)

        common /integerData/ samples,q
        common /realVectorData/ t,y

        dFmin(1:n) = 0.0d0
                
        do j = 1, q
            gx = (model(x,indices(j),n) - y(indices(j)))
            dFmin(1:n) = dFmin(1:n) + gx * (/((t(indices(j)) - t(samples))**i, i = 1, n)/)
        enddo

    end subroutine mount_grad

    function fi(x,i,n) result (res)
        implicit none

        integer,        intent(in) :: n,i
        real(kind=8),   intent(in) :: x(n)
        real(kind=8) :: res

        ! COMMON SCALARS
        integer :: samples,q

        ! COMMON ARRAYS
        real(kind=8),   pointer :: t(:), y(:)

        common /integerData/ samples,q
        common /realVectorData/ t,y
        
        res = model(x,i,n) - y(i)
        res = 0.5d0 * (res**2)

    end function fi

    function model(x,i,n) result(res)
        implicit none 

        integer,        intent(in) :: n,i
        real(kind=8),   intent(in) :: x(n-1)
        real(kind=8) :: res

        ! COMMON SCALARS
        integer :: samples,q

        ! COMMON ARRAYS
        real(kind=8),   pointer :: t(:), y(:)

        common /integerData/ samples,q
        common /realVectorData/ t,y

        res = (1.0d0 / x(3)) * exp(-x(3) * t(i)) * ((x(1) * t(i)) + (x(1) / x(3)) - x(2))
        res = 1.0d0 - exp(res - (x(2) * t(i)) - (x(1) / x(3)**2) + (x(2) / x(3))) 

    end function model

    !==============================================================================
    ! READ THE DATA CORRESPONDING TO THE NUMBER OF days DESIRED
    !==============================================================================
    subroutine read_data()
        implicit none

        ! COMMON SCALARS
        integer :: samples,q

        ! COMMON ARRAYS
        real(kind=8),   pointer :: t(:),y(:)

        ! COMMON BLOCKS
        common /integerData/ samples,q
        common /realVectorData/ t,y

        ! SCALARS
        integer :: i

        Open(Unit = 10, File = "output/zika.txt", ACCESS = "SEQUENTIAL")

        do i = 1, samples
            read(10,*) t(i), y(i)
        enddo

        close(10)
    end subroutine read_data
end Program main
