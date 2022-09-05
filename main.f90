Program main
    use sort

    implicit none 
    
    integer :: allocerr,iter,intIter,predictedDays,i,j,r,i4_choose,vk,dimImin
    integer, parameter :: n = 3,maxIter = 1000,maxIntIter = 1000,nzMaxImin = 10
    real, parameter :: alpha = 0.5d0,epsilon = 1.0d-7
    real(kind=8), allocatable :: sk(:),xk(:),xtrial(:),dFmin(:),faux(:),gaux(:),indices(:)
    real(kind=8) :: lambda,fxk,fxtrial,aux,optCond,sigma

    ! COMMON SCALARS
    integer :: m,q

    ! COMMON ARRAYS
    real(kind=8),   pointer :: t(:), y(:)

    ! COMMON BLOCKS
    common /integerData/ m,q
    common /realVectorData/ t,y

    m = 34; q = 1

    allocate(t(m),y(m),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error in main program'
        stop
    end if
  
    call readData()

    print*, y

    CONTAINS

    !==============================================================================
    ! EXPORT RESULT TO PLOT
    !==============================================================================
    subroutine export(x,n,x1,ntrain,nval)
        implicit none

        integer,        intent(in) :: n,ntrain,nval
        real(kind=8),   intent(in) :: x(n),x1
        real(kind=8) :: aux
        integer :: i,j

        Open(Unit = 100, File = "output/xstarlovo.txt", ACCESS = "SEQUENTIAL")
        Open(Unit = 110, File = "output/training.txt", ACCESS = "SEQUENTIAL")
        Open(Unit = 120, File = "output/validation.txt", ACCESS = "SEQUENTIAL")
        Open(Unit = 130, File = "output/data.txt", ACCESS = "SEQUENTIAL")
        Open(Unit = 140, File = "output/data2.txt", ACCESS = "SEQUENTIAL")

        write(100,*) x1
        write(100,*) x(1)
        write(100,*) x(2)
        write(100,*) x(3)

        do i = 1, ntrain
            read(130,*) aux
            write(110,*) i, aux
        enddo

        j = i

        do i = 1, nval
            read(140,*) aux
            write(120,*) j, aux
            j = j + 1
        enddo
        
        close(100)

    end subroutine export

    !==============================================================================
    !
    !==============================================================================

    subroutine mountGrad(dFmin,x,n,indices)
        implicit none

        integer,        intent(in) :: n,indices(q)
        real(kind=8),   intent(in) :: x(n)
        real(kind=8),   intent(out) :: dFmin(n)
        integer :: i,j
        real(kind=8) :: gx

        ! COMMON SCALARS
        integer :: m,p,q

        ! COMMON ARRAYS
        real(kind=8),   pointer :: t(:), y(:)

        common /integerData/ m,q
        common /realVectorData/ t,y

        dFmin(1:n) = 0.0d0
                
        do j = 1, q
            gx = (model(x,indices(j),n) - y(indices(j)))
            dFmin(1:n) = dFmin(1:n) + gx * (/((t(indices(j)) - t(m))**i, i = 1, n)/)
        enddo

    end subroutine mountGrad

    function fi(x,i,n) result (res)
        implicit none

        integer,        intent(in) :: n,i
        real(kind=8),   intent(in) :: x(n)
        real(kind=8) :: res
        integer :: k

        ! COMMON SCALARS
        integer :: m,q

        ! COMMON ARRAYS
        real(kind=8),   pointer :: t(:), y(:)

        common /integerData/ m,q
        common /realVectorData/ t,y
        
        res = model(x,i,n) - y(i)
        res = 0.5d0 * (res**2)

    end function fi

    function model(x,i,n) result(res)
        implicit none 

        integer,        intent(in) :: n,i
        real(kind=8),   intent(in) :: x(n)
        real(kind=8) :: res
        integer :: k

        ! COMMON SCALARS
        integer :: m,q

        ! COMMON ARRAYS
        real(kind=8),   pointer :: t(:), y(:)

        common /integerData/ m,q
        common /realVectorData/ t,y

        res = dot_product(x,(/((t(i) - t(m))**k, k = 1, n)/))
        res = res + y(m)

    end function model

    !==============================================================================
    ! READ THE DATA CORRESPONDING TO THE NUMBER OF days DESIRED
    !==============================================================================
    subroutine readData()
        implicit none

        ! COMMON SCALARS
        integer :: m,q

        ! COMMON ARRAYS
        real(kind=8),   pointer :: t(:),y(:)

        ! COMMON BLOCKS
        common /integerData/ m,q
        common /realVectorData/ t,y

        ! SCALARS
        integer :: i

        Open(Unit = 10, File = "output/zika.txt", ACCESS = "SEQUENTIAL")

        do i = 1, m
            read(10,*) t(i), y(i)
        enddo

        close(10)
    end subroutine readData
end Program main
