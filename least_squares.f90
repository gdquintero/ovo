program least_squares
    implicit none

    integer                     :: allocerr,i,k,info,lwork,nrhs,lda,ldb,m,n
	real(kind=8),   allocatable :: data(:),A(:,:),b(:),tempo(:),work(:)
    integer,        allocatable :: jpvt(:)
	real(kind=8)                :: rank,rcond = 1.0d-8

    n = 3
    Open(Unit = 200, File = "output/days.txt", ACCESS = "SEQUENTIAL")
    read(200,*) m
    close(200)

    lwork = max(min(m,n) + 3 * n + 1, 2 * min(m,n) + nrhs)
    nrhs = 1
    lda = max(1,m)
    ldb = max(1,m,n)

    allocate(data(m),tempo(m),A(lda,n),b(m-1),jpvt(n),work(max(1,lwork)),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error in main program'
        stop
    end if

    call read_data(m,data)

    tempo(1:m) = (/(i,i=1,m)/)

    do i = 1, lda - 1
        A(i,1:n) = (/((tempo(i) - tempo(m))**k,k = 1, n)/)
    enddo

    b(:) = data(1:m-1) - data(m)

    call dgelsy(m,n,nrhs,A,lda,b,ldb,jpvt,rcond,rank,work,lwork,info)

    call export(b(1:n),data(m),n)

    CONTAINS

    !==============================================================================
    ! EXPORT RESULT TO PLOT
    !==============================================================================
    subroutine export(xtrial,n)
        implicit none

        integer,        intent(in) :: n
        real(kind=8),   intent(in) :: xtrial(n),x1

        Open(Unit = 10, File = "output/xstarleastsquares.txt", ACCESS = "SEQUENTIAL")

        write(10,*) xtrial(1)
        write(10,*) xtrial(2)
        write(10,*) xtrial(3)

        close(10)

    end subroutine export

    !==============================================================================
    ! READ THE DATA CORRESPONDING TO THE NUMBER OF days DESIRED
    !==============================================================================
    subroutine read_data(days,y)
        implicit none

        ! SCALARS
        integer, intent(in) :: days
        integer :: i, j = 1
        real(kind=8) :: dat

        ! ARRAYS
        real(kind=8), intent(out) :: y(days)

        Open(Unit = 20, File = "output/data.txt", ACCESS = "SEQUENTIAL")

        do i = 1, days
            read(20,*) dat

            if (i .ge. 1) then
                y(j) = dat
                j = j + 1
            endif
        enddo

        close(20)
    end subroutine read_data
end program least_squares