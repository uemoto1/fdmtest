program main
    use omp_lib
    ! use mpi ! MPI
    implicit none
    integer, parameter :: nprocx = 1 ! MPI procs
    integer, parameter :: nprocy = 1 ! MPI procs
    integer, parameter :: nprocz = 1 ! MPI procs
    integer, parameter :: nx = 64
    integer, parameter :: ny = 64
    integer, parameter :: nz = 64
    integer, parameter :: nloop = 4
    real(8), parameter :: dx = 0.01d0
    real(8), parameter :: dy = 0.01d0
    real(8), parameter :: dz = 0.01d0
    integer, parameter :: nd = 4
    real(8), parameter :: lap_coeff(0:nd) = &
        & (/-205d0/72d0, 8d0/5d0, -1d0/5d0, 8d0/315d0, -1d0/560d0/)
    real(8), parameter :: pi = 4.0d0 * atan(1.0d0)

    ! Floating point operation per single grid point
    real(8), parameter :: flop1 = 2 * (3 * (2 * nd + 1) - 2) - 1

    real(8) :: lap_x(0:nd), lap_y(0:nd), lap_z(0:nd), lap0
    real(8), allocatable :: f(:, :, :)
    real(8), allocatable :: g(:, :, :)
    integer :: i, ix, iy, iz
    real(8) :: qx, qy, qz
    real(8) :: x, y, z, t1, t2
    real(8) :: flops, g_ref, g_err, g_max

    ! MPI variables
    integer :: irank, nproc, ierr
    integer, allocatable :: itbl_rank(:, :, :), itbl_xnum(:), itbl_ynum(:), itbl_znum(:)
    real(8), allocatable :: buf_xy(:, :, :), buf_yz(:, :, :), buf_zx(:, :, :)

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr) 
    call MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)

    if (nproc .ne. nprocx * nproxy * nprocz) &
        & stop " nprocx * nprocy * nprocz is not equal nproc!"

    allocate(itbl_rank(nprocx, nprocy, nprocz))
    allocate(itbl_xnum(0:nproc-1))
    allocate(itbl_ynum(0:nproc-1))
    allocate(itbl_znum(0:nproc-1))
    ! allocate(buf_x(nx, ny, nd))
    ! allocate(buf_y(nx, nd, nz))
    ! allocate(buf_z(nd, ny, nz))
    
    i = 0
    do iz = 1, nprocz
        do iy = 1, nprocy
            do ix = 1, nprocx
                itbl_rank(ix, iy, iz) = i
                itbl_xnum(i) = ix
                itbl_ynum(i) = iy
                itbl_znum(i) = iz
                i = i + 1
            end do
        end do
    end do

    lap_x(0:nd) = lap_coeff(0:nd) * (1.0d0 / dx ** 2)
    lap_y(0:nd) = lap_coeff(0:nd) * (1.0d0 / dy ** 2)
    lap_z(0:nd) = lap_coeff(0:nd) * (1.0d0 / dz ** 2)
    lap0 = lap_x(0) + lap_y(0) + lap_z(0)

    allocate(f(1-nd:nx+nd, 1-nd:ny+nd, 1-nd:nz+nd))
    allocate(g(1:nx, 1:ny, 1:nz))

    g(1:nx, 1:ny, 1:nz) = 0d0 
    f(1-nd:nx+nd, 1-nd:ny+nd, 1-nd:nz+nd) = 0d0 

    ! Initial setting
    qx = 2.0d0 * pi / (nx * nprocx * dx) * 1
    qy = 2.0d0 * pi / (ny * nprocy * dy) * 2
    qz = 2.0d0 * pi / (nz * nprocz * dz) * 3
    do iz = 1-nd, nz+nd
        do iy = 1-nd, ny+nd
            do ix = 1-nd, nx+nd
                x = ix * dx
                y = iy * dy
                z = iz * dz
                f(ix, iy, iz) = sin(qx*x) * sin(qy*y) * sin(qz*z)
            end do
        end do
    end do

    ! Stencil calculation
    t1 =  omp_get_wtime()
    do i = 1, nloop

        

        do iz = 1, nz
            do iy = 1, ny
                do ix = 1, nx
                    g(ix, iy, iz) = lap0 * f(ix, iy, iz) &
                        & + lap_x(4) * f(ix-4, iy, iz) &
                        & + lap_x(3) * f(ix-3, iy, iz) &
                        & + lap_x(2) * f(ix-2, iy, iz) &
                        & + lap_x(1) * f(ix-1, iy, iz) &
                        & + lap_x(1) * f(ix+1, iy, iz) &
                        & + lap_x(2) * f(ix+2, iy, iz) &
                        & + lap_x(3) * f(ix+3, iy, iz) &
                        & + lap_x(4) * f(ix+4, iy, iz) &
                        & + lap_y(4) * f(ix, iy-4, iz) &
                        & + lap_y(3) * f(ix, iy-3, iz) &
                        & + lap_y(2) * f(ix, iy-2, iz) &
                        & + lap_y(1) * f(ix, iy-1, iz) &
                        & + lap_y(1) * f(ix, iy+1, iz) &
                        & + lap_y(2) * f(ix, iy+2, iz) &
                        & + lap_y(3) * f(ix, iy+3, iz) &
                        & + lap_y(4) * f(ix, iy+4, iz) &
                        & + lap_z(4) * f(ix, iy, iz-4) &
                        & + lap_z(3) * f(ix, iy, iz-3) &
                        & + lap_z(2) * f(ix, iy, iz-2) &
                        & + lap_z(1) * f(ix, iy, iz-1) &
                        & + lap_z(1) * f(ix, iy, iz+1) &
                        & + lap_z(2) * f(ix, iy, iz+2) &
                        & + lap_z(3) * f(ix, iy, iz+3) &
                        & + lap_z(4) * f(ix, iy, iz+4)
                end do
            end do
        end do
    end do
    t2 =  omp_get_wtime()

    ! Flops rating
    flops = dble(flop1) * dble(nx) * dble(ny) * dble(nz) * dble(nloop) / (t2 - t1)

    write(*,*) "Size:", nx, ny, nz
    write(*,*) "Repeat:", nloop
    write(*,*) "Total time [s]:", t2-t1
    write(*, *) "Performance [GFlop/s]:", flops * 1d-9

    ! Result error check
    g_err = 0.0d0
    g_max = 0.0d0
    do iz = 1, nz
        do iy = 1, ny
            do ix = 1, nx
                x = ix * dx
                y = iy * dy
                z = iz * dz
                g_ref = (-qx*qx-qy*qy-qz*qz) * sin(qx*x)*sin(qy*y)*sin(qz*z)
                g_err = g_err + abs(g(ix, iy, iz) - g_ref)
                g_max = max(g_max, abs(g(ix, iy, iz) - g_ref))
            end do
        end do
    end do

    write(*,*) "Sum(g):", sum(g)
    write(*,*) "AvgErr(g-g_ref):",  g_err / (nx * ny * nz)
    write(*,*) "MaxErr(g-g_ref):",  g_max

    stop
end program main

