program helicity

!claculate helicty use Devore gauge in FV method
implicit none
integer nx, ny, nz,max_n,file_num ,file_num_min,file_num_max,ios,error_estimation
CHARACTER(LEN=200) :: file_path
real(kind=8) h_v,h_j,energy,free_energy,b_flux,Kp,delta_x,pi,energy_s,div_error,C_erg 

!================================parameter that can be modified  =======================
parameter(file_num_min=0,file_num_max=79,error_estimation=0)  !0.146515 for 1 arcsec,2*3.6795d7

real(kind=8),allocatable :: fp(:,:,:)   !,E,E_free,H_J_3,H_V_3
real(kind=8),allocatable  :: Bx(:,:,:),By(:,:,:), Bz(:,:,:),Bp_x(:,:,:), Bp_y(:,:,:), Bp_z(:,:,:) 
real(kind=8), allocatable :: Ax(:,:,:),Ay(:,:,:), Az(:,:,:),Ap_x(:,:,:), Ap_y(:,:,:), Ap_z(:,:,:)
CHARACTER(LEN=200) :: file_bd,file_nlfff,file_out,file_potential,flie_Ap,file_A,formatted_num
character(6) BCtype
parameter(pi=3.1415926D0)
real(kind=8),allocatable :: bd_ax(:,:),bd_bx(:,:),bd_ay(:,:),bd_by(:,:),bd_az(:,:),bd_bz(:,:)
real(kind=8) :: start_time, end_time, elapsed_time
!parameters for energy accuracy
real(kind=8), allocatable :: grad_Bx(:,:,:),grad_By(:,:,:),grad_Bz(:,:,:),Bpn_y(:,:,:),Bpn_x(:,:,:),Bpn_z(:,:,:)
real(kind=8), allocatable :: Bj_x(:,:,:),Bj_y(:,:,:),Bj_z(:,:,:),Bjn_y(:,:,:),Bjn_x(:,:,:),Bjn_z(:,:,:)

! time 
call cpu_time(start_time)  
!====================================path===============================================
!file_path = '/share/home/guoyang1/zhouchang/helicity_tmf/11261/datamr/out_data/'
file_path = '/share/home/guoyang1/zhouchang/gjh/datamr/out_data/'

file_out = trim(file_path) // 'out.binary'

OPEN(UNIT=5, FILE=file_out, FORM='UNFORMATTED', STATUS='replace',access='stream', ACTION='WRITE')

! for many file
do file_num=file_num_min,file_num_max
    write(formatted_num, '(I4.4)') file_num
    
file_bd = trim(file_path) // 'bd/bd_' // trim(formatted_num) // '.binary'
file_nlfff = trim(file_path) // 'nlfff/nlfff_'// TRIM(formatted_num) //  '.binary'

OPEN(UNIT=4,FILE=file_nlfff,STATUS='OLD',ACTION='READ',FORM='UNFORMATTED',access='stream')
REWIND(4)
    READ(4) delta_x
    READ(4) nx,ny,nz
    !nx = 320
    !ny = 160
    !nz = 160
    if (file_num==file_num_min) then
        ALLOCATE(Bx(nx,ny,nz),By(nx,ny,nz),Bz(nx,ny,nz),fp(nx,ny,nz))
        ALLOCATE(Ax(nx,ny,nz),Ay(nx,ny,nz),Az(nx,ny,nz))
        ALLOCATE(Bp_x(nx,ny,nz),Bp_y(nx,ny,nz),Bp_z(nx,ny,nz))
        ALLOCATE(Ap_x(nx,ny,nz),Ap_y(nx,ny,nz),Ap_z(nx,ny,nz))
        ALLOCATE(bd_ax(ny,nz),bd_bx(ny,nz),bd_ay(nx,nz),bd_by(nx,nz),bd_az(nx,ny),bd_bz(nx,ny))
    endif
	READ(4) Bx,By,Bz
CLOSE(4)

if (error_estimation==1) then
if (file_num==file_num_min) then
ALLOCATE(grad_Bx(nx,ny,nz),grad_By(nx,ny,nz),grad_Bz(nx,ny,nz),Bpn_x(nx,ny,nz),Bpn_y(nx,ny,nz),Bpn_z(nx,ny,nz))
ALLOCATE(Bj_x(nx,ny,nz),Bj_y(nx,ny,nz),Bj_z(nx,ny,nz),Bjn_x(nx,ny,nz),Bjn_y(nx,ny,nz),Bjn_z(nx,ny,nz))
endif
endif

OPEN(UNIT=3,FILE=file_bd,STATUS='OLD',ACTION='READ',Access='stream',FORM='UNFORMATTED',POSITION='REWIND')
	READ(3) bd_ax,bd_bx,bd_ay,bd_by,bd_az,bd_bz
CLOSE(3)

! Setting the type of the boundary conditions on each surface
! Dirichlet boundary condition or Neumann boundary condition
BCtype = 'NNNNNN'
fp = 0.0D0


call poisson_function(BCtype,nx-1,ny-1,nz-1,bd_ax,bd_bx,bd_ay,bd_by,bd_az,bd_bz,fp) 
call gradient(fp,1,Bp_x)
call gradient(fp,2,Bp_y)
call gradient(fp,3,Bp_z)
!let the normal component of Bp equal to B  ,valori 2013,p3
Bp_x(1,:,:) = Bx(1,:,:)
Bp_x(nx,:,:) = Bx(nx,:,:)
Bp_y(:,1,:) = By(:,1,:)
Bp_y(:,ny,:) = By(:,ny,:)
Bp_z(:,:,1) = Bz(:,:,1)
Bp_z(:,:,nz) = Bz(:,:,nz)
Bj_x = Bx - Bp_x
Bj_y = By - Bp_y
Bj_z = Bz - Bp_z

!calculate A and Ap
call compute_A(Bp_x, Bp_y, Bp_z,Ap_x,Ap_y,Ap_z,delta_x)
call compute_A(Bx,By,Bz,Ax, Ay, Az,delta_x)

!energy
fp = (bx**2+by**2+bz**2)/(8*pi)
call volume_integral(fp,delta_x,energy)
!free energy
fp = (Bj_x**2+Bj_y**2+Bj_z**2)/(8*pi)
call volume_integral(fp,delta_x,free_energy)
!Hv
fp = (Ax+Ap_x)*(Bx-Bp_x)+(Ay+Ap_y)*(By-Bp_y)+(Az+Ap_z)*(Bz-Bp_z)
call volume_integral(fp,delta_x,h_v)
!Hj
fp = (Ax-Ap_x)*(Bx-Bp_x)+(Ay-Ap_y)*(By-Bp_y)+(Az-Ap_z)*(Bz-Bp_z)
call volume_integral(fp,delta_x,h_j)

b_flux = sum(abs(Bz(:,:,1)))*delta_x**2

!write useful parameters
write(5) energy
write(5) free_energy
write(5) h_v
write(5) h_j
write(5) b_flux



if (error_estimation==1) then
    call gradient(Bp_x,1,grad_Bx)
    call gradient(Bp_y,2,grad_By)
    call gradient(Bp_z,3,grad_Bz)
    fp = -(grad_Bx + grad_By +grad_Bz)
    bd_ax = 0
    bd_bx = 0
    bd_ay = 0
    bd_by = 0
    bd_az = 0
    bd_bz = 0
    call poisson_function(BCtype,nx-1,ny-1,nz-1,bd_ax,bd_bx,bd_ay,bd_by,bd_az,bd_bz,fp)
    call gradient(fp,1,Bpn_x)
    call gradient(fp,2,Bpn_y)
    call gradient(fp,3,Bpn_z)
    !valori 2013, let normal component to 0
    Bpn_x(1,:,:) = 0
    Bpn_x(nx,:,:) = 0
    Bpn_y(:,1,:) = 0
    Bpn_y(:,ny,:) = 0
    Bpn_z(:,:,1) = 0
    Bpn_z(:,:,nz) = 0
    
    !Bjn
    call gradient(Bj_x,1,grad_Bx)
    call gradient(Bj_y,2,grad_By)
    call gradient(Bj_z,3,grad_Bz)
    fp = -(grad_Bx + grad_By +grad_Bz)
    bd_ax = 0
    bd_bx = 0
    bd_ay = 0
    bd_by = 0
    bd_az = 0
    bd_bz = 0
    call poisson_function(BCtype,nx-1,ny-1,nz-1,bd_ax,bd_bx,bd_ay,bd_by,bd_az,bd_bz,fp)
    call gradient(fp,1,Bjn_x)
    call gradient(fp,2,Bjn_y)
    call gradient(fp,3,Bjn_z)
    Bjn_x(1,:,:) = 0
    Bjn_x(nx,:,:) = 0
    Bjn_y(:,1,:) = 0
    Bjn_y(:,ny,:) = 0
    Bjn_z(:,:,1) = 0
    Bjn_z(:,:,nz) = 0

    fp  = (Bjn_x**2+Bjn_y**2+Bjn_z**2+Bpn_x**2+Bpn_y**2+Bpn_z**2)/(8*pi) + &
          ((Bp_x-Bpn_x)*Bjn_x+(Bp_y-Bpn_y)*Bjn_y+(Bp_z-Bpn_z)*Bjn_z+ &
          (Bp_x-Bpn_x)*Bpn_x+(Bp_y-Bpn_y)*Bpn_y+(Bp_z-Bpn_z)*Bpn_z+ &
          (Bp_x-Bpn_x)*(Bj_x-Bjn_x)+(Bp_y-Bpn_y)*(Bj_y-Bjn_y)+(Bp_z-Bpn_z)*(Bj_z-Bjn_z)+&
          Bpn_x*(Bj_x-Bjn_x)+Bpn_y*(Bj_y-Bjn_y)+Bpn_z*(Bj_z-Bjn_z)+&
          Bjn_x*(Bj_x-Bjn_x)+Bjn_y*(Bj_y-Bjn_y)+Bjn_z*(Bj_z-Bjn_z)+&
          Bjn_x*Bpn_x+Bjn_y*Bpn_y+Bjn_z*Bpn_z )/(4*pi)         
          
    call volume_integral(fp,delta_x,energy_s)
    div_error = energy_s/energy
    write(5) div_error
    
!error from nabla A
    call rot(Ax,Ay,Az,grad_Bx,grad_By,grad_Bz,delta_x)
   ! call rot(Ap_x,Ap_y,Ap_z,grad_Bx,grad_By,grad_Bz,delta_x)
    C_erg = sum(grad_Bx**2+grad_By**2+grad_Bz**2)/sum(Bx**2+By**2+Bz**2)
    write(5) C_erg 

    file_potential = trim(file_path)// 'potential/potential_'// TRIM(formatted_num) //  '.binary'
	OPEN(UNIT=4, FILE=file_potential, FORM='UNFORMATTED', STATUS='replace',access='stream', ACTION='WRITE')
    write(4) nx,ny,nz
	write(4) Bp_x
	write(4) Bp_y
	write(4) Bp_z
	CLOSE(4)
!save A and Ap

	file_A = trim(file_path) // 'A/A_'// TRIM(formatted_num) //  '.binary'
	OPEN(UNIT=4, FILE=file_A, FORM='UNFORMATTED', STATUS='replace',access='stream', ACTION='WRITE')
    write(4) nx,ny,nz
	write(4) Ax
	write(4) Ay
	write(4) Az
	CLOSE(4)
	file_A = trim(file_path) // 'A/Ap_'// TRIM(formatted_num) //  '.binary'
	OPEN(UNIT=4, FILE=file_A, FORM='UNFORMATTED', STATUS='replace',access='stream', ACTION='WRITE')
    write(4) nx,ny,nz
	write(4) Ap_x
	write(4) Ap_y
	write(4) Ap_z
	CLOSE(4)

    
endif


call cpu_time(end_time)  
elapsed_time = end_time - start_time
print *, 'running time: ', elapsed_time
enddo

CLOSE(5)
! Free MKL memory if any was allocated
! Success message to print if everything is OK
! End of the example code
contains
   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!INTERFACES FOR 3D CASE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

! zhouchang,2024.12.
subroutine gradient(f,axis,grad_f,delta_x)
    implicit none
    real(kind=8), dimension(:, :, :), intent(in)  :: f
    integer, intent(in)                          :: axis
    real(kind=8), dimension(:, :, :), intent(out) :: grad_f
    real(kind=8), intent(in),optional :: delta_x
    real(kind=8) :: h
    integer :: nx, ny, nz

    nx = size(f, 1)   
    ny = size(f, 2)  
    nz = size(f, 3)   
    
    if (present(delta_x)) then
        h = delta_x
    else
        h = 1.0d0
    end if
    
    select case (axis)
        case (1)
            grad_f(2:nx-1, :, :) = (f(3:nx, :, :) - f(1:nx-2, :, :)) * 0.5d0/h
            grad_f(1, :, :) = (-3.0d0*f(1,:,:) + 4.0d0*f(2,:,:) - f(3,:,:)) * 0.5d0/h
            grad_f(nx, :, :) = (3.0d0*f(nx,:,:) - 4.0d0*f(nx-1,:,:) + f(nx-2,:,:)) * 0.5d0/h
        case (2)
            grad_f(:, 2:ny-1, :) = (f(:, 3:ny, :) - f(:, 1:ny-2, :)) * 0.5d0/h
            grad_f(:, 1, :) = (-3.0d0*f(:,1,:) + 4.0d0*f(:,2,:) - f(:,3,:)) * 0.5d0/h
            grad_f(:, ny, :) = (3.0d0*f(:,ny,:) - 4.0d0*f(:,ny-1,:) + f(:,ny-2,:)) * 0.5d0/h
        case (3)
            grad_f(:, :, 2:nz-1) = (f(:, :, 3:nz) - f(:, :, 1:nz-2)) * 0.5d0/h
            grad_f(:, :, 1) = (-3.0d0*f(:,:,1) + 4.0d0*f(:,:,2) - f(:,:,3)) * 0.5d0/h
            grad_f(:, :, nz) = (3.0d0*f(:,:,nz) - 4.0d0*f(:,:,nz-1) + f(:,:,nz-2)) * 0.5d0/h
        case default
            print *, "错误：无效的轴向值 (axis=", axis, ")"
            print *, "允许的值：1(X), 2(Y), 3(Z)"
    end select
end subroutine gradient

subroutine rot(Ax, Ay, Az, Bx, By, Bz,delta_x)
    implicit none
    real(kind=8), dimension(:, :, :), intent(in)  :: Ax, Ay, Az  ! 输入向量场分量
    real(kind=8), dimension(:, :, :), intent(out) :: Bx, By, Bz   ! 输出旋度分量
    real(kind=8), intent(in) :: delta_x
    integer :: nx, ny, nz
    ! 临时数组存储偏导数
    real(kind=8), dimension(:, :, :), allocatable :: dAx_dy, dAx_dz
    real(kind=8), dimension(:, :, :), allocatable :: dAy_dx, dAy_dz
    real(kind=8), dimension(:, :, :), allocatable :: dAz_dx, dAz_dy
    
    ! 获取数组维度
    nx = size(Ax, 1)
    ny = size(Ax, 2)
    nz = size(Ax, 3)
    
    ! 验证所有输入数组尺寸一致
    if (size(Ay,1) /= nx .or. size(Ay,2) /= ny .or. size(Ay,3) /= nz .or. &
        size(Az,1) /= nx .or. size(Az,2) /= ny .or. size(Az,3) /= nz) then
        print *, "错误：输入向量场分量尺寸不一致！"
        return
    end if
    
    ! 验证输出数组尺寸匹配
    if (size(Bx,1) /= nx .or. size(Bx,2) /= ny .or. size(Bx,3) /= nz .or. &
        size(By,1) /= nx .or. size(By,2) /= ny .or. size(By,3) /= nz .or. &
        size(Bz,1) /= nx .or. size(Bz,2) /= ny .or. size(Bz,3) /= nz) then
        print *, "错误：输出旋度分量尺寸不匹配输入！"
        return
    end if
    
    ! 分配临时数组存储偏导数
    allocate(dAx_dy(nx, ny, nz), dAx_dz(nx, ny, nz))
    allocate(dAy_dx(nx, ny, nz), dAy_dz(nx, ny, nz))
    allocate(dAz_dx(nx, ny, nz), dAz_dy(nx, ny, nz))
    
    ! 计算所有需要的偏导数
    call gradient(Ax,2,dAx_dy)  ! ∂Ax/∂y
    call gradient(Ax,3,dAx_dz)  ! ∂Ax/∂z
    
    call gradient(Ay,1,dAy_dx)  ! ∂Ay/∂x
    call gradient(Ay,3,dAy_dz)  ! ∂Ay/∂z
    
    call gradient(Az,1,dAz_dx)  ! ∂Az/∂x
    call gradient(Az,2,dAz_dy)  ! ∂Az/∂y
    
    ! 计算旋度的各个分量
    ! Bx = ∂Az/∂y - ∂Ay/∂z
    Bx = (dAz_dy - dAy_dz)/delta_x
    
    ! By = ∂Ax/∂z - ∂Az/∂x
    By = (dAx_dz - dAz_dx)/delta_x
    
    ! Bz = ∂Ay/∂x - ∂Ax/∂y
    Bz = (dAy_dx - dAx_dy)/delta_x
    
    ! 清理临时数组
    deallocate(dAx_dy, dAx_dz, dAy_dx, dAy_dz, dAz_dx, dAz_dy)
end subroutine rot

subroutine compute_A(Bx,By,Bz,Ax, Ay, Az,delta_x)
	IMPLICIT NONE  
    INTEGER :: i,j,k,nx,ny,nz
    real(kind=8):: delta_x
    real(kind=8), dimension(:, :, :),intent(in) :: Bx,By,Bz
    real(kind=8), dimension(:, :, :),intent(out):: Ax,Ay,Az
    real(kind=8), dimension(:, :),allocatable   :: b_x, b_y
    
    nx = size(Bx, 1)
    ny = size(Bx, 2)
    nz = size(Bx, 3)
    allocate(b_x(nx, ny),b_y(nx, ny))
    
    Az = 0
!Az is 0 in DV gauge,Ax,Ay have two part integral,
!integral method is trapezoidal method o(h2)
    do i=1,nx
	   do j=1,ny
		   do k=nz,1,-1
		    If (k==nz) THEN
				Ax(i,j,k) = 0
				Ay(i,j,k) = 0
			ELSE if(k==nz-1) THEN
				Ax(i,j,k) = (By(i,j,nz-1)+By(i,j,nz))/2
				Ay(i,j,k) = (Bx(i,j,nz-1)+Bx(i,j,nz))/2
            else
                Ax(i,j,k) = Ax(i,j,k+2)+2*By(i,j,k+1)
                Ay(i,j,k) = Ay(i,j,k+2)+2*Bx(i,j,k+1) 
            END If
		   end do
	   end do
    end do
	do i=1,nx
	   do j=1,ny
		    if (j==1) THEN
				b_x(i,j) = 0
			ELSE if(j==2) then
				b_x(i,j) = (Bz(i,1,nz)+Bz(i,2,nz))/2
            else
                b_x(i,j) = b_x(i,j-2)+2*Bz(i,j-1,nz)  !valori 2012 eq(37)
			END If
			If (i ==1) THEN
				b_y(i,j) = 0
            else if(i==2) then
                b_y(i,j) =  (Bz(1,j,nz)+Bz(2,j,nz))/2
			ELSE
				b_y(i,j) =  b_y(i-2,j)+2*Bz(i-1,j,nz)
            END If
		   end do
	   end do
	do k=1,nz
        Ax(:,:,k) = -Ax(:,:,k)-0.5*b_x
        Ay(:,:,k) = Ay(:,:,k)+0.5*b_y
    end do
    Ax = delta_x*Ax
    Ay = delta_x*Ay
    
END subroutine compute_A
    

subroutine poisson_function(BCtype,nx,ny,nz,bd_ax,bd_bx,bd_ay,bd_by,bd_az,bd_bz,fp)
    use mkl_poisson
    implicit none
    integer,INTENT(IN) :: nx, ny, nz
    real(kind=8),INTENT(IN):: bd_ax(ny+1,nz+1), bd_bx(ny+1,nz+1)
    real(kind=8),INTENT(IN):: bd_ay(nx+1,nz+1), bd_by(nx+1,nz+1)
    real(kind=8),INTENT(IN):: bd_az(nx+1,ny+1), bd_bz(nx+1,ny+1)
    character(6),INTENT(IN):: BCtype
    real(kind=8), INTENT(inout) :: fp(nx+1,ny+1,nz+1)
    integer stat,i,j,k,ios,max_n
    integer ipar(128)
    real(kind=8) domain_ax, domain_bx, domain_ay, domain_by, domain_az, domain_bz, lx, ly, lz, hx, hy, hz,pi
    real(kind=8) dpar(13*(nx+ny)/2+9)
    parameter(pi=3.14159265358979324D0)
    CHARACTER(LEN=10) :: status
    real(kind=8) q
    type(DFTI_DESCRIPTOR), pointer :: xhandle, yhandle
    
    ! Defining the parallelepiped domain 0<x<1, 0<y<1, 0<z<1 for 3D Poisson Solver
    domain_ax=0.0D0
    domain_bx= real(nx)
    domain_ay=0.0D0
    domain_by= real(ny)
    domain_az=0.0D0
    domain_bz= real(nz)

    ipar = 0.0D0
    !*******************************************************************************
    ! Setting the coefficient q to 0.
    ! Note that this is the way to use Helmholtz Solver to solve Poisson problem!
    !*******************************************************************************
    q=0.0D0
    ! Initializing simple data structures of Poisson Library for 3D Poisson Solver
    call D_INIT_HELMHOLTZ_3D(domain_ax, domain_bx, domain_ay, domain_by, domain_az, &
                                domain_bz, nx, ny, nz, BCtype, q, ipar, dpar,stat)

    call D_COMMIT_HELMHOLTZ_3D(fp, bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, xhandle,&
                                                                yhandle, ipar, dpar, stat)

    call d_Helmholtz_3D(fp, bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, xhandle,&
                                                                yhandle, ipar, dpar, stat)

    ! Cleaning the memory used by xhandle and yhandle
    call free_Helmholtz_3D(xhandle, yhandle, ipar, stat)

END SUBROUTINE poisson_function

    
subroutine volume_integral(array,delta_x,sum_array)
    implicit none
    integer :: i, j, k, nx, ny, nz
    real(kind=8), dimension(:,:,:), intent(in) :: array
    real(kind=8), intent(in) :: delta_x
    real(kind=8), intent(out) :: sum_array
    real(kind=8) :: e1, e2, e3

    ! 获取数组的维度
    nx = size(array, 1)
    ny = size(array, 2)
    nz = size(array, 3)

    ! 初始化 e1, e2, e3
    e1 = 0.0
    e2 = 0.0
    e3 = 0.0

    ! 计算 e1: 8个角点
    e1 = 0.125 * (array(1, 1, 1) + array(1, 1, nz) + array(1, ny, 1) + array(1, ny, nz) + &
                  array(nx, 1, 1) + array(nx, 1, nz) + array(nx, ny, 1) + array(nx, ny, nz))

    ! 计算 e2: 12条边
    do i = 1, nz-2
        e2 = e2 + 0.25 * (array(1, 1, i+1) + array(1, ny, i+1) + array(nx, 1, i+1) + array(nx, ny, i+1))
    end do
    do j = 1, ny-2
        e2 = e2 + 0.25 * (array(1, j+1, 1) + array(1, j+1, nz) + array(nx, j+1, 1) + array(nx, j+1, nz))
    end do
    do k = 1, nx-2
        e2 = e2 + 0.25 * (array(k+1, 1, 1) + array(k+1, 1, nz) + array(k+1, ny, 1) + array(k+1, ny, nz))
    end do

    ! 计算 e3: 6个面
    do j = 1, ny-2
        do k = 1, nz-2
            e3 = e3 + 0.5 * (array(1, j+1, k+1) + array(nx, j+1, k+1))
        end do
    end do
    do i = 1, nx-2
        do k = 1, nz-2
            e3 = e3 + 0.5 * (array(i+1, 1, k+1) + array(i+1, ny, k+1))
        end do
    end do
    do i = 1, nx-2
        do j = 1, ny-2
            e3 = e3 + 0.5 * (array(i+1, j+1, 1) + array(i+1, j+1, nz))
        end do
    end do

    ! 计算内部单元
    sum_array = (e1 + e2 + e3 + sum(array(2:nx-1, 2:ny-1, 2:nz-1))) * delta_x*delta_x*delta_x
end subroutine volume_integral

end program helicity