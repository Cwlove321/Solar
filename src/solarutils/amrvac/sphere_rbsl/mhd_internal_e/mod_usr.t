module mod_usr   ! add filament simulation in spherical geometry
  use mod_mhd
  implicit none
  ! variables for boundary condition
  double precision, allocatable, save :: Bnl(:,:,:)
  double precision, allocatable, save :: density0(:)
  ! some global parameters
  double precision, save :: k_B,miu0,mass_H
  double precision, save :: grav_sun,radius_sun
  double precision, save :: h_trans,h_trans2,T_photo,T_corona
  double precision, save :: rho_photo,rho_tr_top,rho_tr_top2
  double precision, save :: k0,T0,const
  double precision, save :: densitymax,densitymin
  logical, save :: firstusrglobaldata=.true.
  double precision, save :: mu_density, dr, small_density0

  contains

  !==============================================================================
  ! Purpose: to include global parameters, set user methods, set coordinate 
  !          system and activate physics module.
  !==============================================================================
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_special_resistivity => specialres_usr
    usr_gravity         => gravity
    usr_source          => specialsource_usr
    !usr_get_dt          => special_getdt
    usr_refine_grid     => special_refine_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_transform_w    => transform_w_fill_plasma
    call set_coordinate_system("spherical")
    call mhd_activate()
  end subroutine usr_init

  !==============================================================================
  ! Purpose: to initialize user public parameters and reset global parameters.
  !          Input data are also read here.
  !==============================================================================
  subroutine initglobaldata_usr()
    use mod_global_parameters
    integer :: ix1
    double precision :: rr

    ! normalization unit in CGS Unit
    k_B = 1.3806d-16          ! erg*K^-1
    miu0 = 4.d0*dpi           ! Gauss^2 cm^2 dyne^-1
    mass_H = 1.67262d-24      ! g
    unit_length        = 6.955d10 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3
    unit_density       = 1.4d0*mass_H*unit_numberdensity               ! 2.341668000000000E-015 g*cm^-3
    unit_pressure      = 2.3d0*unit_numberdensity*k_B*unit_temperature ! 0.317538000000000 erg*cm^-3
    unit_magneticfield = dsqrt(miu0*unit_pressure)                     ! 1.99757357615242 Gauss
    unit_velocity      = unit_magneticfield/dsqrt(miu0*unit_density)   ! 1.16448846777562E007 cm/s = 116.45 km/s
    unit_time          = unit_length/unit_velocity                     ! 5972.5794 s = 99.543 min
    
    mu_density = 0.3d0
    small_density0 = 1.0d0    !The small density to set initial and boundary conditions 模拟暗条慢相可以设高

    grav_sun=-2.74d4*unit_length/unit_velocity**2     ! solar gravity
    radius_sun= 6.955d10/unit_length                  ! solar radius
    h_trans   = 3.5d8/unit_length         ! lower height of the transition region
    h_trans2  = 1.0d9/unit_length         ! upper height of the transition region
    T_photo   = 8.5d3/unit_temperature    ! temperature in the photosphere and chromosppere
    T_corona  = 1.0d6/unit_temperature    ! temperature in the corona
    rho_photo = 2.341668d-7/unit_density  ! density in the photosphere, which is 1.0d8 

    k0 = (T_corona-T_photo)/(h_trans2-h_trans)
    T0= T_photo-k0*h_trans    ! T = T_photo                   if     0.0  <= h < h_trans
                              ! T = k(h-h_trans) + T_photo    if h_trans  <= h < h_trans2
                              ! T = T_corona                  if h_trans2 <= h
                              ! so between h_trans and h_trans2, T = k*h - k*h_trans + T_photo. T0 is defined
                              ! as T_photo-k0*h_trans, and T = k*h + T0
    rho_tr_top=rho_photo*dexp(grav_sun*radius_sun**2/T_photo*&               
                     (1.d0/radius_sun-1.d0/(h_trans+radius_sun)))
    const = rho_tr_top/dexp(-(-grav_sun*radius_sun**2*(k0*radius_sun-T0)+&
            grav_sun*k0*radius_sun**2*(radius_sun+h_trans)*log(abs(radius_sun+h_trans))+&
            (-grav_sun*k0*radius_sun**2+(-k0*radius_sun+T0)**2)*(radius_sun+h_trans)*&
            log(abs(T0+k0*h_trans)))/((-k0*radius_sun+T0)**2*(radius_sun+h_trans)))
    rho_tr_top2 = const*dexp(-(-grav_sun*radius_sun**2*(k0*radius_sun-T0)+&
            grav_sun*k0*radius_sun**2*(radius_sun+h_trans2)*log(abs(radius_sun+h_trans2))+&
            (-grav_sun*k0*radius_sun**2+(-k0*radius_sun+T0)**2)*(radius_sun+h_trans2)*&
            log(abs(T0+k0*h_trans2)))/((-k0*radius_sun+T0)**2*(radius_sun+h_trans2)))

    allocate(density0(-1:domain_nx1+2))
    dr = (xprobmax1-xprobmin1)/domain_nx1
    do ix1 = -1, domain_nx1+2
      rr = xprobmin1 + (ix1-0.5d0)*dr
      if ((rr-radius_sun) <= h_trans) then
        density0(ix1)=rho_photo*dexp(grav_sun*radius_sun**2/T_photo*(1.d0/radius_sun-1.d0/rr))
      else
        if ((rr-radius_sun) <= h_trans2) then
          density0(ix1)=const*exp(-(-grav_sun*radius_sun**2*(k0*radius_sun-T0)+&
                grav_sun*k0*radius_sun**2*rr*log(abs(rr))+&
              (-grav_sun*k0*radius_sun**2+(-k0*radius_sun+T0)**2)*rr*&
                log(abs(T0+k0*(rr-radius_sun))))/((-k0*radius_sun+T0)**2*rr))
        else
          density0(ix1)=rho_tr_top2*dexp(grav_sun*radius_sun**2/T_corona*(1.d0/(radius_sun+h_trans2)-1.d0/rr))
        endif
      endif
      if (density0(ix1) < small_density0) density0(ix1) = small_density0
    end do
    densitymax=maxval(density0)
    densitymin=minval(density0)
    !print*, densitymax, densitymin
    
    ! prepare magnetogram at bottom
    if(firstusrglobaldata) then
      call init_b_data('../mf/nlfff_boundary/04field0_cutout.dat',unit_length,unit_magneticfield)
      firstusrglobaldata=.false.
    endif



  end subroutine initglobaldata_usr

  !==============================================================================
  ! Purpose: to read vector magnetic field from an input file
  !==============================================================================
  subroutine init_b_data(magnetogramname,qLunit,qBunit)
    use mod_global_parameters

    character(len=*), intent(in) :: magnetogramname
    double precision, intent(in) :: qLunit,qBunit
    integer, dimension(MPI_STATUS_SIZE) :: statuss
    integer :: file_handle,nx1,nx2
    logical :: aexist

    ! nx1,nx2 are numbers of cells for each direction
    ! Bnl(:,:,1), Bnl(:,:,2), and Bnl(:,:,3) are the radial, co-latitudinal, and 
    ! longitudinal components of the 2D vector magnetic field on the bottom.
    inquire(file=magnetogramname,exist=aexist)
    if(.not. aexist) then
      if(mype==0) write(*,'(2a)') "can not find file:",magnetogramname
      call mpistop("no input magnetogram----init_b_data")
    end if
    call MPI_FILE_OPEN(icomm,magnetogramname,MPI_MODE_RDONLY,MPI_INFO_NULL,&
                       file_handle,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,nx1,1,MPI_INTEGER,statuss,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,nx2,1,MPI_INTEGER,statuss,ierrmpi)

    allocate(Bnl(nx1,nx2,3))
    call MPI_FILE_READ_ALL(file_handle,Bnl,nx1*nx2*3,MPI_DOUBLE_PRECISION,&
                           statuss,ierrmpi)
    call MPI_FILE_CLOSE(file_handle,ierrmpi)
    ! declare and define global variables Lunit and Bunit to be your length unit in
    ! cm and magnetic strength unit in Gauss first
    Bnl=Bnl/qBunit
    
    if(mype==0) then
      print*,'zero-beta simulation using a bottom boundary from an observed Brad,Bthe, and Bphi'
      print*,'magnetogram of',nx1,'by',nx2,'pixels, which includes the ghost layers.'
    endif
  end subroutine init_b_data

  !==============================================================================
  ! Purpose: to initialize the initial condition
  !==============================================================================
  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: ix^D
    double precision :: temperature_local
    logical, save :: first=.true.

    if(first)then
      if(mype==0) then
      write(*,*)'Simulating 3D active region'
      print*,mhd_gravity
    endif
      first=.false.
    endif
   ! w(ixO^S,mom(1):mom(3))=zero
   ! w(ixO^S,mag(1:3))=w(ixO^S,mag(1:3))
   !w(ixO^S,mag(1):mag(3))=zero      ! to test the balance between the gravity and pressure gradient
   !w(ixO^S,rho_)=0.005d0*(^C&w(ixO^S,b^C_)**2+) + 0.005d0*0.125781**2    ! rho = c*(B^2+min(b)^2)
   !w(ixO^S,rho_)=rho_photo*dexp(grav_sun*radius_sun**2/T_corona*&               
   !               (1.d0/radius_sun-1.d0/(x(ixO^S,1)+radius_sun)))+&
   !                0.005d0*(^C&w(ixO^S,b^C_)**2+)                        ! rho = stratified density + c*B^2
   !w(ixO^S,rho_)=rho_photo*dexp(grav_sun*radius_sun**2/T_corona*&               
   !               (1.d0/radius_sun-1.d0/(x(ixO^S,1)+radius_sun)))              ! rho = stratified density

{do ix^DB=ixOmin^DB,ixOmax^DB\}

  if ((x(ix^D,1)-radius_sun) <= h_trans) then

    temperature_local = T_photo

    w(ix^D,rho_)=rho_photo*dexp(grav_sun*radius_sun**2/T_photo* &
         (1.d0/radius_sun-1.d0/(x(ix^D,1))))

  else

    if ((x(ix^D,1)-radius_sun) <= h_trans2) then

      temperature_local = k0*(x(ix^D,1)-radius_sun) + T_photo

      w(ix^D,rho_)=const*exp(-(-grav_sun*radius_sun**2*(k0*radius_sun-T0)+&
            grav_sun*k0*radius_sun**2*x(ix^D,1)*log(abs(x(ix^D,1)))+&
           (-grav_sun*k0*radius_sun**2+(-k0*radius_sun+T0)**2)*x(ix^D,1)*&
            log(abs(T0+k0*(x(ix^D,1)-radius_sun))))/((-k0*radius_sun+T0)**2*x(ix^D,1)))

    else

      temperature_local = T_corona

      w(ix^D,rho_)=rho_tr_top2*dexp(grav_sun*radius_sun**2/T_corona* &
           (1.d0/(radius_sun+h_trans2)-1.d0/x(ix^D,1)))

    endif

  endif

  ! finite-beta pressure
  w(ix^D,e_) = w(ix^D,rho_) * temperature_local/(mhd_gamma-1.d0)

{end do\}
    
  !  write(*,*)'use the atmosphere guo 2019'
    call add_filament_density_from_path(ixI^L, ixO^L, w, x)
  end subroutine initonegrid_usr

!==============================================================================
! Purpose:
!   Add cool dense filament/prominence material around ../pfss/rbsl/path.dat.
!
! Usage:
!   call add_filament_density_from_path(ixI^L, ixO^L, w, x)
!
! Assumption:
!   - path.dat is a binary stream file with 6000 points:
!       lon(deg), lat(deg), r(Rsun)
!   - mhd_internal_e = .true.
!   - w(e_) is internal energy density.
!
! Loading:
!   rho -> F*rho
!   mom -> F*mom
!   e_  unchanged
!
! Therefore p = (gamma-1)*e_ is unchanged and T = p/rho decreases.
!==============================================================================
subroutine add_filament_density_from_path(ixI^L, ixO^L, w, x)
  use mod_global_parameters
  implicit none

  integer, intent(in) :: ixI^L, ixO^L
  double precision, intent(inout) :: w(ixI^S,1:nw)
  double precision, intent(in)    :: x(ixI^S,1:ndim)

  integer, parameter :: np_file = 6000
  character(len=*), parameter :: pathfile = '../pfss/rbsl/path.dat'

  logical, save :: initialized = .false.
  double precision, allocatable, save :: x_axis(:,:)

  integer :: ix^D
  integer :: ip, imin, ios
  double precision :: lon_deg, lat_deg, rr
  double precision :: lon, lat
  double precision :: Xc(3)
  double precision :: d2, d2min, dist
  double precision :: q, qmin, qmax, qw
  double precision :: taper_s, taper_r
  double precision :: F
  double precision :: prom_Crho, prom_Fmax
  double precision :: prom_radius, prom_width

  !------------------------------------------------------------
  ! Parameters for filament material.
  ! First test: moderate density enhancement.
  !------------------------------------------------------------
  prom_Crho  = 5.0d0
  prom_Fmax  = 20.0d0

  prom_radius = 8.0d0 * 1.0d8 / unit_length   ! 8 Mm
  prom_width  = 3.0d0 * 1.0d8 / unit_length   ! 3 Mm

  ! Along-axis loading range.
  ! q = 0 at one end, q = 1 at the other end.
  ! This avoids loading too much mass near the two footpoints.
  qmin = 0.08d0
  qmax = 0.92d0
  qw   = 0.03d0

  !------------------------------------------------------------
  ! Read path.dat only once.
  !------------------------------------------------------------
  if (.not. initialized) then

    allocate(x_axis(np_file,3))

    open(unit=993, file=pathfile, action='read', status='old', &
         access='stream', iostat=ios)

    if (ios /= 0) then
      if (mype == 0) write(*,*) 'ERROR: cannot open filament path file: ', pathfile
      call mpistop('cannot open filament path file')
    endif

    do ip = 1, np_file

      read(993, iostat=ios) lon_deg, lat_deg, rr

      if (ios /= 0) then
        if (mype == 0) write(*,*) 'ERROR: failed to read path.dat at ip = ', ip
        call mpistop('failed to read filament path file')
      endif

      lon = lon_deg*dpi/180.0d0
      lat = lat_deg*dpi/180.0d0

      x_axis(ip,1) = rr*cos(lat)*cos(lon)
      x_axis(ip,2) = rr*cos(lat)*sin(lon)
      x_axis(ip,3) = rr*sin(lat)

    enddo

    close(993)

    initialized = .true.

    if (mype == 0) then
      write(*,*) 'Filament density loading enabled in initonegrid_usr.'
      write(*,*) 'Path file: ', pathfile
      write(*,*) 'Number of path points = ', np_file
      write(*,*) 'prom_Crho, prom_Fmax = ', prom_Crho, prom_Fmax
      write(*,*) 'prom_radius(Mm), prom_width(Mm) = ', &
                 prom_radius*unit_length/1.0d8, prom_width*unit_length/1.0d8
    endif

  endif

  !------------------------------------------------------------
  ! Add density around the path.
  !------------------------------------------------------------
{do ix^DB=ixOmin^DB,ixOmax^DB\}

    ! Spherical grid coordinate -> Cartesian coordinate.
    Xc(1) = x(ix^D,1) * sin(x(ix^D,2)) * cos(x(ix^D,3))
    Xc(2) = x(ix^D,1) * sin(x(ix^D,2)) * sin(x(ix^D,3))
    Xc(3) = x(ix^D,1) * cos(x(ix^D,2))

    d2min = 1.0d99
    imin  = 1

    do ip = 1, np_file

      d2 = (Xc(1)-x_axis(ip,1))**2 + &
           (Xc(2)-x_axis(ip,2))**2 + &
           (Xc(3)-x_axis(ip,3))**2

      if (d2 < d2min) then
        d2min = d2
        imin  = ip
      endif

    enddo

    dist = sqrt(d2min)

    ! Skip cells far from the filament tube.
    if (dist > prom_radius + 3.0d0*prom_width) cycle

    ! Along-axis coordinate normalized to [0,1].
    q = dble(imin-1)/dble(np_file-1)

    ! Smooth taper along the path.
    taper_s = 0.5d0*(tanh((q-qmin)/qw) - tanh((q-qmax)/qw))

    if (taper_s <= 0.0d0) cycle

    ! Smooth transverse tube profile.
    taper_r = 0.5d0*(1.0d0 + tanh((prom_radius-dist)/prom_width))

    F = 1.0d0 + prom_Crho*taper_s*taper_r

    if (F > prom_Fmax) F = prom_Fmax

    !----------------------------------------------------------
    ! Isobaric density loading.
    ! For internal-energy equation, do not change e_.
    !----------------------------------------------------------
    w(ix^D,rho_) = F*w(ix^D,rho_)

    w(ix^D,mom(1)) = F*w(ix^D,mom(1))
    w(ix^D,mom(2)) = F*w(ix^D,mom(2))
    w(ix^D,mom(3)) = F*w(ix^D,mom(3))

    ! Do NOT change w(ix^D,e_).

{end do\}

end subroutine add_filament_density_from_path

  subroutine transform_w_fill_plasma(ixI^L,ixO^L,nw_in,w_in,x,w_out)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, nw_in
    double precision, intent(in)  :: w_in(ixI^S,1:nw_in)
    double precision, intent(in)  :: x(ixI^S, 1:ndim)
    double precision, intent(out) :: w_out(ixI^S,1:nw)

    double precision :: res,Bf(ixI^S,1:ndir),xC(ixI^S,1:ndim)
    integer :: ix^D,na,idir,ixC^L
    w_out=0.d0
    ! only preserve magnetic field
    w_out(ixO^S,mag(1:3)) = w_in(ixO^S,5:7)
    w_out(ixO^S,mom(1:3)) = 0.0d0
    if(mhd_glm)  w_out(ixO^S,psi_)=0.d0
  end subroutine transform_w_fill_plasma
  !==============================================================================
  ! Purpose: to provide special boundary conditions set by users.
  !==============================================================================
  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters
    use mod_physics
    integer, intent(in) :: ixI^L, ixO^L, iB
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: drad(0:2),dthe,dphi,delt_x1,delt_x2
    double precision :: xlen^D,dxb^D,startpos^D
    integer :: ngl,ix^D,ixbc^D,af, ind_nx1

    ngl=nghostcells*(2**(refine_max_level-1))
    drad(0)=1.d0
    !drad(1)=drad(0)*qstretch_baselevel(1)
    !drad(2)=drad(1)*qstretch_baselevel(1)
    drad(1)=1.d0
    drad(2)=1.d0
    dthe=1.d0
    dphi=1.d0

    select case(iB)
    case(1)    ! r mininum surface
      w(ixO^S,rho_)=rho_photo*dexp(grav_sun*radius_sun**2/T_photo*(1.d0/radius_sun-1.d0/(x(ixO^S,1))))
      w(ixO^S,mom(1):mom(3))=zero
      w(ixO^S,e_) = w(ixO^S,rho_) * T_photo/(mhd_gamma-1.d0)
      dxb^D=dx(^D,refine_max_level);
      af=nint(dxlevel(2)/dxb2)
      ! The first layer (with the resolution of the level 1 grid) are assigned with the fixed boundary data
      ! The outer ghost layers are assigned with the following one-sided difference scheme
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        xlen^D=x({ix^D,},^D)-xprobmin^D+dble(ngl)*dxb^D;\
        if(af>=2) then
          call mpistop("Bottom boundary is not defined ")
        else
          ixbc^D=nint(xlen^D/dxb^D);\
          w(ix^D,mag(1))=Bnl(ixbc2,ixbc3,1)
          w(ix^D,mag(2))=Bnl(ixbc2,ixbc3,2)
          w(ix^D,mag(3))=Bnl(ixbc2,ixbc3,3)
        endif
      {end do\}
      ! af is the number of ghost cells that is filled with the fixed data. When the 
      ! bottom outer ghost layers are derived constant extrapolation, we assume that
      ! only 1 inner most layer are provided by the fixed data. In real observations,
      ! it requires that the lowest block should always have the highest resolution.
      af=1
      ! 2nd order accuacy constant value extrapolation
      do ix1=ixOmax1-af,ixOmin1,-1
        delt_x1 = x(ix1+1,ixOmin2,ixOmin3,1)-x(ix1+0,ixOmin2,ixOmin3,1)
        delt_x2 = x(ix1+2,ixOmin2,ixOmin3,1)-x(ix1+1,ixOmin2,ixOmin3,1)
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))=                   &
                (1.0d0/(delt_x2**2+2.0d0*delt_x1*delt_x2))*(-delt_x1**2* &
                     w(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))*    &
                    (x(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2*      &
                 sin(x(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*        &
                    dthe*dphi) + (delt_x1+delt_x2)**2*                   &
                     w(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))*    &
                    (x(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2*      &
                 sin(x(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*        &
                     dthe*dphi))/                                        &
                    (x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2*      &
                 sin(x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*        &
                     dthe*dphi)

        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2))=                   &
                (1.0d0/(delt_x2**2+2.0d0*delt_x1*delt_x2))*(-delt_x1**2* &
                     w(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2))*    &
                    (x(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*         &
                 sin(x(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*        &
                     drad(2)*dphi) + (delt_x1+delt_x2)**2*               &
                     w(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2))*    &
                    (x(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*         &
                 sin(x(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*        &
                     drad(1)*dphi))/                                     &
                    (x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*         &
                 sin(x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*        &
                     drad(0)*dphi)

        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3))=                   &
                (1.0d0/(delt_x2**2+2.0d0*delt_x1*delt_x2))*(-delt_x1**2* &
                     w(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3))*    &
                    (x(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*         &
                     drad(2)*dthe) + (delt_x1+delt_x2)**2*               &
                     w(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3))*    &
                    (x(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*         &
                     drad(1)*dthe))/                                     &
                    (x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*         &
                     drad(0)*dthe)
      enddo
    case(2)    ! r maxinum surface
      !w(ixO^S,rho_)=rho_tr_top2*dexp(grav_sun*radius_sun**2/T_corona*(1.d0/(radius_sun+h_trans2)-1.d0/x(ixO^S,1)))
      w(ixO^S,rho_)=w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
      w(ixO^S,e_)=w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)
      w(ixO^S,mom(1))=w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))
      w(ixO^S,mom(2))=w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))
      w(ixO^S,mom(3))=w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))
      ! 2nd order accuacy constant value extrapolation
      do ix1=ixOmin1,ixOmax1
        delt_x1 = x(ix1+0,ixOmin2,ixOmin3,1)-x(ix1-1,ixOmin2,ixOmin3,1)
        delt_x2 = x(ix1-1,ixOmin2,ixOmin3,1)-x(ix1-2,ixOmin2,ixOmin3,1)
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))=                   &
                (1.0d0/(delt_x2**2+2.0d0*delt_x1*delt_x2))*(-delt_x1**2* &
                    w(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))*     &
                   (x(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2*       &
                sin(x(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*         &
                    dthe*dphi) + (delt_x1+delt_x2)**2*                   &
                    w(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))*     &
                   (x(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2*       &
                sin(x(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*         &
                    dthe*dphi))/                                         &
                   (x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2*       &
                sin(x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*         &
                   dthe*dphi)

        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2))=                   &
                (1.0d0/(delt_x2**2+2.0d0*delt_x1*delt_x2))*(-delt_x1**2* &
                    w(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2))*     &
                   (x(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*          &
                sin(x(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*         &
                    drad(0)*dphi) + (delt_x1+delt_x2)**2*                &
                  w(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2))*       &
                   (x(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*          &
                sin(x(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*         &
                    drad(1)*dphi))/                                      &
                   (x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*          &
                sin(x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*         &
                    drad(2)*dphi)

        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3))=                   &
                (1.0d0/(delt_x2**2+2.0d0*delt_x1*delt_x2))*(-delt_x1**2* &
                    w(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3))*     &
                   (x(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*          &
                    drad(0)*dthe) + (delt_x1+delt_x2)**2*                &
                    w(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3))*     &
                   (x(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*          &
                    drad(1)*dthe))/                                      &
                   (x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*          &
                    drad(2)*dthe)
      enddo
    case(3)    ! theta mininum surface
      w(ixO^S,rho_)  = w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
      w(ixO^S,p_)=w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,p_)
      w(ixO^S,mom(1))=w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(1))
      w(ixO^S,mom(2))=w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(2))
      w(ixO^S,mom(3))=w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(3))
      ! 2nd order accuacy constant value extrapolation
      do ix2=ixOmax2,ixOmin2,-1
       w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1))=(1.0d0/3.0d0)*     &
                  (-w(ixOmin1:ixOmax1,ix2+2,ixOmin3:ixOmax3,mag(1))*    &
                   (x(ixOmin1:ixOmax1,ix2+2,ixOmin3:ixOmax3,1)**2*   &
                sin(x(ixOmin1:ixOmax1,ix2+2,ixOmin3:ixOmax3,2)))     &
             +4.0d0*w(ixOmin1:ixOmax1,ix2+1,ixOmin3:ixOmax3,mag(1))*    &
                   (x(ixOmin1:ixOmax1,ix2+1,ixOmin3:ixOmax3,1)**2*   &
                sin(x(ixOmin1:ixOmax1,ix2+1,ixOmin3:ixOmax3,2))))/   &
                   (x(ixOmin1:ixOmax1,ix2+0,ixOmin3:ixOmax3,1)**2*   &
                sin(x(ixOmin1:ixOmax1,ix2+0,ixOmin3:ixOmax3,2)))
       w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(2))=(1.0d0/3.0d0)*     &
                  (-w(ixOmin1:ixOmax1,ix2+2,ixOmin3:ixOmax3,mag(2))*    &
                   (x(ixOmin1:ixOmax1,ix2+2,ixOmin3:ixOmax3,1)*      &
                sin(x(ixOmin1:ixOmax1,ix2+2,ixOmin3:ixOmax3,2)))     &
             +4.0d0*w(ixOmin1:ixOmax1,ix2+1,ixOmin3:ixOmax3,mag(2))*    &
                   (x(ixOmin1:ixOmax1,ix2+1,ixOmin3:ixOmax3,1)*      &
                sin(x(ixOmin1:ixOmax1,ix2+1,ixOmin3:ixOmax3,2))))/   &
                   (x(ixOmin1:ixOmax1,ix2+0,ixOmin3:ixOmax3,1)*      &
                sin(x(ixOmin1:ixOmax1,ix2+0,ixOmin3:ixOmax3,2)))
       w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(3))=(1.0d0/3.0d0)*     &
                  (-w(ixOmin1:ixOmax1,ix2+2,ixOmin3:ixOmax3,mag(3))*    &
                   (x(ixOmin1:ixOmax1,ix2+2,ixOmin3:ixOmax3,1))      &
             +4.0d0*w(ixOmin1:ixOmax1,ix2+1,ixOmin3:ixOmax3,mag(3))*    &
                   (x(ixOmin1:ixOmax1,ix2+1,ixOmin3:ixOmax3,1)))/    &
                   (x(ixOmin1:ixOmax1,ix2+0,ixOmin3:ixOmax3,1))
     enddo
    case(4)    ! theta maxinum surface
      w(ixO^S,rho_)  = w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,rho_)
      w(ixO^S,p_)=     w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,p_)
      w(ixO^S,mom(1))=w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(1))
      w(ixO^S,mom(2))=w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(2))
      w(ixO^S,mom(3))=w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(3))
      ! 2nd order accuacy constant value extrapolation
      do ix2=ixOmin2,ixOmax2
       w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1))=(1.0d0/3.0d0)*     &
                  (-w(ixOmin1:ixOmax1,ix2-2,ixOmin3:ixOmax3,mag(1))*    &
                   (x(ixOmin1:ixOmax1,ix2-2,ixOmin3:ixOmax3,1)**2*   &
                sin(x(ixOmin1:ixOmax1,ix2-2,ixOmin3:ixOmax3,2)))     &
             +4.0d0*w(ixOmin1:ixOmax1,ix2-1,ixOmin3:ixOmax3,mag(1))*    &
                   (x(ixOmin1:ixOmax1,ix2-1,ixOmin3:ixOmax3,1)**2*   &
                sin(x(ixOmin1:ixOmax1,ix2-1,ixOmin3:ixOmax3,2))))/   &
                   (x(ixOmin1:ixOmax1,ix2+0,ixOmin3:ixOmax3,1)**2*   &
                sin(x(ixOmin1:ixOmax1,ix2+0,ixOmin3:ixOmax3,2)))
       w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(2))=(1.0d0/3.0d0)*     &
                  (-w(ixOmin1:ixOmax1,ix2-2,ixOmin3:ixOmax3,mag(2))*    &
                   (x(ixOmin1:ixOmax1,ix2-2,ixOmin3:ixOmax3,1)*      &
                sin(x(ixOmin1:ixOmax1,ix2-2,ixOmin3:ixOmax3,2)))     &
             +4.0d0*w(ixOmin1:ixOmax1,ix2-1,ixOmin3:ixOmax3,mag(2))*    &
                   (x(ixOmin1:ixOmax1,ix2-1,ixOmin3:ixOmax3,1)*      &
                sin(x(ixOmin1:ixOmax1,ix2-1,ixOmin3:ixOmax3,2))))/   &
                   (x(ixOmin1:ixOmax1,ix2+0,ixOmin3:ixOmax3,1)*      &
                sin(x(ixOmin1:ixOmax1,ix2+0,ixOmin3:ixOmax3,2)))
       w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(3))=(1.0d0/3.0d0)*     &
                  (-w(ixOmin1:ixOmax1,ix2-2,ixOmin3:ixOmax3,mag(3))*    &
                   (x(ixOmin1:ixOmax1,ix2-2,ixOmin3:ixOmax3,1))      &
             +4.0d0*w(ixOmin1:ixOmax1,ix2-1,ixOmin3:ixOmax3,mag(3))*    &
                   (x(ixOmin1:ixOmax1,ix2-1,ixOmin3:ixOmax3,1)))/    &
                   (x(ixOmin1:ixOmax1,ix2+0,ixOmin3:ixOmax3,1))
     enddo
    case(5)    ! phi mininum surface
      w(ixO^S,rho_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,rho_)
      w(ixO^S,p_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,p_)
      w(ixO^S,mom(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(1))
      w(ixO^S,mom(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(2)) 
      w(ixO^S,mom(3))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(3)) 
      ! 2nd order accuacy constant value extrapolation
      do ix3=ixOmax3,ixOmin3,-1
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(1):mag(3))=(1.0d0/3.0d0)* &
              (     -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+2,mag(1):mag(3)) &
               +4.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+1,mag(1):mag(3)))
      enddo
    case(6)    ! phi maxinum surface
      w(ixO^S,rho_)  = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,rho_)
      w(ixO^S,p_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,p_)
      w(ixO^S,mom(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(1))
      w(ixO^S,mom(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(2))
      w(ixO^S,mom(3))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(3))
      ! 2nd order accuacy constant value extrapolation
      do ix3=ixOmin3,ixOmax3
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(1):mag(3))=(1.0d0/3.0d0)* &
             (     -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-2,mag(1):mag(3)) &
              +4.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-1,mag(1):mag(3)))
      enddo
    case default
      call mpistop("Special boundary is not defined for this region")
    end select
  end subroutine specialbound_usr

  !==============================================================================
  !> Set the "eta" array for resistive MHD based on w or the
  !> "current" variable which has components between idirmin and 3.
  !==============================================================================
  subroutine specialres_usr(w,ixI^L,ixO^L,idirmin,x,current,eta)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L, idirmin
    double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision             :: current(ixI^S,7-2*ndir:3), eta(ixI^S)
    double precision             :: eta_0, eta_2, j_c, j_cr, j_mag(ixI^S), eta_max
    integer                      :: ix^D

    eta_0 = 0.0d0
    eta_2 = 5.0d-4       ! The normalization factor eta_norm = 8.099d17 cm^2 s^-1(1.0d-3)
    eta_max = 10.0*eta_2
    j_c = 1000         ! Equivalent to 1.6d-6 A/cm^2. The normalization factor J_0 = 2.28835d-11 A/cm^2. 6.99d4
    j_cr = 1.d0/j_c
    j_mag = dsqrt(sum((current(ixI^S,:))**2,dim=ndim+1))
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      if (j_mag(ix^D) < j_c) then
        eta(ix^D) = eta_0
      else
        eta(ix^D) = eta_2*((j_mag(ix^D)-j_c)*j_cr)**2 + eta_0
      endif
    {end do\}

    where(eta(ixO^S)>eta_max)
      eta(ixO^S)=eta_max
    endwhere
  end subroutine specialres_usr

  !==============================================================================
  ! Purpose: get gravity field
  !==============================================================================
  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=grav_sun*(radius_sun/x(ixO^S,1))**2
  end subroutine

  !==============================================================================
  ! Purpose: get gravity field
  !==============================================================================
  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)
    double precision                :: ggrid(ixI^S)

    gravity_field=0.d0
    call getggrav(ggrid,ixI^L,ixO^L,x)
    gravity_field(ixO^S,1)=ggrid(ixO^S)

  end subroutine gravity

  !==============================================================================
  ! Purpose: add density diffusion as a special source
  !==============================================================================
  subroutine specialsource_usr(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    !double precision :: density_deviation(ixI^S), density_diffusion(ixI^S), gradd(ixI^S), &
    !                    graddensity(ixI^S,ndim), density03(ixI^S)
    !integer :: ix^D, idim, ixA^L, ind_nx1
    integer :: ix^D, ind_nx1
    
    ix2 = 1
    ix3 = 1
    do ix1=ixOmin1,ixOmax1
      ind_nx1 = ceiling((x(ix^D,1) - xprobmin1)/dr)
      where(w(ix1,:,:,rho_) < small_density0 ) w(ix1,:,:,rho_) = small_density0 
    end do

    !ix2 = 1
    !ix3 = 1
    !do ix1=ixImin1,ixImax1
    !  ind_nx1 = ceiling((x(ix^D,1) - xprobmin1)/dr)
    !  density03(ix1,:,:) = density0(ind_nx1)
    !end do
      
    !density_deviation(ixI^S)=wCT(ixI^S,rho_)-density03
    !ixA^L=ixO^L^LADD1;
    !if (ixImin^D>ixAmin^D.or.ixImax^D<ixAmax^D|.or.) &
    !     call mpistop("Error in specialsource_usr: Non-conforming input limits")

    !do idim=1,ndim
    !  call gradient(density_deviation,ixI^L,ixA^L,idim,gradd)
    !  graddensity(ixA^S,idim)=gradd(ixA^S)
    !end do
    !call divvector(graddensity, ixI^L, ixO^L, density_diffusion)
    ! !density_diffusion(ixO^S)=density_diffusion(ixO^S)*mu_density &
    ! !                       /(^D&1.0d0/block%ds(ixO^S,^D)**2+)
    !density_diffusion(ixO^S)=density_diffusion(ixO^S)*mu_density*qdt
    !w(ixO^S,rho_)=w(ixO^S,rho_)+density_diffusion(ixO^S)

  end subroutine specialsource_usr

  !==============================================================================
  ! Purpose: Get new dt restricted by density diffusion in specialsource_usr
  !==============================================================================
  subroutine special_getdt(w,ixG^L,ix^L,dtnew,dx^D,x)

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
    double precision, intent(in) :: w(ixG^S,1:nw)
    double precision, intent(inout) :: dtnew
    integer :: idims

    dtnew=bigdouble
    !do idims=1,ndim
    !   dtnew=min(dtnew,0.4d0*minval(block%dx(ix^S,idims)/vinf))
    !enddo
    dtnew=min(dtnew,minval(1.0d0/(^D&1.0d0/block%ds(ix^S,^D)**2+)/mu_density))

  end subroutine special_getdt

  !==============================================================================
  ! Purpose: Enforce additional refinement or coarsening. One can use the
  !          coordinate info in x and/or time qt=t_n and w(t_n) values w.
  !==============================================================================
 subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
   use mod_global_parameters

   integer, intent(in) :: igrid, level, ixI^L, ixO^L
   double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
   integer, intent(inout) :: refine, coarsen
   double precision :: rho(ixI^S), current(ixI^S,1:ndir), btotal(ixI^S,1:ndir)
   double precision :: b2(ixI^S), absj(ixO^S), absb(ixO^S)
   double precision :: local_refine_threshold(ixO^S), coarsen_threshold(ixO^S)
   double precision :: one_over_dx
   integer :: idirmin, idim, idir

   ! fix the bottom layer to the highest level
   if (any(x(ixO^S,1)<=xprobmin1+0.05d0) ) then
       refine=1
       coarsen=-1
   else 
       btotal = w(ixI^S,mag(:))
       rho = w(ixI^S,rho_)
       call curlvector(btotal,ixI^L,ixO^L,current,idirmin,1,ndir)
       b2 = sum(btotal(ixI^S,:)**2,dim=ndim+1)

       absj(ixO^S) = dsqrt(sum(current(ixO^S,:)**2,dim=ndim+1))
       absb(ixO^S) = dsqrt(b2(ixO^S))

       one_over_dx  = 1/min({rnode(rpdx^D_,igrid)},bigdouble)
       
       ! 50/25 can be set as a parameter in amrvac.par for testing
       local_refine_threshold = absj/absb-one_over_dx*50 ! 50 is a adjustable parameter
       coarsen_threshold = absj/absb-one_over_dx*25 ! 25 is a adjustable parameter

       if (any(local_refine_threshold>0) .and. level<refine_max_level) then
           refine = 1
           coarsen = -1
       else if (any(coarsen_threshold<0)) then
           refine = -1
           coarsen = 1
       else
           refine = 0
           coarsen = 0
       end if
   end if
  end subroutine special_refine_grid

  !==============================================================================
  ! Purpose: 
  !   this subroutine can be used in convert, to add auxiliary variables to the
  !   converted output file, for further analysis using tecplot, paraview, ....
  !   these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  !   the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  !   corresponding normalization values (default value 1)
  !==============================================================================
  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision                   :: tmp(ixI^S),dip(ixI^S),divb(ixI^S),B2(ixI^S)
    double precision, dimension(ixI^S,1:ndir) :: Btotal,qvec,curlvec
    integer                            :: ix^D,idirmin,idims,idir,jdir,kdir

    ! Btotal & B^2
    Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1):mag(ndir))
    B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)
    ! output Alfven wave speed B/sqrt(rho)
  !  w(ixO^S,nw+1)=dsqrt(B2(ixO^S)/w(ixO^S,rho_))
    ! output divB1
  !  call divvector(Btotal,ixI^L,ixO^L,divb)
  !  w(ixO^S,nw+2)=0.5d0*divb(ixO^S)/dsqrt(B2(ixO^S))/(^D&1.0d0/dxlevel(^D)+)
    ! output the plasma beta p*2/B**2
  !  w(ixO^S,nw+3)=w(ixO^S,rho_)*T_corona*two/B2(ixO^S)
    ! store current
    call curlvector(Btotal,ixI^L,ixO^L,curlvec,idirmin,1,ndir)
    do idir=1,ndir
      w(ixO^S,nw+idir)=curlvec(ixO^S,idir)
    end do
    w(ixO^S,nw+4)=w(ixO^S,mag(1))
    
  !  w(ixO^S,nw+8)=w(ixO^S,mag(2))
! w(ixO^S,nw+9)=w(ixO^S,mag(3))
  end subroutine specialvar_output

  !==============================================================================
  ! Purpose: names for special variable output
  !==============================================================================
  subroutine specialvarnames_output(varnames)
    use mod_global_parameters
    character(len=*) :: varnames

 !   varnames='Alfv divB beta j1 j2 j3 br bth bph'
    varnames='j1 j2 j3 br'
  end subroutine specialvarnames_output

end module mod_usr

