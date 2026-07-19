module mod_usr   ! computing potential field source surface (PFSS) model and RBSL model
  use mod_mhd
  use mod_pfss
  implicit none
  ! some global parameters
  double precision :: k_B,miu0,mass_H,usr_grav,SRadius,rhob,Tiso
  logical, save :: firstusrglobaldata=.true.

  ! global parameters for RBSL
  integer :: np,nnp
  double precision, allocatable :: x_axis(:,:)
  double precision :: F_flx, a0
  logical          :: I_Helix
  ! global parameters for RBSL boundary
  integer          :: ixgI^L,ixgO^L 
  integer          :: m,n,p
  double precision, allocatable :: br_rbsl(:,:)

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
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output

    call set_coordinate_system("spherical")
    call mhd_activate()
  end subroutine usr_init

  !==============================================================================
  ! Purpose: to initialize user public parameters and reset global parameters.
  !          Input data are also read here.
  !==============================================================================
  subroutine initglobaldata_usr()
    use mod_global_parameters
    character(len=70) :: file_path
    logical :: flip
    double precision, allocatable :: b_r0(:,:),theta(:),phi(:)
    integer :: xm,ym,amode,file_handle
    integer :: lat0,lat1,lon0,lon1
    integer :: iphi, itheta
    integer, dimension(MPI_STATUS_SIZE) :: statuss
    character(len=70) :: mapname

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

    usr_grav=-2.74d4*unit_length/unit_velocity**2  ! solar gravity
    SRadius=6.955d10/unit_length                   ! Solar radius
    rhob=1.0d0                                     ! bottom density
    Tiso=1.0d6/unit_temperature                    ! uniform temperature
    
    np = 6000
    nnp = np*2
    allocate(x_axis(nnp,ndim))

    F_flx = 10.0d20/(unit_magneticfield*unit_length**2)  !---------------------------------------
    a0   = 4.0d9/unit_length
    I_Helix = .false.
    file_path = './rbsl/path.dat'
    flip = .true.
    if(mype==0 .and. .not.convert) then
      print*,' '
      print*,'===== RBSL parameters ====='
      print*,'Number of axis points:',np
      print*,'Flux (Mx):',F_flx*(unit_magneticfield*unit_length**2)
      print*,'Minor radius (Mm):',a0*unit_length/1.0d8              !1.2Mm
      print*,'Positive helicity:',I_Helix
      print*,'Axis file path:   ',file_path
      print*,'Flip axis curve:',flip
      print*,' '
    endif
    call flux_rope_insertion(file_path,flip,x_axis)
 
    if(mype==0 .and. .not.convert) then 
      ! Colat-lon system (tp CE) : m * n
      m=domain_nx2*2**(refine_max_level-1)           !Colat-lon system (Dim 1)
      n=domain_nx3*2**(refine_max_level-1)           !Colat-lon system (Dim 1) 
      p=domain_nx1*2**(refine_max_level-1)           ! n cells from solar surface to source surface
      ixgOmin1=1;
      ixgOmin2=1;
      ixgOmin3=1;
      ixgOmax1=1;
      ixgOmax2=m;
      ixgOmax3=n;
      ixgI^L=ixgO^L^LADD1;
      allocate(br_rbsl(m,n))
      call prep_rbsl(ixgI^L,ixgO^L,br_rbsl)

      mapname = './potential_boundary/pfss_boundary.dat'
      call MPI_FILE_OPEN(MPI_COMM_SELF,mapname,MPI_MODE_RDONLY,MPI_INFO_NULL,&
                         file_handle,ierrmpi)
      call MPI_FILE_READ(file_handle,xm,1,MPI_INTEGER,statuss,ierrmpi)
      call MPI_FILE_READ(file_handle,ym,1,MPI_INTEGER,statuss,ierrmpi)
      allocate(b_r0(xm,ym))
      allocate(theta(ym))
      allocate(phi(xm))
      call MPI_FILE_READ(file_handle,theta,ym,MPI_DOUBLE_PRECISION,&
                         statuss,ierrmpi)
      call MPI_FILE_READ(file_handle,phi,xm,MPI_DOUBLE_PRECISION,&
                         statuss,ierrmpi)
      call MPI_FILE_READ(file_handle,b_r0,xm*ym,MPI_DOUBLE_PRECISION,&
                         statuss,ierrmpi)
      call MPI_FILE_CLOSE(file_handle,ierrmpi)
      print*,' '
      print*,'Radial magnetic field parameters before RBSL field subtraction:'
      print*,'nphi,ntheta',xm,ym
      print*,'theta range:',minval(theta),maxval(theta)
      print*,'phi range:',minval(phi),maxval(phi)
      print*,'Brmax,Brmin',maxval(b_r0),minval(b_r0)
      print*,' '
      print*,'Radial magnetic field parameters after RBSL field subtraction:'

      ! Subregion positions provided by '20120123_0300_hmi/02preprocess/02preprocess_largeFOV/02prep_boundary.pro'
      lat0 = 580
      lat1 = 819
      lon0 = 1930
      lon1 = 2409

      do itheta=1,domain_nx2 
      do iphi=1,domain_nx3
        b_r0(iphi+lon0,itheta-1+ym-lat1) = b_r0(iphi+lon0,itheta-1+ym-lat1) - br_rbsl(domain_nx2+1-itheta,iphi)*unit_magneticfield
      end do
      end do

      mapname = './potential_boundary/pfss_boundary1.dat'
      amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
      call MPI_FILE_OPEN(MPI_COMM_SELF,mapname,amode, &
                           MPI_INFO_NULL,file_handle,ierrmpi)
      call MPI_FILE_WRITE(file_handle,xm,1,MPI_INTEGER,statuss,ierrmpi)
      call MPI_FILE_WRITE(file_handle,ym,1,MPI_INTEGER,statuss,ierrmpi)
      call MPI_FILE_WRITE(file_handle,theta,ym,MPI_DOUBLE_PRECISION,&
                         statuss,ierrmpi)
      call MPI_FILE_WRITE(file_handle,phi,xm,MPI_DOUBLE_PRECISION,&
                         statuss,ierrmpi)
      call MPI_FILE_WRITE(file_handle,b_r0,xm*ym,MPI_DOUBLE_PRECISION,&
                         statuss,ierrmpi)
      call MPI_FILE_CLOSE(file_handle,ierrmpi)
    endif

    call MPI_BARRIER(icomm,ierrmpi)
 
    ! prepare magnetogram at bottom
    trunc=.true.
    lmax=360   
    call harm_coef('potential_boundary/pfss_boundary1.dat')
    firstusrglobaldata=.false.
 !   print*,'lmax?:',lmax
  end subroutine initglobaldata_usr

  !==============================================================================
  ! Purpose: to initialize the initial condition
  !==============================================================================
  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: bpf(ixI^S,1:ndir)
    logical, save :: first=.true.
    integer :: ix^D, idirmin
    double precision :: xglC(ixI^S,1:ndir)
    double precision :: AfrC(ixI^S,1:ndir), AfrS(ixI^S,1:ndir)
    double precision :: BfrS(ixI^S,1:ndir)

    ! +++++++++ Special Regular :: 1=2, 2=3, 3=1 +++++++++
    xglC(ixI^S,1) = x(ixI^S,1) * cos(0.50d0*dpi-x(ixI^S,2)) * cos(x(ixI^S,3))
    xglC(ixI^S,2) = x(ixI^S,1) * cos(0.50d0*dpi-x(ixI^S,2)) * sin(x(ixI^S,3))
    xglC(ixI^S,3) = x(ixI^S,1) * sin(0.50d0*dpi-x(ixI^S,2))

    if(mype==0 .and. first) then
      write(*,*)' '
      write(*,*)'initializing grids with the PFSS and RBSL model...'
      first=.false.
    endif

    ! use RBSL model to compute B  
    ! Note that we have to compute AfrC within the range of ixI^L, but not ixO^L,
    ! since we need one grid more to compute br_mfr as follows.
    call RBSL_flux_rope(ixI^L,ixI^L,nnp,a0,F_flx,I_Helix,xglC,x_axis,AfrC)
    call Cart2SphereVector(ixI^L,ixI^L,x,AfrC,AfrS)
    call curlvector(AfrS,ixI^L,ixO^L,BfrS,idirmin,1,ndir)
    
    ! use PFSS model to compute B  
    call pfss(ixI^L,ixO^L,bpf,x)

    !w(ixO^S,mag(1):mag(ndir)) = zero, hereeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
    w(ixO^S,mag(1):mag(ndir)) = bpf(ixO^S,1:ndir)/unit_magneticfield + BfrS(ixO^S,1:ndir)
    !print*,maxval(BfrS(ixO^S,1:ndir))
    !w(ixO^S,mag(1):mag(ndir)) = bpf(ixO^S,1:ndir)/unit_magneticfield
    w(ixO^S,mom(:))=zero
    w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                  (1.d0/SRadius-1.d0/(x(ixO^S,1)+SRadius)))
                  
  end subroutine initonegrid_usr

  !==============================================================================
  ! Purpose: read path.dat, convert the axis to Cartesian coordinates
  !==============================================================================
  subroutine flux_rope_insertion(file_path,flip,x_axis)

  character(len=*), intent(in) :: file_path
  logical, intent(in) :: flip
  double precision, intent(out) :: x_axis(nnp,ndim)

  integer          :: iax, ixp
  double precision :: lonax(np),latax(np),rax(np)
  double precision ::  lonm(np), latm(np), rm(np)

  double precision :: lon0,lat0,lon1,lat1
  double precision :: lonc,latc,rc
  double precision :: x0,y0,z0
  double precision :: x1,y1,z1
  double precision :: xc,yc,zc
  double precision ::  v1(nnp),  v2(nnp),  v3(nnp)
  double precision ::  xs(nnp),  ys(nnp),  zs(nnp)

    open(unit=62,file=file_path,action='read',status='old',access='stream')
    do iax=1,np
      read(62) lonax(iax),latax(iax),rax(iax)
    !read(62,'(3f11.5)') lonax(iax),latax(iax),rax(iax)
    end do
    close(62)

  if(mype==0) then
    write(6,'(a)') ' ===== The spherical coordinates of the filament from SCC Measure ====='
    write(6,'(a)') '          ===== Only the two footprints are printed here ====='
    write(6,'(a)') '            Longitude(deg) , Latitude(deg) , Radius( x Rsun)            '
    !do iax=1,np
    !  write(6,'(a,3f16.5,a)') '       ',lonax(iax),latax(iax),rax(iax),'          '
    !end do
    write(6,'(a,3f16.5,a)') '       ',lonax(1),latax(1),rax(1),'          '
    write(6,'(a,3f16.5,a)') '       ',lonax(np),latax(np),rax(np),'          '
    write(6,'(a)') ' '
  end if

    lonax=lonax*dpi/180.0d0
    latax=latax*dpi/180.0d0

    lon0=lonax(1)
    lat0=latax(1)
    lon1=lonax(np)
    lat1=latax(np)

    x0=cos(lat0)*cos(lon0)
    y0=cos(lat0)*sin(lon0)
    z0=sin(lat0)

    x1=cos(lat1)*cos(lon1)
    y1=cos(lat1)*sin(lon1)
    z1=sin(lat1)

    xc=0.50d0*(x0+x1)
    yc=0.50d0*(y0+y1)
    zc=0.50d0*(z0+z1)

    call get_lonlat(xc,yc,zc,lonc,latc)
    ! Note that get_lonlat returs longitude in the range of [-0.5*Pi, 0.5*Pi],
    ! While in this case the longitude is measured from the backside of the Sun.
    lonc = lonc + dpi
    rc=dsqrt(xc**2+yc**2+zc**2)

    call get_mirror_path(lonc,latc,rc,np,&
                         lonax,latax,rax,lonm,latm,rm)
    do ixp=1,np
      v1(ixp)   =   rax(ixp)
      v2(ixp)   = latax(ixp)
      v3(ixp)   = lonax(ixp)
      v1(np+ixp)=    rm(np-ixp+1)
      v2(np+ixp)=  latm(np-ixp+1)
      v3(np+ixp)=  lonm(np-ixp+1)
    end do

    if (flip) then
      ! to flip the axis curve, thus the path starts from the positive polarity
      do ixp=1,np
        xs(ixp)    = v1(np-ixp+1)   * cos(v2(np-ixp+1))   * cos(v3(np-ixp+1))
        ys(ixp)    = v1(np-ixp+1)   * cos(v2(np-ixp+1))   * sin(v3(np-ixp+1))
        zs(ixp)    = v1(np-ixp+1)   * sin(v2(np-ixp+1))
        xs(np+ixp) = v1(2*np-ixp+1) * cos(v2(2*np-ixp+1)) * cos(v3(2*np-ixp+1))
        ys(np+ixp) = v1(2*np-ixp+1) * cos(v2(2*np-ixp+1)) * sin(v3(2*np-ixp+1))
        zs(np+ixp) = v1(2*np-ixp+1) * sin(v2(2*np-ixp+1))
      end do
    else
      do ixp=1,np
        xs(ixp)    = v1(ixp)   * cos(v2(ixp))   * cos(v3(ixp))
        ys(ixp)    = v1(ixp)   * cos(v2(ixp))   * sin(v3(ixp))
        zs(ixp)    = v1(ixp)   * sin(v2(ixp))
        xs(np+ixp) = v1(np+ixp) * cos(v2(np+ixp)) * cos(v3(np+ixp))
        ys(np+ixp) = v1(np+ixp) * cos(v2(np+ixp)) * sin(v3(np+ixp))
        zs(np+ixp) = v1(np+ixp) * sin(v2(np+ixp))
      end do
    endif

    x_axis(:,1) = xs
    x_axis(:,2) = ys
    x_axis(:,3) = zs

  end subroutine flux_rope_insertion

  !==============================================================================
  ! Purpose: convert Caretesian coordinates to spherical ones
  !==============================================================================
  subroutine get_lonlat(x_in,y_in,z_in,lon,lat)

  real*8 :: x_in,y_in,z_in
  real*8 :: x,y,z
  real*8 :: lon,lat
  real*8 :: Ebm

    Ebm=dsqrt(x_in**2+y_in**2+z_in**2)
    x=x_in/Ebm
    y=y_in/Ebm
    z=z_in/Ebm

    if(x==0.and.y==0) then
      if(z.ge.0) then
        lat=0.5d0*dpi
      else
        lat=-0.5d0*dpi
      end if
        lon=0.0d0
    else if(z==0) then
           lat=0.0d0
           if(y==0) then
             lon=0.0d0
           else if(x==0) then
             if(y.ge.0) then
               lon=0.5d0*dpi
             else
               lon=-0.5d0*dpi
             end if
           else
             lon=atan(y/x)
           end if
         else
           lat=asin(z)
           if(y==0) then
             lon=0.0d0
           else if(x==0) then
             if(y.ge.0) then
               lon=0.5d0*dpi
             else
               lon=-0.5d0*dpi
             end if
           else
             lon=atan(y/x)
           end if
    end if

    return

  end subroutine get_lonlat

  !==============================================================================
  ! Purpose: get mirror path. And the mirror plane is perpendicular to the radius
  !          vector of the middle point of the two footprints, and crossing it.
  !==============================================================================
  subroutine get_mirror_path(lonc,latc,rc,np,&
                             lon_in,lat_in,r_in,lon_out,lat_out,r_out)

  integer,intent(in)            :: np
  double precision, intent(in)  :: lonc,latc,rc
  double precision, intent(in)  :: lon_in(np),lat_in(np),r_in(np)
  double precision, intent(out) :: lon_out(np),lat_out(np),r_out(np)

  integer          :: inp
  double precision :: lonp,latp
  double precision :: xp,yp,zp
  double precision :: xtmp,ytmp,ztmp
  double precision :: lontmp,lattmp
  double precision :: lona(np),lata(np),ra(np)
  double precision :: lonm(np),latm(np),rm(np)
  double precision :: xa(np),ya(np),za(np)
  double precision :: xm(np),ym(np),zm(np)

    lona=lon_in-lonc
    lata=lat_in-latc
    ra=r_in
    xa=ra*cos(lata)*cos(lona)
    ya=ra*cos(lata)*sin(lona)
    za=ra*sin(lata)

    lonp=lonc-lonc
    latp=latc-latc
    xp=rc*cos(latp)*cos(lonp)

    xm=2*xp-xa
    ym=ya
    zm=za

    do inp=1,np
      xtmp=xm(inp)
      ytmp=ym(inp)
      ztmp=zm(inp)
      lontmp=0.0d0
      lattmp=0.0d0
      call get_lonlat(xtmp,ytmp,ztmp,lontmp,lattmp)
      lonm(inp)=lontmp
      latm(inp)=lattmp
    end do

    lon_out=lonm+lonc
    lat_out=latm+latc
    rm=dsqrt(xm**2+ym**2+zm**2)
    r_out=rm

    return

  end subroutine get_mirror_path

  !==============================================================================
  ! Purpose: calculate RBSL magnetic field on the bottom boundary
  !==============================================================================
  subroutine prep_rbsl(ixI^L,ixO^L,br_mfr)

    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(out) :: br_mfr(m,n)

    integer :: ix^D
    double precision :: dx01,dx02,dx03
    double precision :: xgl(ixI^S,1:ndim), xglC(ixI^S,1:ndim)
    double precision :: AfrC(ixI^S,1:ndim), AfrS(ixI^S,1:ndim)

    dx01=(xprobmax1-xprobmin1)/p;
    dx02=(xprobmax2-xprobmin2)/m;
    dx03=(xprobmax3-xprobmin3)/n;
    !write(*,*) 'xprobmin1,xprobmax1:',xprobmin1,xprobmax1
    !write(*,*) 'xprobmin2,xprobmax2:',xprobmin2,xprobmax2
    !write(*,*) 'xprobmin3,xprobmax3:',xprobmin3,xprobmax3
    {do ix^DB=ixgImin^DB,ixgImax^DB\}
      xgl(ix^D,1)=(ix1-1)*dx01+xprobmin1
      xgl(ix^D,2)=(ix2-1)*dx02+xprobmin2
      xgl(ix^D,3)=(ix3-1)*dx03+xprobmin3
    {end do\}
    ! +++++++++ Special Regular :: 1=2, 2=3, 3=1 +++++++++
    xglC(ixI^S,1) = xgl(ixI^S,1) * cos(0.50d0*dpi-xgl(ixI^S,2)) * cos(xgl(ixI^S,3))
    xglC(ixI^S,2) = xgl(ixI^S,1) * cos(0.50d0*dpi-xgl(ixI^S,2)) * sin(xgl(ixI^S,3))
    xglC(ixI^S,3) = xgl(ixI^S,1) * sin(0.50d0*dpi-xgl(ixI^S,2))

    print*,'Call prep_rbsl: Calculating Br_rbsl on the bottom boundary ...'   
    ! Note that we cannot use RBSL_flux_rope to calculate BfrC directly, since the 
    ! magnetic field is located on the bottom boundary and we compute the whole 
    ! boundary, but do not divide it into blocks. So, the curlvector() cannot be 
    ! applied here.
    ! And note that we have to compute AfrC within the range of ixI^L, but not ixO^L,
    ! since we need one grid more to compute br_mfr as follows.
    call RBSL_flux_rope(ixI^L,ixI^L,nnp,a0,F_flx,I_Helix,xglC,x_axis,AfrC)
    call Cart2SphereVector(ixI^L,ixI^L,xgl,AfrC,AfrS)

    ix1=1
    do ix3=ixOmin3,ixOmax3
    do ix2=ixOmin2,ixOmax2
      br_mfr(ix2,ix3) = 1.d0/(xgl(ix^D,1)*sin(xgl(ix^D,2)))*((sin(xgl(ix1,ix2+1,ix3,2)) &
                     *AfrS(ix1,ix2+1,ix3,3)-sin(xgl(ix1,ix2-1,ix3,2))*AfrS(ix1,ix2-1,ix3,3)) &
                     /(2.0d0*dx02) - (AfrS(ix1,ix2,ix3+1,2)-AfrS(ix1,ix2,ix3-1,2)) &
                     /(2.0d0*dx03))
    end do
    end do

  end subroutine prep_rbsl

  !==============================================================================
  ! Purpose: calculate RBSL magnetic field
  !==============================================================================
  subroutine RBSL_flux_rope(ixI^L,ixO^L,np,a,F_flx,positive_helicity,x,x_axis,Atotal,Bfr)

    integer, intent(in)             :: ixI^L,ixO^L,np
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    ! coordinates of flux rope axis (integral path)
    double precision, intent(in)    :: x_axis(np,1:ndim)
    ! cross-sectional radius of flux rope
    double precision, intent(in)    :: a
    ! net magnetic flux along flux rope axis
    double precision, intent(in)    :: F_flx
    ! is positive helicity
    logical, intent(in) :: positive_helicity
    ! vector potential
    double precision, intent(out)   :: Atotal(ixI^S,1:ndim)
    ! magnetic field
    double precision, optional, intent(out)   :: Bfr(ixI^S,1:ndim)

    ! net current along flux rope axis for azimuthal magnetic field
    double precision   :: I_cur
    ! vector potential for azimuthal magnetic field
    double precision   :: AIx(ixI^S,1:ndim)
    ! vector potential for axial magnetic field
    double precision   :: AFx(ixI^S,1:ndim)
    integer :: ix^D, ixp, idirmin,ixM^L
    double precision :: r_mag, KIr, KFr, dl, re_pi, sqrt1r, f52r,fsqrt6
    double precision :: Rpl(1:ndim), r_vec(1:ndim), Rcr(1:ndim)

    if(positive_helicity) then
      I_cur = 5.d0*sqrt(2.d0)*F_flx/(3.d0*4.d0*dpi*a)
    else
      I_cur =-5.d0*sqrt(2.d0)*F_flx/(3.d0*4.d0*dpi*a)
    end if

    re_pi=1.d0/dpi
    AIx = 0.d0
    AFx = 0.d0
    fsqrt6=1.d0/sqrt(6.d0)
    if(present(Bfr)) then
      ixM^L=ixI^L;
    else
      ixM^L=ixO^L;
    end if
    {do ix^DB=ixMmin^DB,ixMmax^DB\}
      do ixp=1,np
        ! position vector from source point to field point
        r_vec(:) = (x(ix^D,:) - x_axis(ixp,:))/a
        r_mag = sqrt(sum(r_vec(:)**2))
        ! calculate tangential vector of the axis which is a circle
        if (ixp == 1) then
          Rpl(:) = 0.5d0*(x_axis(ixp+1,:)-x_axis(np,:))
        else if (ixp==np) then
          Rpl(:) = 0.5d0*(x_axis(1,:)-x_axis(ixp-1,:))
        else
          Rpl(:) = 0.5d0*(x_axis(ixp+1,:)-x_axis(ixp-1,:))
        end if
        ! Rpl X r_vec
        Rcr(1) = Rpl(2)*r_vec(3) - Rpl(3)*r_vec(2)
        Rcr(2) = Rpl(3)*r_vec(1) - Rpl(1)*r_vec(3)
        Rcr(3) = Rpl(1)*r_vec(2) - Rpl(2)*r_vec(1)
        if (r_mag < 1.d0) then
          sqrt1r=sqrt(1.d0-r_mag**2)
          f52r=5.d0-2.d0*r_mag**2
          KIr = 2.d0*re_pi*(asin(r_mag)/r_mag + f52r*third*sqrt1r)
          KFr = 2.d0*re_pi/r_mag**2*(asin(r_mag)/r_mag-sqrt1r) + &
                2.d0*re_pi*sqrt1r + f52r*0.5d0*fsqrt6*(1.d0 - &
                2.d0*re_pi*asin((1.d0+2.d0*r_mag**2)/f52r))
        else
          KIr = 1.d0/r_mag
          KFr = KIr**3
        endif
        AIx(ix^D,:) = AIx(ix^D,:) + KIr*Rpl(:)
        AFx(ix^D,:) = AFx(ix^D,:) + KFr*Rcr(:)
      end do
      AIx(ix^D,:) = AIx(ix^D,:)*I_cur/a
      AFx(ix^D,:) = AFx(ix^D,:)*F_flx*0.25d0*re_pi/a**2
    {end do\}
    Atotal=AIx+AFx
    if(present(Bfr)) call curlvector(Atotal,ixI^L,ixO^L,Bfr,idirmin,1,ndir)

  end subroutine RBSL_flux_rope

  !==============================================================================
  ! Purpose: convert vectors in Cartesian coordinates to spherical ones
  !==============================================================================
  subroutine Cart2SphereVector(ixI^L,ixO^L,x,A_in,A_out)

    integer,intent(in)           :: ixI^L,ixO^L
    double precision,intent(in)  :: x(ixI^S,1:ndim)
    double precision,intent(in)  :: A_in(ixI^S,1:ndim)
    double precision,intent(out) :: A_out(ixI^S,1:ndim)

    double precision :: raddeg = dpi/ 180.
    double precision :: lon(ixI^S),lat(ixI^S)
    double precision :: bxCart(ixI^S),byCart(ixI^S),bzCart(ixI^S)
    double precision :: br(ixI^S),bth(ixI^S),bph(ixI^S)
    double precision :: a11(ixI^S),a12(ixI^S),a13(ixI^S)
    double precision :: a21(ixI^S),a22(ixI^S),a23(ixI^S)
    double precision :: a31(ixI^S),a32(ixI^S),a33(ixI^S)
    double precision :: latc(ixI^S),lonc(ixI^S),pAng(ixI^S)

    bxCart=zero
    byCart=zero
    bzCart=zero
       lon=zero
       lat=zero
       bph=zero
       bth=zero
        br=zero
     A_out=zero

    bxCart(ixO^S) = A_in(ixO^S,2)
    byCart(ixO^S) = A_in(ixO^S,3)
    bzCart(ixO^S) = A_in(ixO^S,1)
    lon(ixO^S) =  x(ixO^S,3)
    lat(ixO^S) =  0.5d0*dpi - x(ixO^S,2)
    latc = 0.0d0
    lonc = 0.0d0
    pAng = 0.0d0

    a11 = -sin(latc) * sin(pAng) * sin(lon - lonc) + cos(pAng) * cos(lon - lonc)
    a12 =  sin(latc) * cos(pAng) * sin(lon - lonc) + sin(pAng) * cos(lon - lonc)
    a13 = -cos(latc) * sin(lon - lonc)
    a21 = -sin(lat) * (sin(latc) * sin(pAng) * cos(lon - lonc) + cos(pAng) * sin(lon - lonc)) - cos(lat) * cos(latc) * sin(pAng)
    a22 =  sin(lat) * (sin(latc) * cos(pAng) * cos(lon - lonc) - sin(pAng) * sin(lon - lonc)) + cos(lat) * cos(latc) * cos(pAng)
    a23 = -cos(latc) * sin(lat) * cos(lon - lonc) + sin(latc) * cos(lat)
    a31 =  cos(lat) * (sin(latc) * sin(pAng) * cos(lon - lonc) + cos(pAng) * sin(lon - lonc)) - sin(lat) * cos(latc) * sin(pAng)
    a32 = -cos(lat) * (sin(latc) * cos(pAng) * cos(lon - lonc) - sin(pAng) * sin(lon - lonc)) + sin(lat) * cos(latc) * cos(pAng)
    a33 =  cos(lat) * cos(latc) * cos(lon - lonc) + sin(lat) * sin(latc)

    bph(ixO^S) = a11(ixO^S) * bxCart(ixO^S) +& 
                 a12(ixO^S) * byCart(ixO^S) +& 
                 a13(ixO^S) * bzCart(ixO^S)
    bth(ixO^S) = a21(ixO^S) * bxCart(ixO^S) +&
                 a22(ixO^S) * byCart(ixO^S) +&
                 a23(ixO^S) * bzCart(ixO^S)
     br(ixO^S) = a31(ixO^S) * bxCart(ixO^S) +&
                 a32(ixO^S) * byCart(ixO^S) +&
                 a33(ixO^S) * bzCart(ixO^S)

    A_out(ixO^S,1) =         br(ixO^S)
    A_out(ixO^S,2) = -1.0d0*bth(ixO^S)
    A_out(ixO^S,3) =        bph(ixO^S)

  end subroutine Cart2SphereVector

  !==============================================================================
  ! Purpose: to provide special boundary conditions set by users.
  !==============================================================================
  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters
    use mod_physics
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    !double precision :: ft,tfstop,tramp1,tramp2,coeffrho,vlimit,vsign
    !double precision :: delxdely,delxdelz,delydelx,delydelz,delzdelx,delzdely
    !double precision :: xlen^D,dxa^D,startpos^D
    !integer :: ix^D,ixIM^L,ixbc^D,af

    select case(iB)
    case(1)

    case(2)

    case(3)

    case(4)

    case(5)

    case(6)

    case default
      call mpistop("Special boundary is not defined for this region")
    end select
  end subroutine specialbound_usr

  !==============================================================================
  ! Purpose: get gravity field
  !==============================================================================
  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=usr_grav*(SRadius/(SRadius+x(ixO^S,3)))**2
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
    gravity_field(ixO^S,3)=ggrid(ixO^S)
  end subroutine gravity

  !==============================================================================
  ! Purpose: Enforce additional refinement or coarsening. One can use the
  !          coordinate info in x and/or time qt=t_n and w(t_n) values w.
  !==============================================================================
  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! fix the bottom layer to the highest level
    if (any(x(ixO^S,1)<=xprobmin1+0.075d0)) then
      refine=1
      coarsen=-1
    endif
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
    !Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1):mag(ndir))
    !B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)
      ! output Alfven wave speed B/sqrt(rho)
    !w(ixO^S,nw+1)=dsqrt(B2(ixO^S)/w(ixO^S,rho_))
      ! output divB1
    !call divvector(Btotal,ixI^L,ixO^L,divb)
    !w(ixO^S,nw+2)=0.5d0*divb(ixO^S)/dsqrt(B2(ixO^S))/(^D&1.0d0/dxlevel(^D)+)
      ! output the plasma beta p*2/B**2
    !w(ixO^S,nw+3)=w(ixO^S,rho_)*Tiso*two/B2(ixO^S)
      ! store current
    !call curlvector(Btotal,ixI^L,ixO^L,curlvec,idirmin,1,ndir)
    !do idir=1,ndir
    !  w(ixO^S,nw+3+idir)=curlvec(ixO^S,idir)
    !end do
    !w(ixO^S,nw+7)=w(ixO^S,mag(1))
    !w(ixO^S,nw+8)=w(ixO^S,mag(2))
    !w(ixO^S,nw+9)=w(ixO^S,mag(3))
    w(ixO^S,nw+1)=w(ixO^S,mag(1))
    w(ixO^S,nw+2)=w(ixO^S,mag(2))
    w(ixO^S,nw+3)=w(ixO^S,mag(3))
    ! find magnetic dips
!    dip=0.d0
!    do idir=1,ndir
!      call gradient(w(ixI^S,mag(3)),ixI^L,ixO^L,idir,tmp)
!      dip(ixO^S)=dip(ixO^S)+w(ixO^S,b0_+idir)*tmp(ixO^S)
!    end do
!    where(dabs(w(ixO^S,mag(3)))<0.08d0 .and. dip(ixO^S)>=0.d0)
!      w(ixO^S,nw+8)=1.d0
!    elsewhere
!      w(ixO^S,nw+8)=0.d0
!    end where
  end subroutine specialvar_output

  !==============================================================================
  ! Purpose: names for special variable output
  !==============================================================================
  subroutine specialvarnames_output(varnames)
    use mod_global_parameters
    character(len=*) :: varnames

    !varnames='Alfv divB beta j1 j2 j3 br bth bph'
    varnames='br bth bph'
  end subroutine specialvarnames_output

end module mod_usr





