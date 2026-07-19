module mod_usr   ! computing nonlinear force-free field in spherical geometry
  use mod_mhd
  implicit none
  ! variables for boundary condition
  double precision, allocatable, save :: Bnl(:,:,:)
  ! some global parameters
  double precision :: k_B,miu0,mass_H,usr_grav,SRadius,rhob,Tiso
  logical, save :: firstusrglobaldata=.true.

contains

  !==============================================================================
  ! Purpose: to include global parameters, set user methods, set coordinate 
  !          system and activate physics module.
  !==============================================================================
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    usr_set_parameters  => initglobaldata_usr
    ! since we restart from a potential field, we don't need usr_init_one_grid
    !usr_init_one_grid   => initonegrid_usr  
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
 
    ! prepare magnetogram at bottom
    if(firstusrglobaldata) then
      call init_b_nlfff_data('./nlfff_boundary/cutout.dat',unit_length,unit_magneticfield)
      firstusrglobaldata=.false.
    endif

  end subroutine initglobaldata_usr

  !==============================================================================
  ! Purpose: to read vector magnetic field from an input file
  !==============================================================================
  subroutine init_b_nlfff_data(magnetogramname,qLunit,qBunit)
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
      call mpistop("no input magnetogram----init_b_nlfff_data")
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
      print*,'extrapolating 3D nonlinear force-free field from an observed Brad,Bthe, and Bphi'
      print*,'magnetogram of',nx1,'by',nx2,'pixels, which includes the ghost layers.'
    endif
  end subroutine init_b_nlfff_data

  !==============================================================================
  ! Purpose: to initialize the initial condition
  !==============================================================================
  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S, 1:ndim)
  !  since we restart from a potential field, we don't need this subroutine 
  end subroutine initonegrid_usr

  !==============================================================================
  ! Purpose: to provide special boundary conditions set by users.
  !==============================================================================
  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters
    use mod_physics
    integer, intent(in) :: ixI^L, ixO^L, iB
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: drad(0:2),dthe,dphi
    double precision :: xlen^D,dxb^D,startpos^D
    integer :: ngl,ix^D,ixbc^D,af

    ngl=nghostcells*(2**(refine_max_level-1))
    drad(0)=1.d0
    drad(1)=drad(0) !*qstretch_baselevel(1)
    drad(2)=drad(1) !*qstretch_baselevel(1)
    dthe=1.d0
    dphi=1.d0

    select case(iB)
    case(1)    ! r mininum surface
      w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                    (1.d0/SRadius-1.d0/(x(ixO^S,1)+SRadius)))
      w(ixO^S,mom(1))=-w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))
      w(ixO^S,mom(2))=-w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))
      w(ixO^S,mom(3))=-w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))
      ! The first layer (with the resolution of the level 1 grid) are assigned with the fixed boundary data
      ! The outer ghost layers are assigned with the following one-sided difference scheme
      dxb^D=dx(^D,refine_max_level);
      af=nint(dxlevel(2)/dxb2)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        xlen^D=x({ix^D,},^D)-xprobmin^D+dble(ngl)*dxb^D;\
        if(af>=2) then
          call mpistop("Bottom boundary is not defined ")
        else
          ixbc^D=ceiling(xlen^D/dxb^D);\
          !if (ixbc1 <= 0 .or. ixbc2 <=0 .or. ixbc1 >=245 .or. ixbc2 >=245) then
          !  write(*,*) 'ixbc1,ixbc2:',ixbc1,ixbc2
          !endif
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
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))=(1.0d0/3.0d0)*     &
                   (-w(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))*    &
                    (x(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2*      &
                 sin(x(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*        &
                     dthe*dphi)                                          &
              +4.0d0*w(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))*    &
                    (x(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2*      &
                 sin(x(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*        &
                     dthe*dphi))/                                        &
                    (x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2*      &
                 sin(x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*        &
                     dthe*dphi)

        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2))=(1.0d0/3.0d0)*     &
                   (-w(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2))*    &
                    (x(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*         &
                 sin(x(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*        &
                     drad(2)*dphi)                                       &
              +4.0d0*w(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2))*    &
                    (x(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*         &
                 sin(x(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*        &
                     drad(1)*dphi))/                                     &
                    (x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*         &
                 sin(x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*        &
                     drad(0)*dphi)

        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3))=(1.0d0/3.0d0)*     &
                   (-w(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3))*    &
                    (x(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*         &
                     drad(2)*dthe)                                       &
              +4.0d0*w(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3))*    &
                    (x(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*         &
                     drad(1)*dthe))/                                     &
                    (x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*         &
                     drad(0)*dthe)
      enddo
    case(2)    ! r maxinum surface
      w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                    (1.d0/SRadius-1.d0/(x(ixO^S,1)+SRadius)))
      w(ixO^S,mom(1))=-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))
      w(ixO^S,mom(2))=-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))
      w(ixO^S,mom(3))=-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))
      ! 2nd order accuacy constant value extrapolation
      do ix1=ixOmin1,ixOmax1
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))=(1.0d0/3.0d0)*     &
                  (-w(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))*     &
                   (x(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2*       &
                sin(x(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*         &
                    dthe*dphi)                                           &
             +4.0d0*w(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))*     &
                   (x(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2*       &
                sin(x(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*         &
                    dthe*dphi))/                                         &
                   (x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2*       &
                sin(x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*         &
                   dthe*dphi)

        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2))=(1.0d0/3.0d0)*     &
                  (-w(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2))*     &
                   (x(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*          &
                sin(x(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*         &
                    drad(0)*dphi)                                        &
             +4.0d0*w(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2))*     &
                   (x(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*          &
                sin(x(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*         &
                    drad(1)*dphi))/                                      &
                   (x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*          &
                sin(x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*         &
                    drad(2)*dphi)

        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3))=(1.0d0/3.0d0)*     &
                  (-w(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3))*     &
                   (x(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*          &
                    drad(0)*dthe)                                        &
             +4.0d0*w(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3))*     &
                   (x(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*          &
                    drad(1)*dthe))/                                      &
                   (x(ix1+0,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*          &
                    drad(2)*dthe)
      enddo
    case(3)    ! theta mininum surface
      w(ixO^S,rho_)  = w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
      w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(1))
      w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(2))
      w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(3))
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
      w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(1))
      w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(2))
      w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(3))
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
      w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(1))
      w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(2)) 
      w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(3)) 
      ! 2nd order accuacy constant value extrapolation
      do ix3=ixOmax3,ixOmin3,-1
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(1):mag(3))=(1.0d0/3.0d0)* &
              (     -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+2,mag(1):mag(3)) &
               +4.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+1,mag(1):mag(3)))
      enddo
    case(6)    ! phi maxinum surface
      w(ixO^S,rho_)  = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,rho_)
      w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(1))
      w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(2))
      w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(3))
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
    if (any(x(ixO^S,3)<=xprobmin3+0.3d0)) then
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
    Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1):mag(ndir))
    B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)
    ! output Alfven wave speed B/sqrt(rho)
w(ixO^S,nw+1)=dsqrt(B2(ixO^S)/w(ixO^S,rho_))
    ! output divB1
    call divvector(Btotal,ixI^L,ixO^L,divb)
    w(ixO^S,nw+2)=0.5d0*divb(ixO^S)/dsqrt(B2(ixO^S))/(^D&1.0d0/dxlevel(^D)+)
    ! output the plasma beta p*2/B**2
    w(ixO^S,nw+3)=w(ixO^S,rho_)*Tiso*two/B2(ixO^S)
    ! store current
    w(ixO^S,nw+4)=w(ixO^S,mag(1))
    w(ixO^S,nw+5)=w(ixO^S,mag(2))
    w(ixO^S,nw+6)=w(ixO^S,mag(3))
    call curlvector(Btotal,ixI^L,ixO^L,curlvec,idirmin,1,ndir)
    do idir=1,ndir
      w(ixO^S,nw+6+idir)=curlvec(ixO^S,idir)
    end do
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

    varnames='Alfv divB beta br bth bph j1 j2 j3'
   !varnames='br bth bph j1 j2 j3'
  end subroutine specialvarnames_output

end module mod_usr





