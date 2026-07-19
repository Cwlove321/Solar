!> Magnetic bipolar field
module mod_usr
  use mod_mhd

  implicit none
  integer, parameter :: jmax=8000
  double precision, allocatable :: pbc(:),rbc(:)
  double precision :: t_driven, t_relax, a0,q_para,d_para,L_para,L1_para,d1_para,q1_para,v_0
  ! 1D solar atmosphere table for pressure, density, and height
  double precision :: pa(jmax),ra(jmax),ya(jmax)
  double precision :: usr_grav,SRadius,rhob,Tiso,dr,joverb_cr,gzone,bQ0,lQ0,startpos^D,tdstop
  ! variables for boundary condition
  double precision, allocatable, save :: Bbt(:,:,:,:),vbt(:,:,:,:)
  integer, save :: nxbc^D,ngl
  ! x,y,z coordinates of the flux rope axis (integral path)
  double precision, allocatable :: x_axis(:,:), f(:), s(:)
  integer :: np,ixp

contains

  !==============================================================================
  ! Purpose: to include global parameters, set user methods, set coordinate 
  !          system and activate physics module.
  !==============================================================================
  subroutine usr_init()

    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_refine_threshold=> specialthreshold
    usr_set_electric_field => boundary_electric_field
    usr_source => specialsource_density
    usr_transform_w    => transform_w_fill_plasma
    !usr_internal_bc     => driven
    !usr_var_for_errest      => p_for_errest


    call set_coordinate_system("Cartesian_3D")
    call mhd_activate()

  end subroutine usr_init

  !==============================================================================
  ! Purpose: to initialize user public parameters and reset global parameters.
  !          Input data are also read here.
  !==============================================================================
  subroutine initglobaldata_usr()
    double precision :: x0, y0, z0, r0, h, theta, L, Bperp
    integer :: ixp, nphalf

    !unit_density       = 1.4d0*mass_H*unit_numberdensity               ! 2.341668000000000E-015 g*cm^-3
    !unit_pressure      = 2.3d0*unit_numberdensity*k_B*unit_temperature ! 0.317538000000000 erg*cm^-3
    !unit_magneticfield = dsqrt(miu0*unit_pressure)                     ! 1.99757357615242 Gauss
    !unit_velocity      = unit_magneticfield/dsqrt(miu0*unit_density)   ! 1.16448846777562E007 cm/s = 116.45 km/s
    !unit_time          = unit_length/unit_velocity                     ! 85.8746159942810 s 
    rhob=10.d0
    ! isothermal uniform temperature
    Tiso=1.0d0

    !driven velocity
    t_relax = 5.0d0 
    t_driven= 20.0d0
    v_0= 7.5*1.0d5/unit_velocity
    !v_0=20.0d0*1.0d5/unit_velocity
    gzone=0.2d0
    ! cell size in 1D solar atmosphere table
    dr=(2.d0*gzone+xprobmax3-xprobmin3)/dble(jmax)
    usr_grav=-1.5*2.74d4*unit_length/unit_velocity**2  ! solar gravity
    SRadius=6.955d10/unit_length                   ! Solar radius
    d_para=1.0d9/unit_length ! depth of magnetic charges
    L_para=1.0d9/unit_length ! half distance between magnetic charges

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: xC(ixI^S,1:ndim)
    double precision :: Bf(ixI^S,1:ndir)
    double precision :: res,qalpha,llift
    integer :: ix^D,na,idir,ixC^L
    logical, save :: first=.true.

    w(ixO^S,mag(1:3))=w(ixO^S,mag(1:3))
    w(ixO^S,rho_)=w(ixO^S,rho_)
    w(ixO^S,mom(1:3))=0.0d0
    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  ! restart full mhd from isothermal mhd data and fill solar atmosphere
  subroutine transform_w_fill_plasma(ixI^L,ixO^L,nw_in,w_in,x,w_out)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, nw_in
    double precision, intent(in)  :: w_in(ixI^S,1:nw_in)
    double precision, intent(in)  :: x(ixI^S, 1:ndim)
    double precision, intent(out) :: w_out(ixI^S,1:nw)

    double precision :: d_f,res,Bf(ixI^S,1:ndir),xC(ixI^S,1:ndim)
    integer :: ix^D,na,idir,ixC^L

    w_out=0.d0
    ! only preserve magnetic field
    w_out(ixO^S,mag(1:3))=w_in(ixO^S,4:6)
    w_out(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                    (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
    !print *, rhob*dexp(usr_grav*SRadius**2/Tiso*&
    !                (1.d0/SRadius-1.d0/(xprobmax3+SRadius)))
    w_out(ixO^S,mom(1:3)) = 0.0d0
    call mhd_to_conserved(ixI^L,ixO^L,w_out,x)

  end subroutine transform_w_fill_plasma

  ! allow user to change inductive electric field, especially for boundary driven applications
  subroutine boundary_electric_field(ixI^L,ixO^L,qt,qdt,fE,s)
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt,qdt
    type(state)                        :: s
    double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)

    integer :: ixC^L,ixA^L


    if (qt > t_driven+t_relax .or. qt < t_relax) then
      if(s%is_physical_boundary(5)) then
        ixCmin^D=ixOmin^D-kr(3,^D);
        ixCmax^D=ixOmax^D;
        ixAmin^D=ixCmin^D-kr(2,^D);
        ixAmax^D=ixCmax^D;
        fE(ixAmin3^%3ixA^S,1)=0.d0
        ixAmin^D=ixCmin^D-kr(1,^D);
        ixAmax^D=ixCmax^D;
        fE(ixAmin3^%3ixA^S,2)=0.d0
      end if
    endif


    ! fix normal B at side boundaries
    if(s%is_physical_boundary(1)) then
      ixCmin^D=ixOmin^D-kr(1,^D);
      ixCmax^D=ixOmax^D;
      ixAmin^D=ixCmin^D-kr(3,^D);
      ixAmax^D=ixCmax^D;
      fE(ixAmin1^%1ixA^S,2)=0.d0
      ixAmin^D=ixCmin^D-kr(2,^D);
      ixAmax^D=ixCmax^D;
      fE(ixAmin1^%1ixA^S,3)=0.d0
      if(s%is_physical_boundary(5)) then
        ixCmin^D=ixOmin^D-kr(3,^D);
        ixCmax^D=ixOmax^D;
        ixAmin^D=ixCmin^D-kr(2,^D);
        ixAmax^D=ixCmax^D-kr(1,^D);
        fE(ixAmin3^%3ixA^S,1)=0.d0
        ixAmin^D=ixCmin^D-kr(1,^D);
        ixAmax^D=ixCmax^D-kr(1,^D);
        fE(ixAmin3^%3ixA^S,2)=0.d0
      end if
    end if
    if(s%is_physical_boundary(2)) then
      ixCmin^D=ixOmin^D-kr(1,^D);
      ixCmax^D=ixOmax^D;
      ixAmin^D=ixCmin^D-kr(3,^D);
      ixAmax^D=ixCmax^D;
      fE(ixAmax1^%1ixA^S,2)=0.d0
      ixAmin^D=ixCmin^D-kr(2,^D);
      ixAmax^D=ixCmax^D;
      fE(ixAmax1^%1ixA^S,3)=0.d0
      if(s%is_physical_boundary(5)) then
        ixCmin^D=ixOmin^D-kr(3,^D);
        ixCmax^D=ixOmax^D;
        ixAmin^D=ixCmin^D-kr(2,^D)+kr(1,^D);
        ixAmax^D=ixCmax^D;
        fE(ixAmin3^%3ixA^S,1)=0.d0
        ixAmin^D=ixCmin^D;
        ixAmax^D=ixCmax^D;
        fE(ixAmin3^%3ixA^S,2)=0.d0
      end if
    end if
    if(s%is_physical_boundary(3)) then
      ixCmin^D=ixOmin^D-kr(2,^D);
      ixCmax^D=ixOmax^D;
      ixAmin^D=ixCmin^D-kr(3,^D);
      ixAmax^D=ixCmax^D;
      fE(ixAmin2^%2ixA^S,1)=0.d0
      ixAmin^D=ixCmin^D-kr(1,^D);
      ixAmax^D=ixCmax^D;
      fE(ixAmin2^%2ixA^S,3)=0.d0
      if(s%is_physical_boundary(5)) then
        ixCmin^D=ixOmin^D-kr(3,^D);
        ixCmax^D=ixOmax^D;
        ixAmin^D=ixCmin^D-kr(2,^D);
        ixAmax^D=ixCmax^D-kr(2,^D);
        fE(ixAmin3^%3ixA^S,1)=0.d0
        ixAmin^D=ixCmin^D-kr(1,^D);
        ixAmax^D=ixCmax^D-kr(2,^D);
        fE(ixAmin3^%3ixA^S,2)=0.d0
      end if
    end if

    if(s%is_physical_boundary(4)) then
      ixCmin^D=ixOmin^D-kr(2,^D);
      ixCmax^D=ixOmax^D;
      ixAmin^D=ixCmin^D-kr(3,^D);
      ixAmax^D=ixCmax^D;
      fE(ixAmax2^%2ixA^S,1)=0.d0
      ixAmin^D=ixCmin^D-kr(1,^D);
      ixAmax^D=ixCmax^D;
      fE(ixAmax2^%2ixA^S,3)=0.d0
      if(s%is_physical_boundary(5)) then
        ixCmin^D=ixOmin^D-kr(3,^D);
        ixCmax^D=ixOmax^D;
        ixAmin^D=ixCmin^D;
        ixAmax^D=ixCmax^D;
        fE(ixAmin3^%3ixA^S,1)=0.d0
        ixAmin^D=ixCmin^D-kr(1,^D)+kr(2,^D);
        ixAmax^D=ixCmax^D;
        fE(ixAmin3^%3ixA^S,2)=0.d0
      end if
    end if

  end subroutine boundary_electric_field

  subroutine driven(level,qt,ixI^L,ixO^L,w,x)
    integer, intent(in) :: ixI^L,ixO^L,level
    double precision, intent(in) :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    type(state)                        :: s
    double precision :: mu_density, density0 
    integer :: ix^D


    {do ix^DB=ixOmin^DB,ixOmax^DB\}
       density0=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/(x(ix^D,3)+SRadius)))
       if (w(ix^D,p_)<0) then
          w(ix^D,rho_)=density0
          w(ix^D,p_)=density0
        endif
    {end do\}

  end subroutine driven 

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: tmp1(ixI^S),tmp2(ixI^S),pth(ixI^S),Qp(ixI^S)
    double precision :: xlen^D,coeffrho,ft,tfstop,tramp1,tramp2,vlimit
    integer :: idir,ix^D,ixIM^L,ixOs^L,ixbc^D

    select case(iB)
     case(1)
       w(ixO^S,mom(1))=w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))/&
                        w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       w(ixO^S,mom(2))=w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))/&
                        w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       w(ixO^S,mom(3))=w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))/&
                        w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       if(stagger_grid) then
         do idir=1,nws
           if(idir==1) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix1=ixOsmax1,ixOsmin1,-1
             block%ws(ix1^%1ixOs^S,idir) = 1.d0/3.d0*&
                    (-block%ws(ix1+2^%1ixOs^S,idir)&
                +4.d0*block%ws(ix1+1^%1ixOs^S,idir))
             !block%ws(ix1^%1ixOs^S,idir) = 1.d0/3.d0*&
             !       ( block%ws(ix1+3^%1ixOs^S,idir)&
             !   -5.d0*block%ws(ix1+2^%1ixOs^S,idir)&
             !   +7.d0*block%ws(ix1+1^%1ixOs^S,idir))
             !block%ws(ix1^%1ixOs^S,idir) = 1.d0/11.d0*&
             !    ( -2.d0*block%ws(ix1+4^%1ixOs^S,idir)&
             !     +11.d0*block%ws(ix1+3^%1ixOs^S,idir)&
             !     -27.d0*block%ws(ix1+2^%1ixOs^S,idir)&
             !     +29.d0*block%ws(ix1+1^%1ixOs^S,idir))
           end do
         end do
         ixOs^L=ixO^L-kr(1,^D);
         block%ws(ixOs^S,1)=zero
         do ix1=ixOsmax1,ixOsmin1,-1
           call get_divb(w,ixI^L,ixO^L,Qp)
           block%ws(ix1^%1ixOs^S,1)=Qp(ix1+1^%1ixO^S)*block%dvolume(ix1+1^%1ixO^S)&
             /block%surfaceC(ix1^%1ixOs^S,1)
         end do
         call mhd_face_to_center(ixO^L,block)
       else
         ! 2nd order accuacy constant value extrapolation
         do ix1=ixOmax1,ixOmin1,-1
           w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                      (-w(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)) &
                 +4.0d0*w(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)))
         enddo
       end if
       if(mhd_energy) then
         !ixIM^L=ixO^L;
         !ixIMmin1=ixOmax1+1;ixIMmax1=ixOmax1+nghostcells;
         !call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
         !w(ixO^S,rho_)=w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
         !w(ixO^S,p_)=pth(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
         w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
         w(ixO^S,p_)= w(ixO^S,rho_)*Tiso
       else if(mhd_adiab==0) then
         ! zero beta
         w(ixO^S,rho_)=w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       else
         ! isothermal
         w(ixO^S,rho_)=w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       end if
       call mhd_to_conserved(ixI^L,ixO^L,w,x)
     case(2)
       w(ixO^S,mom(1))=w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))/&
                        w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       w(ixO^S,mom(2))=w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))/&
                        w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       w(ixO^S,mom(3))=w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))/&
                        w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       if(stagger_grid) then
         do idir=1,nws
           if(idir==1) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix1=ixOsmin1,ixOsmax1
             block%ws(ix1^%1ixOs^S,idir) = 1.d0/3.d0*&
                    (-block%ws(ix1-2^%1ixOs^S,idir)&
                +4.d0*block%ws(ix1-1^%1ixOs^S,idir))
             !block%ws(ix1^%1ixOs^S,idir) = 1.d0/3.d0*&
             !       ( block%ws(ix1-3^%1ixOs^S,idir)&
             !   -5.d0*block%ws(ix1-2^%1ixOs^S,idir)&
             !   +7.d0*block%ws(ix1-1^%1ixOs^S,idir))
             !block%ws(ix1^%1ixOs^S,idir) = 1.d0/11.d0*&
             !    ( -2.d0*block%ws(ix1-4^%1ixOs^S,idir)&
             !     +11.d0*block%ws(ix1-3^%1ixOs^S,idir)&
             !     -27.d0*block%ws(ix1-2^%1ixOs^S,idir)&
             !     +29.d0*block%ws(ix1-1^%1ixOs^S,idir))
           end do
         end do
         ixOs^L=ixO^L;
         block%ws(ixOs^S,1)=zero
         do ix1=ixOsmin1,ixOsmax1
           call get_divb(w,ixI^L,ixO^L,Qp)
           block%ws(ix1^%1ixOs^S,1)=-Qp(ix1^%1ixO^S)*block%dvolume(ix1^%1ixO^S)&
             /block%surfaceC(ix1^%1ixOs^S,1)
         end do
         call mhd_face_to_center(ixO^L,block)
       else
         ! 2nd order accuacy constant value extrapolation
         do ix1=ixOmin1,ixOmax1
           w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                      (-w(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)) &
                 +4.0d0*w(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)))
         enddo
       end if
       if(mhd_energy) then
         !ixIM^L=ixO^L;
         !ixIMmin1=ixOmin1-nghostcells;ixIMmax1=ixOmin1-1;
         !call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
         !w(ixO^S,rho_)=w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
         !w(ixO^S,p_)=pth(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
         w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
         w(ixO^S,p_)= w(ixO^S,rho_)*Tiso
       else if(mhd_adiab==0) then
         ! zero beta
         w(ixO^S,rho_)=w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       else
         !w(ixO^S,rho_)=w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
         w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
       end if
       call mhd_to_conserved(ixI^L,ixO^L,w,x)
     case(3)
       w(ixO^S,mom(1))=w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(1))/&
                        w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
       w(ixO^S,mom(2))=w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(2))/&
                        w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
       w(ixO^S,mom(3))=w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(3))/&
                        w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
       if(stagger_grid) then
         do idir=1,nws
           if(idir==2) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix2=ixOsmax2,ixOsmin2,-1
             block%ws(ix2^%2ixOs^S,idir) = 1.d0/3.d0*&
                    (-block%ws(ix2+2^%2ixOs^S,idir)&
                +4.d0*block%ws(ix2+1^%2ixOs^S,idir))
             !block%ws(ix2^%2ixOs^S,idir) = 1.d0/3.d0*&
             !       ( block%ws(ix2+3^%2ixOs^S,idir)&
             !   -5.d0*block%ws(ix2+2^%2ixOs^S,idir)&
             !   +7.d0*block%ws(ix2+1^%2ixOs^S,idir))
             !block%ws(ix2^%2ixOs^S,idir) = 1.d0/11.d0*&
             !    ( -2.d0*block%ws(ix2+4^%2ixOs^S,idir)&
             !     +11.d0*block%ws(ix2+3^%2ixOs^S,idir)&
             !     -27.d0*block%ws(ix2+2^%2ixOs^S,idir)&
             !     +29.d0*block%ws(ix2+1^%2ixOs^S,idir))
           end do
         end do
         ixOs^L=ixO^L-kr(2,^D);
         block%ws(ixOs^S,2)=zero
         do ix2=ixOsmax2,ixOsmin2,-1
           call get_divb(w,ixI^L,ixO^L,Qp)
           block%ws(ix2^%2ixOs^S,2)=Qp(ix2+1^%2ixO^S)*block%dvolume(ix2+1^%2ixO^S)&
             /block%surfaceC(ix2^%2ixOs^S,2)
         end do
         call mhd_face_to_center(ixO^L,block)
       else
         ! 2nd order accuacy constant value extrapolation
         do ix2=ixOmax2,ixOmin2,-1
           w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                      (-w(ixOmin1:ixOmax1,ix2+2,ixOmin3:ixOmax3,mag(1):mag(3)) &
                 +4.0d0*w(ixOmin1:ixOmax1,ix2+1,ixOmin3:ixOmax3,mag(1):mag(3)))
         enddo
       end if
       if(mhd_energy) then
         !ixIM^L=ixO^L;
         !ixIMmin2=ixOmax2+1;ixIMmax2=ixOmax2+nghostcells;
         !call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
         !w(ixO^S,p_)=pth(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3)
         !w(ixO^S,rho_)=w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
         w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
         w(ixO^S,p_)= w(ixO^S,rho_)*Tiso
       else if(mhd_adiab==0) then
         ! zero beta
         w(ixO^S,rho_)=w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
       else
         !w(ixO^S,rho_)=w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
         w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
       end if
       call mhd_to_conserved(ixI^L,ixO^L,w,x)
     case(4)
       w(ixO^S,mom(1))=w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(1))/&
                        w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,rho_)
       w(ixO^S,mom(2))=w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(2))/&
                        w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,rho_)
       w(ixO^S,mom(3))=w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(3))/&
                        w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,rho_)
       if(stagger_grid) then
         do idir=1,nws
           if(idir==2) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix2=ixOsmin2,ixOsmax2
             block%ws(ix2^%2ixOs^S,idir) = 1.d0/3.d0*&
                    (-block%ws(ix2-2^%2ixOs^S,idir)&
                +4.d0*block%ws(ix2-1^%2ixOs^S,idir))
             !block%ws(ix2^%2ixOs^S,idir) = 1.d0/3.d0*&
             !       ( block%ws(ix2-3^%2ixOs^S,idir)&
             !   -5.d0*block%ws(ix2-2^%2ixOs^S,idir)&
             !   +7.d0*block%ws(ix2-1^%2ixOs^S,idir))
             !block%ws(ix2^%2ixOs^S,idir) = 1.d0/11.d0*&
             !    ( -2.d0*block%ws(ix2-4^%2ixOs^S,idir)&
             !     +11.d0*block%ws(ix2-3^%2ixOs^S,idir)&
             !     -27.d0*block%ws(ix2-2^%2ixOs^S,idir)&
             !     +29.d0*block%ws(ix2-1^%2ixOs^S,idir))
           end do
         end do
         ixOs^L=ixO^L;
         block%ws(ixOs^S,2)=zero
         do ix2=ixOsmin2,ixOsmax2
           call get_divb(w,ixI^L,ixO^L,Qp)
           block%ws(ix2^%2ixOs^S,2)=-Qp(ix2^%2ixO^S)*block%dvolume(ix2^%2ixO^S)&
             /block%surfaceC(ix2^%2ixOs^S,2)
         end do
         call mhd_face_to_center(ixO^L,block)
       else
         ! 2nd order accuacy constant value extrapolation
         do ix2=ixOmin2,ixOmax2
           w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                       (-w(ixOmin1:ixOmax1,ix2-2,ixOmin3:ixOmax3,mag(1):mag(3))&
                  +4.0d0*w(ixOmin1:ixOmax1,ix2-1,ixOmin3:ixOmax3,mag(1):mag(3)))
         enddo
       end if
       if(mhd_energy) then
         !ixIM^L=ixO^L;
         !ixIMmin2=ixOmin2-nghostcells;ixIMmax2=ixOmin2-1;
         !call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
         !w(ixO^S,p_)=pth(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3)
         !w(ixO^S,rho_)=w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,rho_)
         w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
         w(ixO^S,p_)= w(ixO^S,rho_)*Tiso
       else if(mhd_adiab==0) then
         ! zero beta
         w(ixO^S,rho_)=w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,rho_)
       else
         !w(ixO^S,rho_)=w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,rho_)
         w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
       end if
       call mhd_to_conserved(ixI^L,ixO^L,w,x)
     case(5)
       if(stagger_grid) then
         do idir=1,nws
           if(idir==3) cycle
             ixOsmax^D=ixOmax^D;
             ixOsmin^D=ixOmin^D-kr(^D,idir);
             do ix3=ixOsmax3,ixOsmin3,-1
               ! 3 one-sided equal gradient
               block%ws(ix3^%3ixOs^S,idir) = 1.d0/11.d0*&
                  ( -2.d0*block%ws(ix3+4^%3ixOs^S,idir)&
                   +11.d0*block%ws(ix3+3^%3ixOs^S,idir)&
                   -27.d0*block%ws(ix3+2^%3ixOs^S,idir)&
                   +29.d0*block%ws(ix3+1^%3ixOs^S,idir))
             end do
         end do
         ixOs^L=ixO^L-kr(3,^D);
         block%ws(ixOs^S,3)=zero
         do ix3=ixOsmax3,ixOsmin3,-1
           call get_divb(w,ixI^L,ixO^L,Qp)
           block%ws(ix3^%3ixOs^S,3)=Qp(ix3+1^%3ixO^S)*block%dvolume(ix3+1^%3ixO^S)&
             /block%surfaceC(ix3^%3ixOs^S,3)
         end do
         call mhd_face_to_center(ixO^L,block)
       else
         ! 4th order accuacy constant value extrapolation
         do ix3=ixOmax3,ixOmin3,-1
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(1):mag(3))=(1.0d0/25.0d0)* &
                ( -3.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+4,mag(1):mag(3)) &
                 +16.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+3,mag(1):mag(3)) &
                 -36.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+2,mag(1):mag(3)) &
                 +48.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+1,mag(1):mag(3)))
         enddo
       end if

         w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(1))/&
                          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,rho_)
         w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(2))/& 
                          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,rho_)
         w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(3))/& 
                          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,rho_)
  

         do ix3 = ixOmax3, ixOmax3
           do ix2 = ixOmin2, ixOmax2
              do ix1 = ixOmin1, ixOmax1
                if (abs(x(ix^D,1))< 1 .and. abs(x(ix^D,2))<1) then
                 w(ix^D,mom(1))=-1.0d0*v_0*sin(dpi*x(ix^D,1)/L_para)
                 w(ix^D,mom(2))=0.0d0
                 w(ix^D,mom(3))=0.0d0
                end if
              end do
           end do
        end do
   
         do ix3=ixOmax3-1,ixOmin3,-1
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mom(1:3))= &
                 0.12d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+5,mom(1:3)) &
                -0.76d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+4,mom(1:3)) &
                +2.08d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+3,mom(1:3)) &
                -3.36d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+2,mom(1:3)) &
                +2.92d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+1,mom(1:3))
         enddo
 
        if (qt > t_driven+t_relax .or. qt < t_relax) then
         w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(1))/&
                          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,rho_)
         w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(2))/& 
                          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,rho_)
         w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(3))/& 
                          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,rho_)
        endif 

       if(mhd_energy) then
         w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
         w(ixO^S,p_)= w(ixO^S,rho_)*Tiso
       else if(mhd_adiab==0) then
         ! zero beta
         w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
       else
         w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
       end if
       call mhd_to_conserved(ixI^L,ixO^L,w,x)
     case(6)
       if(stagger_grid) then
         do idir=1,nws
           if(idir==3) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix3=ixOsmin3,ixOsmax3
             block%ws(ix3^%3ixOs^S,idir) = 1.d0/3.d0*&
                    ( block%ws(ix3-3^%3ixOs^S,idir)&
                -5.d0*block%ws(ix3-2^%3ixOs^S,idir)&
                +7.d0*block%ws(ix3-1^%3ixOs^S,idir))
             !block%ws(ix3^%3ixOs^S,idir) = 1.d0/11.d0*&
             !   (+2.d0*block%ws(ix3-3^%3ixOs^S,idir)&
             !    -9.d0*block%ws(ix3-2^%3ixOs^S,idir)&
             !   +18.d0*block%ws(ix3-1^%3ixOs^S,idir))
           end do
         end do
         ixOs^L=ixO^L;
         block%ws(ixOs^S,3)=zero
         do ix3=ixOsmin3,ixOsmax3
           call get_divb(w,ixI^L,ixO^L,Qp)
           block%ws(ix3^%3ixOs^S,3)=-Qp(ix3^%3ixO^S)*block%dvolume(ix3^%3ixO^S)&
             /block%surfaceC(ix3^%3ixOs^S,3)
         end do
         call mhd_face_to_center(ixO^L,block)
       else
         ! 2nd order accuacy constant value extrapolation
         do ix3=ixOmin3,ixOmax3
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                 (     -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-2,mag(1):mag(3)) &
                  +4.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-1,mag(1):mag(3)))
         enddo
       end if
       if(mhd_energy) then
         ixIM^L=ixO^L;
         ixIMmin3=ixOmin3-1;ixIMmax3=ixOmax3;
         call getggrav(tmp1,ixI^L,ixIM^L,x)
         ixIMmin3=ixOmin3-1;ixIMmax3=ixOmin3-1;
         call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
         pth(ixOmin3-1^%3ixO^S)=pth(ixOmin3-1^%3ixO^S)/w(ixOmin3-1^%3ixO^S,rho_)
         where(pth(ixOmin3-1^%3ixO^S)<1.0d0)
           pth(ixOmin3-1^%3ixO^S)=1.0d0
         end where
         tmp2=0.d0
         do ix3=ixOmin3,ixOmax3
           tmp2(ixOmin3^%3ixO^S)=tmp2(ixOmin3^%3ixO^S)+0.5d0*(tmp1(ix3-1^%3ixO^S)+tmp1(ix3^%3ixO^S))/pth(ixOmin3-1^%3ixO^S)
           w(ix3^%3ixO^S,rho_)=w(ixOmin3-1^%3ixO^S,rho_)*dexp(tmp2(ixOmin3^%3ixO^S)*dxlevel(3))
           w(ix3^%3ixO^S,p_)=w(ix3^%3ixO^S,rho_)*pth(ixOmin3-1^%3ixO^S)
           where(w(ix3^%3ixO^S,mom(3))<0.d0)
            w(ix3^%3ixO^S,mom(3))=0.d0
          end where
         enddo
       else if(mhd_adiab==0) then
         ! zero beta
         coeffrho=usr_grav*SRadius**2/Tiso
         do ix3=ixOmin3,ixOmax3
           w(ix3^%3ixO^S,rho_)=w(ixOmin3-1^%3ixO^S,rho_)*dexp(coeffrho*(1.d0/(SRadius+&
           x(ixOmin3-1^%3ixO^S,3))-1.d0/(SRadius+x(ix3^%3ixO^S,3))))
         enddo
       else
         coeffrho=usr_grav*SRadius**2/Tiso
         do ix3=ixOmin3,ixOmax3
           w(ix3^%3ixO^S,rho_)=w(ixOmin3-1^%3ixO^S,rho_)*dexp(coeffrho*(1.d0/(SRadius+&
           x(ixOmin3-1^%3ixO^S,3))-1.d0/(SRadius+x(ix3^%3ixO^S,3))))
         enddo
       end if
       w(ixO^S,mom(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(1))/&
                       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,rho_)
       w(ixO^S,mom(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(2))/&
                       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,rho_)
       w(ixO^S,mom(3))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(3))/&
                       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,rho_)
       !w(ixO^S,mom(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(1))/&
       !                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,rho_)
       !w(ixO^S,mom(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(2))/&
       !                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,rho_)
       !w(ixO^S,mom(3))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(3))/&
       !                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,rho_)
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          if(w(ix^D,mom(3))<0) then
             w(ix^D,mom(3))=0.0d0
          end if
       {end do\}
       call mhd_to_conserved(ixI^L,ixO^L,w,x)
     case default
      call mpistop("Special boundary is not defined for this region")
    end select
  end subroutine specialbound_usr

  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    integer :: ix3
    do ix3=ixOmin3,ixOmax3
      ggrid(ix3^%3ixO^S)=usr_grav*(SRadius/(SRadius+x(ix3^%3ixO^S,3)))**2
    end do
  end subroutine

  !==============================================================================
  ! Purpose: get gravity field
  !==============================================================================
  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)
    double precision                :: ggrid(ixI^S)

    gravity_field=0.d0
    call getggrav(ggrid,ixI^L,ixO^L,x)
    gravity_field(ixO^S,3)=ggrid(ixO^S)
  end subroutine gravity

  subroutine specialsource_density(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: density_deviation(ixI^S), density_diffusion(ixI^S), gradd(ixI^S), &
                        graddensity(ixI^S,ndim),density0, mu_density, res
    integer :: ix^D, idim, ixA^L, na
    !double precision :: density0, mu_density, res
    !integer :: ix^D, na

    mu_density=0.1d0
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
       density0=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/(x(ix^D,3)+SRadius)))
       if (x(ix^D,3) < xprobmin3 + 1.0d0 .or. x(ix^D, 3) > xprobmax3 - 0.1d0) then
        if (w(ix^D,rho_)<mu_density*density0) then
          w(ix^D,rho_)=mu_density*density0
        endif
       endif
    {end do\}
  end subroutine specialsource_density


  !==============================================================================
  ! Purpose: Enforce additional refinement or coarsening. One can use the
  !          coordinate info in x and/or time qt=t_n and w(t_n) values w.
  !==============================================================================
subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen
    double precision  :: current(ixI^S,1:ndir), btotal(ixI^S,1:ndir), absj(ixO^S), absb(ixO^S)
    integer :: idirmin


    joverb_cr=1.0
    btotal = w(ixI^S,mag(:))
    call curlvector(btotal,ixI^L,ixO^L,current,idirmin,1,ndir)
    absj(ixO^S) = dsqrt(sum(current(ixO^S,:)**2,dim=ndim+1))
    absb(ixO^S) = dsqrt(sum(btotal(ixO^S,:)**2,dim=ndim+1))

     if ((level<refine_max_level) .and. any((absj/absb)>joverb_cr)) then
        refine=1
        coarsen=-1
     else if ((level<(refine_max_level-1)) .and. any((absj/absb)>(joverb_cr-0.3))) then
        refine=1
        coarsen=-1
      else if ((level<(refine_max_level-2)) .and. any((absj/absb)>(joverb_cr-0.5))) then
        refine=1
        coarsen=-1
      else if (level<(refine_max_level-3)) then
        refine=0
        coarsen=0
      else
        refine=-1
        coarsen=0
      end if
  end subroutine special_refine_grid


  !PURPOSE: use different tolerance in special regions for AMR to
  !reduce/increase resolution there where nothing/something interesting happens.
  subroutine specialthreshold(wlocal,xlocal,tolerance,qt,level)
    use mod_global_parameters
    
    integer, intent(in)          :: level
    double precision, intent(in) :: wlocal(1:nw),xlocal(1:ndim),qt
    double precision, intent(inout) :: tolerance
    
    double precision :: bczone^D,addtol,tol_add

    tol_add=0.d0
    !amplitude of additional tolerance
    addtol=0.3d0
    ! thickness of near-boundary region
    bczone1=1.5d0
    ! linear changing of additional tolerance
    if(xlocal(1)-xprobmin1 < bczone1 .or. xprobmax1-xlocal(1) < bczone1) then
      tol_add=(1.d0-min(xlocal(1)-xprobmin1,xprobmax1-xlocal(1))/bczone1)*addtol
    endif
    bczone2=1.5d0
    if(xlocal(2)-xprobmin2 < bczone2 .or. xprobmax2-xlocal(2) < bczone2) then
      tol_add=(1.d0-min(xlocal(2)-xprobmin2,xprobmax2-xlocal(2))/bczone2)*addtol
    endif
    bczone3=4.d0
    if(xprobmax3-xlocal(3) < bczone3) then
      tol_add=(1.d0-(xprobmax3-xlocal(3))/bczone3)*addtol
    endif
    tolerance=tolerance+tol_add

  end subroutine specialthreshold

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
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision                   :: tmp(ixI^S),dip(ixI^S),divb(ixI^S),B2(ixI^S),pth(ixI^S)
    double precision, dimension(ixI^S,1:ndir) :: Btotal,qvec,curlvec
    integer                            :: ix^D,idirmin,idims,idir,jdir,kdir

   ! Btotal & B^2
    Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))
    B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)
    ! output divB1
    call get_divb(w,ixI^L,ixO^L,divb)
    w(ixO^S,nw+1)=divb(ixO^S)
    ! store current
    call curlvector(Btotal,ixI^L,ixO^L,curlvec,idirmin,1,ndir)
    do idir=1,ndir
      w(ixO^S,nw+1+idir)=curlvec(ixO^S,idir)
    end do

  end subroutine specialvar_output

  !==============================================================================
  ! Purpose: names for special variable output
  !==============================================================================
  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames

      varnames='divB j1 j2 j3'

  end subroutine specialvarnames_output

end module mod_usr
