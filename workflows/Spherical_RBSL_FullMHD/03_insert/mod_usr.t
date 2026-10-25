! Spherical_RBSL_FullMHD workflow source, archived 2026-07-22.
! Includes PDFI_SS potential-field reconstruction, CT-consistent RBSL
! insertion, persistent staggered restart fields, optional filament loading,
! stretched-grid boundary stencils, dynamic AMR, and parameterized heating.
module mod_usr
  use mod_mhd
  implicit none
  !CTMF GLOBAL PARS
  !!Using V-B driven as Bh,V bot BC
  double precision, allocatable :: Avpt(:,:,:), Avpp(:,:,:)
  double precision, allocatable :: brll0(:,:)
  integer          :: m,n,p,ilevel,ic=0,factor
  integer          :: p_eth,p_eph,p_era
  double precision :: Rsrs,ssd
  !!Stretched mesh boudanry condition parameters
  double precision, allocatable :: ece4(:),ece3(:),ece2(:),ece1(:),ecz2(:),ecz1(:)
  !!Magnetic GLOBAL PARS
  logical                       :: do_rbsl=.false.,weak_field=.true.,I_Helix=.true.
  logical                       :: flip_axis=.true.
  integer                       :: ixgI^L,ixgO^L 
  integer                       :: npoints,np,nnp
  double precision              :: wf_scale
  double precision              :: Bflux,Rmfr
  double precision              :: F_flx,a0
  logical                       :: load_filament_density=.false.
  double precision              :: filament_density_factor=20.d0
  double precision              :: filament_core_fraction=0.6d0
  double precision              :: filament_min_height_Mm=5.d0
  logical                       :: liakh_mass_loading=.false.
  double precision              :: liakh_density_contrast=30.d0
  double precision              :: liakh_t_start_s=0.d0
  double precision              :: liakh_t_load_s=100.d0
  double precision              :: liakh_sigma_horizontal_Mm=3.6d0
  double precision              :: liakh_sigma_radial_Mm=4.8d0
  double precision              :: liakh_min_height_Mm=5.d0
  double precision              :: background_heating_cgs=0.d0
  double precision              :: dtheta,dphi
  double precision, allocatable :: br_rbsl(:,:)
  double precision, allocatable :: x_axis(:,:)
  double precision, allocatable :: liakh_axis_radius(:,:),liakh_axis_angle(:,:)
  !!MHD ATMOSPHERE PARS
  ! BEGIN ATMOSPHERE PARAMETERS
  ! Physical inputs are in Mm, K and cgs; defaults reproduce the successful RBSL15 atmosphere.
  character(len=32) :: atmosphere_model='analytic_tr'
  character(len=32) :: atmosphere_normalization='number_density'
  double precision :: atmosphere_reference_height_Mm=10.d0
  double precision :: atmosphere_reference_nH_cgs=6.2d8
  double precision :: atmosphere_reference_pressure_cgs=1.38365d2
  double precision :: atmosphere_chromosphere_K=8.d3
  double precision :: atmosphere_tanh_top_K=1.5d6
  double precision :: atmosphere_transition_height_Mm=2.d0
  double precision :: atmosphere_transition_width_Mm=0.25d0
  double precision :: atmosphere_transition_K=1.6d5
  double precision :: atmosphere_conductive_flux_cgs=2.d5
  double precision :: atmosphere_floor_height_Mm=-0.2d0
  double precision :: atmosphere_isothermal_K=1.d6
  ! END ATMOSPHERE PARAMETERS
  integer, parameter            :: jmax=8000
  double precision              :: ya(jmax),ra(jmax),pa(jmax)
  double precision              :: usr_grav,SRadius,rhob,Tiso,dr,gzone,rmax
  double precision, allocatable :: rbc(:),pbc(:)
  double precision              :: heatscale,bQ0,lQ0

contains

  subroutine usr_init()

    unit_length        = 6.955d10 !cm
    unit_numberdensity = 1.d9 !cm^-3
    unit_temperature   = 1.d6 !K

    usr_set_parameters        => initglobaldata_usr
    usr_init_one_grid         => initonegrid_usr
    usr_special_bc            => specialbound_usr
    usr_refine_grid           => special_refine_grid
    usr_aux_output            => specialvar_output
    usr_add_aux_names         => specialvarnames_output
    usr_refine_threshold      => specialthreshold
    usr_init_vector_potential => initvecpot_usr
    usr_set_electric_field    => boundary_electric_field
    usr_gravity               => gravity
    usr_source                => specialsource
    usr_improve_initial_condition => preserve_inserted_restart_field
!   usr_special_resistivity   => my_special_resistivity

    call set_coordinate_system('spherical_3D')
    call mhd_activate()

  end subroutine usr_init

  subroutine preserve_inserted_restart_field()
    ! The restart path already adds curl(A_RBSL) to the CT face field.
    ! Do not let the default initializer rebuild B from A_RBSL alone.
  end subroutine preserve_inserted_restart_field

  subroutine initglobaldata_usr 

    integer :: iv
    double precision :: qb,qc,qb23,qc23,q1b3,q1c3
    character(len=20) :: printsettingformat

    ! coefficients for gradient extrapolation in stretched grid
    if(stretched_dim(1)) then
      allocate(ece4(refine_max_level),ece3(refine_max_level),ece2(refine_max_level),ece1(refine_max_level))
      allocate(ecz2(refine_max_level),ecz1(refine_max_level))
      do iv=1,refine_max_level
        qb=1.d0+qstretch(iv,1)
        qc=qb+qstretch(iv,1)**2
        qb23=qb**2-qb**3
        qc23=qc**2-qc**3
        q1b3=1.d0-qb**3
        q1c3=1.d0-qc**3
        ece1(iv)=1.d0/(qstretch(iv,1)*(q1b3*qc23-q1c3*qb23))
        ece4(iv)=qb23*ece1(iv)
        ece3(iv)=-(qstretch(iv,1)*qb23+qc23)*ece1(iv)
        ece2(iv)=(qb**3*qc23-qc**3*qb23+qstretch(iv,1)*qc23)*ece1(iv)
        ece1(iv)=(q1b3*qc23-q1c3*qb23-qstretch(iv,1)*(qb**3*qc23-qc**3*qb23))*ece1(iv)
        ecz2(iv)=1.d0/(1.d0-qb**2)
        ecz1(iv)=-qb**2*ecz2(iv)
      end do
    end if
    !mf_eta=2.4693930682273452d-6

    printsettingformat='(1x,A20,ES15.7,A30)'
    ssd    = 1.d0
    Rsrs   = (1.d0+ssd)*6.955d10/unit_length
    ilevel = refine_max_level
    factor = nint(ssd/(xprobmax1-xprobmin1))
    ! ***********data resolution*******************8
    !unit_density       = 1.4d0*mass_H*unit_numberdensity               ! 2.341668000000000E-015 g*cm^-3
    !unit_pressure      = 2.3d0*unit_numberdensity*k_B*unit_temperature ! 0.317538000000000 erg*cm^-3
    !unit_magneticfield = dsqrt(miu0*unit_pressure)                     ! 1.99757357615242 Gauss
    !unit_velocity      = unit_magneticfield/dsqrt(miu0*unit_density)   ! 1.16448846777562E007 cm/s = 116.45 km/s
    !unit_time          = unit_length/unit_velocity                     ! about 5.973d3 s


    ! 3D Colat-lon system (tp CE) : p * m * n
    p = domain_nx1*(2**(refine_max_level-1))*factor  ! p cells from solar surface to source surface [1.0rsun -> 2.0rsun]
    m = domain_nx2*2**(refine_max_level-1)           ! Colat-lon system (theta) 
    n = domain_nx3*2**(refine_max_level-1)           ! Colat-lon system (phi) 
    p_eth = m*(n+1)
    p_eph = (m+1)*n
    p_era = (m+1)*(n+1)

    ixgOmin1=1;
    ixgOmin2=1;
    ixgOmin3=1;
    ixgOmax1=1;
    ixgOmax2=m+1;
    ixgOmax3=n+1;
    ixgI^L=ixgO^L^LADD1;

    allocate(brll0(m,n))
    allocate(Avpt(m,n+1,p+1))
    allocate(Avpp(m+1,n,p+1))


    ! input data file name
    call usr_params_read(par_files)
    if(do_rbsl.or.liakh_mass_loading) then
      nnp    = np*2
      dtheta = (xprobmax2-xprobmin2)/m
      dphi   = (xprobmax3-xprobmin3)/n
      allocate(x_axis(nnp,ndim))
      call flux_rope_insertion(x_axis)
      if(do_rbsl) allocate(br_rbsl(m,n))
      if(liakh_mass_loading) call prepare_liakh_mass_source_geometry()
    end if

    if(mype==0) then
      write(*,*) "m,n,p",m,n,p
      write(*,*) "Dimensionless units:"
      write(*,printsettingformat) "unit_length ",unit_length,"(cm)"
      write(*,printsettingformat) "unit_density ",unit_density,"(g cm^-3)"
      write(*,printsettingformat) "unit_temperature",unit_temperature,"(K)"
      write(*,*) "Deduced dimensionless units:"
      write(*,printsettingformat) "unit_time ",unit_time,"(s)"
      write(*,printsettingformat) "unit_velocity ",unit_velocity,"(cm s^-1)"
      write(*,printsettingformat)  &
        "unit_pressure ",unit_pressure,"(Ba = g cm^-1 s^-2)"
      write(*,printsettingformat)  &
        "unit_magneticfield ",unit_magneticfield,"(G=g^{1/2} cm^{-1/2} s^-1)"
    end if

    if(restart_from_file == undefined) call extrapolate_potential_field()

    !Setting of density and pressure
    if(.not.mhd_energy) then
      rhob = 1.d0       ! bottom density
      Tiso = mhd_adiab  ! isothermal uniform temperature
    end if
        rmax = xprobmax1 - xprobmin1
       gzone = nghostcells*rmax/(domain_nx1*2**(refine_max_level-1)) ! ghost cell size in 1D solar atmosphere table
!      gzone = 2.d7/unit_length                      ! ghost cell size in 1D solar atmosphere table
          dr = (2.d0*gzone + rmax)/dble(jmax)
    usr_grav = -2.74d4*unit_length/unit_velocity**2  ! solar gravity constant
     SRadius =  6.955d10/unit_length                 ! Solar radius
   heatscale = 6.d9/unit_length                      ! decay scale of background heating [6x10^9 cm]
         bQ0 = background_heating_cgs/unit_pressure*unit_time

    if(mhd_energy) call inithdstatic


  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    !initialize one grid within ixI^L
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    logical, save:: first=.true., last=.true.
    !MHD PARS
    integer :: ix^D,na,ixGs^L
    double precision :: res
    double precision :: pth(ixI^S)
    double precision, allocatable      :: wstmp(:,:,:,:)
    !Filament Insertion
    type(state)                        :: block0
    integer                            :: idirmin,idir,ixCT^L
    double precision                   :: wtmp(ixI^S,1:nw)
    double precision                   :: j_crit,eps_b,rho_c
    double precision                   :: current(ixI^S,3)
    double precision                   :: j_norm(ixI^S),dip(ixI^S),B2(ixI^S),tmp(ixI^S),rho_o(ixI^S)
    double precision                   :: bczone^D


    if(mype==0 .and. first) then
      print*,'initialize filament channel in sphere'
      first=.false.
    end if

    w(ixO^S,mom(:))=zero
    !w(ixO^S,mom(:))=w(ixO^S,mom(:))*0.1

    if(restart_from_file /= undefined.and.firstprocess.and.iprob==2) then
      if(stagger_grid) then
        !block0=block
        !ixCT^L=ixO^L^LADD1;
        !call rbsl_insert(block0%ixGs^L,ixI^L,ixCT^L,block0%ws,x,block0)
        !call mhd_face_to_center(ixO^L,block0)
        ixGs^L=block%ixGs^L;
        allocate(wstmp(ixGs^S,nws))
        wstmp(ixGs^S,1:nws) = block%ws(ixGs^S,1:nws)
        call rbsl_insert(block%ixGs^L,ixI^L,ixO^L,wstmp,x,block)
        block%ws(ixGs^S,1:nws) = wstmp(ixGs^S,1:nws)
        call fake_face_to_center(block%ixGs^L,ixO^L,wstmp,block)
        !call rbsl_insert(block%ixGs^L,ixI^L,ixCT^L,wstmp,x,block)
        !call fake_face_to_center(block%ixGs^L,ixCT^L,wstmp,block)
        deallocate(wstmp)
      end if
      !where(w(ixO^S,rho_)<=small_density)
      !      w(ixO^S,rho_)=small_density
      !      w(ixO^S,p_)  =small_pressure
      !end where
      !rho_o(ixO^S) = w(ixO^S,rho_)
            
      !print*,'mhd_gamma',mhd_gamma
      w(ixO^S,p_) = (mhd_gamma-1.d0)*w(ixO^S,e_)

      if(load_filament_density) then
        call enhance_filament_density(ixI^L,ixO^L,w,x)
      end if

      !!! find where the filament channel
      !wtmp = block%w
      !call get_current(wtmp,ixI^L,ixO^L,idirmin,current) 
      !j_norm(ixO^S) = dsqrt(sum(current(ixO^S,:)**2,ndim+1))
      !j_crit = 2400.d0
      !! find magnetic dips                                     
      !dip = 0.d0
      !do idir=1,ndir
      !  call gradient(wtmp(ixI^S,mag(idir)),ixI^L,ixO^L,1,tmp)
      !  call gradientS(wtmp(ixI^S,mag(idir)),ixI^L,ixM^LL,1,tmp)
      !  dip(ixO^S) = dip(ixO^S) + wtmp(ixO^S,mag(idir))*tmp(ixO^S)  !calculate B dot grad(B)
      !end do
      !! find filament threads angle scale [where |dtheta|<eps. log10(rho/rhoc)=-1]
      !B2(ixO^S) = sum((wtmp(ixO^S,mag(:)))**2,dim=ndim+1)
      !B2(ixO^S) = wtmp(ixO^S,mag(1))/dsqrt(B2(ixO^S))
      !!   eps_b = 0.486765
      !! find where is above the chromosphere [rho <= 1.d11]
      !    rho_c = 1.d11/unit_numberdensity
      !! thickness of near-boundary region
      !bczone1=0.3d0*(xprobmax1-xprobmin1)
      !bczone2=0.15d0*(xprobmax2-xprobmin2)
      !bczone3=0.15d0*(xprobmax3-xprobmin3)
      !where(w(ixO^S,rho_)<=rho_c .and. dip(ixO^S)>=0.d0 .and. j_norm(ixO^S)>=j_crit)
      !! w(ixO^S,rho_)=1.d11/unit_numberdensity*(10.0d0**(2.d0*exp(-dabs(B2(ixO^S))/eps_b*dlog(2.d0))))
      !  w(ixO^S,rho_)=rho_o(ixO^S)*200.0d0
      !end where

    else
      if(B0field) then
        w(ixO^S,mag(:))=0.d0
      else 
        if(stagger_grid) then
          call b_from_vector_potential(block%ixGs^L,ixI^L,ixO^L,block%ws,x)
          call mhd_face_to_center(ixO^L,block)
        else
          w(ixO^S,mag(:))=1.d0
        end if
      end if

      if(mhd_energy) then
        {do ix^DB=ixOmin^DB,ixOmax^DB\}
           na=floor((x(ix^D,1)-xprobmin1+gzone)/dr+0.5d0)
           res=x(ix^D,1)-xprobmin1+gzone-(dble(na)-0.5d0)*dr
           w(ix^D,rho_)=ra(na)+(one-cos(dpi*res/dr))/two*(ra(na+1)-ra(na))
           w(ix^D,p_)  =pa(na)+(one-cos(dpi*res/dr))/two*(pa(na+1)-pa(na))
        {end do\}
      else if(mhd_adiab/=0) then
        ! isothermal
        ! rho/rhob = exp(g rsun^2/gamma(1/rsun-1/r))
        ! ln(rho/rhob) = g rsun^2/gamma(1/rsun-1/r)
        ! drho/rho = (g rsun^2/gamma)*(dr/r^2)
        ! drho/rho = g(r)/gamma*(dr)
        w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                      (1.d0/SRadius-1.d0/x(ixO^S,1)))
      else
        ! zero beta
        w(ixO^S,rho_)=sum(block%w(ixO^S,mag(:))**2,dim=ndim+1)
      end if
    end if

    call mhd_to_conserved(ixI^L,ixO^L,w,x)
    
    if(mype==0 .and. last) then
      print*,'initial phase OK'
      last=.false.
    end if

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    
    double precision :: Qp(ixI^S),xlen^D
    double precision, allocatable :: magCC(:,:,:),xex(:,:,:)
    integer :: ix^D,ixOs^L,jxO^L,idir,ixbc^D,ixA^L
    !! MHD extra parameters
    integer :: ixIM^L
    double precision :: tmp1(ixI^S),tmp2(ixI^S),pth(ixI^S)
    double precision :: coeffrho

    select case(iB)
     case(1)
       !! Non-flow walls
       w(ixO^S,mom(:))=0.d0
       !{do ix^DB=ixOmin^DB,ixOmax^DB\}
       !    ! upflow continuous
       !    if(w(ixOmax1+1,ix2,ix3,mom(1))>zero) w(ix^D,mom(1))=w(ixOmax1+1,ix2,ix3,mom(1))
       !{end do\}
       if(stagger_grid) then
         do idir=1,nws
           if(idir==1) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           if(stretched_dim(1)) then
             do ix1=ixOsmax1,ixOsmin1,-1
             ! 3 one-sided equal gradient in stretched grid
               block%ws(ix1^%1ixOs^S,idir) = &
                   ece4(block%level)*block%ws(ix1+4^%1ixOs^S,idir)+&
                   ece3(block%level)*block%ws(ix1+3^%1ixOs^S,idir)+&
                   ece2(block%level)*block%ws(ix1+2^%1ixOs^S,idir)+&
                   ece1(block%level)*block%ws(ix1+1^%1ixOs^S,idir)
             end do
           else
             do ix1=ixOsmax1,ixOsmin1,-1
                ! 3 one-sided equal gradient
                block%ws(ix1^%1ixOs^S,idir) = 1.d0/11.d0*&
                   ( -2.d0*block%ws(ix1+4^%1ixOs^S,idir)&
                    +11.d0*block%ws(ix1+3^%1ixOs^S,idir)&
                    -27.d0*block%ws(ix1+2^%1ixOs^S,idir)&
                    +29.d0*block%ws(ix1+1^%1ixOs^S,idir))
                !! zero-gradient 
                !block%ws(ix1^%1ixOs^S,idir) = third*&
                !       (-block%ws(ix1+2^%1ixOs^S,idir)&
                !   +4.d0*block%ws(ix1+1^%1ixOs^S,idir))
             end do
           end if
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
         do ix1=ixOmax1,ixOmin1,-1
           w(ix1^%1ixO^S,mag(:))=third* &
                      (-w(ix1+2^%1ixO^S,mag(:)) &
                 +4.0d0*w(ix1+1^%1ixO^S,mag(:)))
         enddo
       end if

       if(mhd_energy) then
         do ix1=ixOmin1,ixOmax1
           w(ix1^%1ixO^S,rho_) = rbc(ix1)
           w(ix1^%1ixO^S,p_)   = pbc(ix1)
         enddo
       else if(mhd_adiab==0) then
         ! zero beta
         w(ixO^S,rho_)=sum(w(ixO^S,mag(:))**2,dim=ndim+1)
       else
         ! isothermal
         w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/x(ixO^S,1)))
       end if
       !w(ixO^S,mom(1))=w(ixO^S,mom(1))/w(ixO^S,rho_)
       call mhd_to_conserved(ixI^L,ixO^L,w,x)

     case(2)
       ixA^L=ixO^L;
       ixAmin1 = ixOmin1-1
       ixAmax1 = ixAmin1
       call mhd_to_primitive(ixI^L,ixA^L,w,x)
       ixAmax1 = ixOmax1
       do ix1 = ixOmin1,ixOmax1
          ! r^2*rho is continuous
          w(ix1^%1ixO^S,rho_) = w(ixOmin1-1^%1ixO^S,rho_)*(x(ixOmin1-1^%1ixO^S,1)/x(ix1^%1ixO^S,1))**2
          ! v_r no-inflow is continuous
          w(ix1^%1ixO^S,mom(1)) = w(ixOmin1-1^%1ixO^S,mom(1))
          where(w(ix1^%1ixO^S,mom(1))<0.d0)
            w(ix1^%1ixO^S,mom(1))=0.d0
          end where
          ! r*v_theta is continuous, angular momentum conservation
          w(ix1^%1ixO^S,mom(2)) = w(ixOmin1-1^%1ixO^S,mom(2))*x(ixOmin1-1^%1ixO^S,1)/x(ix1^%1ixO^S,1)
          ! r*v_phi is continuous, angular momentum conservation
          w(ix1^%1ixO^S,mom(3)) = w(ixOmin1-1^%1ixO^S,mom(3))*x(ixOmin1-1^%1ixO^S,1)/x(ix1^%1ixO^S,1)
       end do
       pth(ixOmin1-1^%1ixO^S) = w(ixOmin1-1^%1ixO^S,p_)/w(ixOmin1-1^%1ixO^S,rho_) ![T]
       do ix1=ixOmin1,ixOmax1
         ! isothermal boundary condition
         w(ix1^%1ixO^S,p_)   = pth(ixOmin1-1^%1ixO^S)*w(ix1^%1ixO^S,rho_)
       enddo

       if(stagger_grid) then
         do idir=1,nws
           if(idir==1) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           if(stretched_dim(1)) then
             do ix1=ixOsmin1,ixOsmax1
                block%ws(ix1^%1ixOs^S,idir) = &
                   ecz2(block%level)*block%ws(ix1-2^%1ixOs^S,idir)+&
                   ecz1(block%level)*block%ws(ix1-1^%1ixOs^S,idir)
             end do
           else
             do ix1=ixOsmin1,ixOsmax1
                !block%ws(ix1^%1ixOs^S,idir) = 0.0d0
                ! 3 one-sided equal gradient
                !block%ws(ix1^%1ixOs^S,idir) = 1.d0/11.d0*&
                !   ( -2.d0*block%ws(ix1-4^%1ixOs^S,idir)&
                !    +11.d0*block%ws(ix1-3^%1ixOs^S,idir)&
                !    -27.d0*block%ws(ix1-2^%1ixOs^S,idir)&
                !    +29.d0*block%ws(ix1-1^%1ixOs^S,idir))
                !! zero-gradient 
                block%ws(ix1^%1ixOs^S,idir) = third*&
                       (-block%ws(ix1-2^%1ixOs^S,idir)&
                   +4.d0*block%ws(ix1-1^%1ixOs^S,idir))
             end do
           end if
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
         do ix1=ixOmin1,ixOmax1
           w(ix1^%1ixO^S,mag(1))=third* &
                      (-w(ix1-2^%1ixO^S,mag(1)) &
                 +4.0d0*w(ix1-1^%1ixO^S,mag(1)))
         enddo
         w(ixO^S,mag(2:3))=0.d0
       end if
       call mhd_to_conserved(ixI^L,ixA^L,w,x)

     case(3)
       !! Non-flow walls
       w(ixO^S,mom(:))=0.d0
       !ixIM^L=ixO^L;
       !ixIMmin2=ixOmax2+1;ixIMmax2=ixOmax2+nghostcells;
       !call mhd_to_primitive(ixI^L,ixIM^L,w,x)
       !w(ixO^S,mom(1))= w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(1))
       !w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(2))
       !w(ixO^S,mom(3))= w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(3))

       if(stagger_grid) then
         do idir=1,nws
           if(idir==2) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix2=ixOsmax2,ixOsmin2,-1
             !! zero-gradient 
             block%ws(ix2^%2ixOs^S,idir)=third*&
                    (-block%ws(ix2+2^%2ixOs^S,idir)&
                +4.d0*block%ws(ix2+1^%2ixOs^S,idir))
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
         do ix2=ixOmax2,ixOmin2,-1
           w(ix2^%2ixO^S,mag(:))=third* &
                      (-w(ix2+2^%2ixO^S,mag(:)) &
                 +4.0d0*w(ix2+1^%2ixO^S,mag(:)))
         enddo
       end if

       if(mhd_energy) then
         ixIM^L=ixO^L;
         ixIMmin2=ixOmax2+1;ixIMmax2=ixOmax2+nghostcells;
         if(iprob==3.and.qt==time_init) then
           if(all(w(ixM^T,mom(:)).eq.zero)) then
             call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
           else
             pth(ixIM^S)=(mhd_gamma-1.d0)*w(ixIM^S,e_)
           end if
         else
           call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
         end if
         !call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
         w(ixO^S,p_)   = pth(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3)
         w(ixO^S,rho_) =   w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
       else if(mhd_adiab==0) then
         ! zero beta
         w(ixO^S,rho_)=sum(w(ixO^S,mag(:))**2,dim=ndim+1)
       else
         ! isothermal
         w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/x(ixO^S,1)))
       end if
       !ixIM^L=ixO^L;
       !ixIMmin2=ixOmin2;ixIMmax2=ixOmax2+nghostcells;
       !call mhd_to_conserved(ixI^L,ixIM^L,w,x)
       call mhd_to_conserved(ixI^L,ixO^L,w,x)

     case(4)
       !! Non-flow walls
       w(ixO^S,mom(:))=0.d0
       !ixIM^L=ixO^L;
       !ixIMmin2=ixOmin2-nghostcells;ixIMmax2=ixOmin2-1;
       !call mhd_to_primitive(ixI^L,ixIM^L,w,x)
       !w(ixO^S,mom(1))= w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(1))
       !w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(2))
       !w(ixO^S,mom(3))= w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(3))

       if(stagger_grid) then
         do idir=1,nws
           if(idir==2) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix2=ixOsmin2,ixOsmax2
             !! zero-gradient 
             block%ws(ix2^%2ixOs^S,idir) = third*&
                    (-block%ws(ix2-2^%2ixOs^S,idir)&
                +4.d0*block%ws(ix2-1^%2ixOs^S,idir))
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
         do ix2=ixOmin2,ixOmax2
           w(ix2^%2ixO^S,mag(:))=third* &
                      (-w(ix2-2^%2ixO^S,mag(:)) &
                 +4.0d0*w(ix2-1^%2ixO^S,mag(:)))
         enddo
       end if

       if(mhd_energy) then
         ixIM^L=ixO^L;
         ixIMmin2=ixOmin2-nghostcells;ixIMmax2=ixOmin2-1;
         if(iprob==3.and.qt==time_init) then
           if(all(w(ixM^T,mom(:)).eq.zero)) then
             call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
           else
             pth(ixIM^S)=(mhd_gamma-1.d0)*w(ixIM^S,e_)
           end if
         else
           call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
         end if
         !call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
         w(ixO^S,p_)   = pth(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3)
         w(ixO^S,rho_) =   w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,rho_)
       else if(mhd_adiab==0) then
         ! zero beta
         w(ixO^S,rho_)=sum(w(ixO^S,mag(:))**2,dim=ndim+1)
       else
         ! isothermal
         w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/x(ixO^S,1)))
       end if
       !ixIM^L=ixO^L;
       !ixIMmin2=ixOmin2-nghostcells;ixIMmax2=ixOmax2;
       !call mhd_to_conserved(ixI^L,ixIM^L,w,x)
       call mhd_to_conserved(ixI^L,ixO^L,w,x)

    case(5)
      !! Non-flow walls
      w(ixO^S,mom(:))=0.d0
      !ixIM^L=ixO^L;
      !ixIMmin3=ixOmax3+1;ixIMmax3=ixOmax3+nghostcells;
      !call mhd_to_primitive(ixI^L,ixIM^L,w,x)
      !w(ixO^S,mom(1))= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(1))
      !w(ixO^S,mom(2))= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(2))
      !w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(3))

      if(stagger_grid) then
        do idir=1,nws
        if(idir==3) cycle
          ixOsmax^D=ixOmax^D;
          ixOsmin^D=ixOmin^D-kr(^D,idir);
          do ix3=ixOsmax3,ixOsmin3,-1
            !! zero-gradient 
            block%ws(ix3^%3ixOs^S,idir)=third*&
                   (-block%ws(ix3+2^%3ixOs^S,idir)&
               +4.d0*block%ws(ix3+1^%3ixOs^S,idir))
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
        do ix3=ixOmax3,ixOmin3,-1
          w(ix3^%3ixO^S,mag(:))=third* &
                     (-w(ix3+2^%3ixO^S,mag(:)) &
                +4.0d0*w(ix3+1^%3ixO^S,mag(:)))
        enddo
      end if

      if(mhd_energy) then
        ixIM^L=ixO^L;
        ixIMmin3=ixOmax3+1;ixIMmax3=ixOmax3+nghostcells;
        if(iprob==3.and.qt==time_init) then
          if(all(w(ixM^T,mom(:)).eq.zero)) then
            call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
          else
            pth(ixIM^S)=(mhd_gamma-1.d0)*w(ixIM^S,e_)
          end if
        else
          call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
        end if
        !call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
        w(ixO^S,p_)   = pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1)
        w(ixO^S,rho_) =   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,rho_)
      else if(mhd_adiab==0) then
        ! zero beta
        w(ixO^S,rho_)=sum(w(ixO^S,mag(:))**2,dim=ndim+1)
      else
        ! isothermal
        w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                      (1.d0/SRadius-1.d0/x(ixO^S,1)))
      end if
      !ixIM^L=ixO^L;
      !ixIMmin3=ixOmin3;ixIMmax3=ixOmax3+nghostcells;
      !call mhd_to_conserved(ixI^L,ixIM^L,w,x)
      call mhd_to_conserved(ixI^L,ixO^L,w,x)

    case(6)
      !! Non-flow walls
      w(ixO^S,mom(:))=0.d0
      !ixIM^L=ixO^L;
      !ixIMmin3=ixOmin3-nghostcells;ixIMmax3=ixOmin3-1;
      !call mhd_to_primitive(ixI^L,ixIM^L,w,x)
      !w(ixO^S,mom(1))= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(1))
      !w(ixO^S,mom(2))= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(2))
      !w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(3))

      if(stagger_grid) then
        do idir=1,nws
          if(idir==3) cycle
          ixOsmax^D=ixOmax^D;
          ixOsmin^D=ixOmin^D-kr(^D,idir);
          do ix3=ixOsmin3,ixOsmax3
            !! zero-gradient 
            block%ws(ix3^%3ixOs^S,idir)=third*&
                   (-block%ws(ix3-2^%3ixOs^S,idir)&
               +4.d0*block%ws(ix3-1^%3ixOs^S,idir))
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
        do ix3=ixOmin3,ixOmax3
          w(ix3^%3ixO^S,mag(:))=third* &
                     (-w(ix3-2^%3ixO^S,mag(:)) &
                +4.0d0*w(ix3-1^%3ixO^S,mag(:)))
        enddo
      end if

      if(mhd_energy) then
        ixIM^L=ixO^L;
        ixIMmin3=ixOmin3-nghostcells;ixIMmax3=ixOmin3-1;
        ! internal energy to pressure
        if(iprob==3.and.qt==time_init) then
          if(all(w(ixM^T,mom(:)).eq.zero)) then
            call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
          else
            pth(ixIM^S)=(mhd_gamma-1.d0)*w(ixIM^S,e_)
          end if
        else
          call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
        end if
        !call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
        w(ixO^S,p_)   = pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1)
        w(ixO^S,rho_) =   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,rho_)
      else if(mhd_adiab==0) then
        ! zero beta
        w(ixO^S,rho_)=sum(w(ixO^S,mag(:))**2,dim=ndim+1)
      else
        ! isothermal
        w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                      (1.d0/SRadius-1.d0/x(ixO^S,1)))
      end if
      !ixIM^L=ixO^L;
      !ixIMmin3=ixOmin3-nghostcells;ixIMmax3=ixOmax3;
      !call mhd_to_conserved(ixI^L,ixIM^L,w,x)
      call mhd_to_conserved(ixI^L,ixO^L,w,x)

     case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
  ! Enforce additional refinement or coarsening
  ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen
    !if (qt>0) then
    !  refine=-1
    !  coarsen=1
    !end if
    ! Keep the bottom boundary at the highest level and allow AMR elsewhere.
    if (ps(igrid)%is_physical_boundary(1)) then
      refine=1
      coarsen=-1
    end if

  end subroutine special_refine_grid

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: w1(ixI^S,nw+nwauxio)
    double precision                   :: qvec(ixI^S,1:ndim)
    double precision                   :: tmp(ixI^S)
    double precision                   :: current(ixI^S,3)
    double precision                   :: xlen^D
    integer :: ix^D,ixbc^D,idirmin,idir,jdir,kdir
    double precision                   :: divb(ixI^S),B2(ixI^S),dip(ixI^S)
    double precision                   :: divv(ixI^S)
    double precision                   :: vels(ixI^S,3)

    w(ixO^S,nw+1:nw+nwauxio)=zero
    w(ixO^S,nw+1)=w(ixO^S,mag(1))
    w(ixO^S,nw+2)=w(ixO^S,mag(2))
    w(ixO^S,nw+3)=w(ixO^S,mag(3))
    call get_current(w,ixI^L,ixO^L,idirmin,current)
    B2(ixO^S)=sum(w(ixO^S,mag(1:ndir))**2,dim=ndim+1)
    !w(ixO^S,nw+4)=current(ixO^S,1)
    !w(ixO^S,nw+5)=current(ixO^S,2)
    !w(ixO^S,nw+6)=current(ixO^S,3)
    w1(ixI^S,1:nw+nwauxio)=zero
    w1(ixI^S,1:nw)=block%w(ixI^S,1:nw)
    call mhd_to_primitive(ixI^L,ixI^L,w1,x)
    !vels(ixI^S,1)=w(ixI^S,mom(1))/w(ixI^S,rho_)
    !vels(ixI^S,2)=w(ixI^S,mom(2))/w(ixI^S,rho_)
    !vels(ixI^S,3)=w(ixI^S,mom(3))/w(ixI^S,rho_)
    vels(ixI^S,1)=w1(ixI^S,mom(1))
    vels(ixI^S,2)=w1(ixI^S,mom(2))
    vels(ixI^S,3)=w1(ixI^S,mom(3))
    w(ixO^S,nw+4)=vels(ixO^S,1)
    w(ixO^S,nw+5)=dabs(vels(ixO^S,1)*w(ixO^S,mag(1))+&
                       vels(ixO^S,2)*w(ixO^S,mag(2))+&
                       vels(ixO^S,3)*w(ixO^S,mag(3)))/dsqrt(B2(ixO^S))
    w(ixO^S,nw+6)=dsqrt(sum(vels(ixO^S,:)**2,dim=ndim+1)-w(ixO^S,nw+5)**2)

    !w(ixO^S,nw+7:nw+10)=zero
    !if(block%is_physical_boundary(1)) then
    !  ix1=ixMlo1
    !  do ix3=ixOmin3,ixOmax3
    !    do ix2=ixOmin2,ixOmax2
    !      xlen2=x(ix^D,2)-xprobmin2
    !      xlen3=x(ix^D,3)-xprobmin3
    !      ixbc2=ceiling(xlen2/dx(2,ilevel))
    !      ixbc3=ceiling(xlen3/dx(3,ilevel))
    !      ixbc2=max(ixbc2,1)
    !      ixbc2=min(ixbc2,m)
    !      ixbc3=max(ixbc3,1)
    !      ixbc3=min(ixbc3,n)
    !      w(ix^D,nw+7)=brllt(ixbc2,ixbc3)
    !      w(ix^D,nw+8)=blatllt(ixbc2,ixbc3)
    !      w(ix^D,nw+9)=blonllt(ixbc2,ixbc3)
    !      w(ix^D,nw+10)=block%ws(nghostcells,ix2,ix3,1)
    !    end do
    !  end do
    !end if

    !! output Alfven wave speed B/sqrt(rho)
    !w(ixO^S,nw+7)=dsqrt(B2(ixO^S)/w(ixO^S,rho_))
    w(ixO^S,nw+7)=dsqrt(sum(current(ixO^S,:)**2,ndim+1))
    ! output the plasma beta p*2/B**2
    call mhd_get_pthermal(w,x,ixI^L,ixO^L,tmp)
    w(ixO^S,nw+8)=2.d0*tmp(ixO^S)/B2(ixO^S)
    ! store temperature
    call mhd_get_pthermal(w,x,ixI^L,ixO^L,tmp)
    w(ixO^S,nw+9)=tmp(ixO^S)/w(ixO^S,rho_)
    !! output divB1
    !call get_divb(w,ixI^L,ixO^L,divb)
    !w(ixO^S,nw+10)=divb(ixO^S)
    call divvector(vels(ixI^S,1:ndir),ixI^L,ixO^L,divv)
    w(ixO^S,nw+10)=divv(ixO^S)

    !! find magnetic dips                                     
    ! dip=0.d0
    ! do idir=1,ndir
    !   call gradient(w(ixI^S,mag(idir)),ixI^L,ixO^L,1,tmp)
    !   call gradientS(w(ixI^S,mag(idir)),ixI^L,ixM^LL,1,tmp)
    !   dip(ixO^S)=dip(ixO^S)+w(ixO^S,mag(idir))*tmp(ixO^S)
    ! end do
    ! B2(ixO^S)=sum((w(ixO^S,mag(:)))**2,dim=ndim+1)
    ! B2(ixO^S)=w(ixO^S,mag(1))/dsqrt(B2(ixO^S))
    ! where(dip(ixO^S)>=0.d0)
    !   w(ixO^S,nw+10)=dabs(B2(ixO^S))
    ! !where(dabs(B2(ixO^S))<0.1d0 .and. dip(ixO^S)>=0.d0)
    ! !  w(ixO^S,nw+10)=2.d0
    ! !elsewhere(dabs(B2(ixO^S))<0.01d0 .and. dip(ixO^S)>=0.d0)
    ! !  w(ixO^S,nw+10)=1.d0
    ! elsewhere
    !   w(ixO^S,nw+10)=0.d0
    ! end where

    ! JXB
    qvec(ixI^S,1:ndir)=zero
    do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
      if(lvc(idir,jdir,kdir)==1) then
        qvec(ixO^S,idir)=qvec(ixO^S,idir)+current(ixO^S,jdir)*w(ixO^S,mag(kdir))
      else if(lvc(idir,jdir,kdir)==-1) then
        qvec(ixO^S,idir)=qvec(ixO^S,idir)-current(ixO^S,jdir)*w(ixO^S,mag(kdir))
      end if
    end do; end do; end do

    ! sin_theta = |JXB|/(JB)
    w(ixO^S,nw+11)=sqrt(sum(qvec(ixO^S,:)**2,ndim+1))/&
                   sqrt(sum(w(ixO^S,mag(:))**2,ndim+1)*sum(current(ixO^S,:)**2,ndim+1))

    call get_normalized_divb(w,ixI^L,ixO^L,tmp)
    w(ixO^S,nw+12)=tmp(ixO^S)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    use mod_global_parameters

    character(len=*) varnames

!   varnames='br bt bp J1 J2 J3 Alfv beta Te divB CW fi'
!   varnames='br bt bp J1 J2 J3 Jn beta Te dip CW fi'
    varnames='br bt bp vr vpar vver Jn beta Te divv CW fi'

  end subroutine specialvarnames_output

  subroutine specialthreshold(wlocal,xlocal,tolerance,qt,level)
    !PURPOSE: use different tolerance in special regions for AMR to
    !reduce/increase resolution there where nothing/something interesting happens.
    use mod_global_parameters

    double precision, intent(in) :: wlocal(1:nw),xlocal(1:ndim),qt
    double precision, intent(inout) :: tolerance
    integer, intent(in) :: level
    
    double precision :: bczone^D,addtol,tol_add

    tol_add=0.d0
    !amplitude of additional tolerance
    addtol=0.5d0
    ! thickness of near-boundary region
    bczone1=0.3d0*(xprobmax1-xprobmin1)
    bczone2=0.15d0*(xprobmax2-xprobmin2)
    bczone3=0.15d0*(xprobmax3-xprobmin3)
    ! linear changing of additional tolerance
    if(xlocal(3)-xprobmin3 < bczone3 .or. xprobmax3-xlocal(3) < bczone3) then
      tol_add=(1.d0-min(xlocal(3)-xprobmin3,xprobmax3-xlocal(3))/bczone3)*addtol
    end if
    if(xlocal(2)-xprobmin2 < bczone2 .or. xprobmax2-xlocal(2) < bczone2) then
      tol_add=(1.d0-min(xlocal(2)-xprobmin2,xprobmax2-xlocal(2))/bczone2)*addtol
    end if
    if(xprobmax1-xlocal(1) < bczone1) then
      tol_add=(1.d0-(xprobmax1-xlocal(1))/bczone1)*addtol
    end if
    tolerance=tolerance+tol_add

  end subroutine specialthreshold

  subroutine pfss_radial_stencil(radius, radial_index, radial_weight)
    ! PDFI_SS uses uniformly spaced log(r) shells when rsun is negative.
    double precision, intent(in)  :: radius
    integer, intent(out)          :: radial_index(4)
    double precision, intent(out) :: radial_weight(4)
    integer                       :: ibase
    double precision              :: shell_index,t,t2,t3

    if(radius<=xprobmin1) then
      radial_index=(/1,1,1,1/)
      radial_weight=(/one,zero,zero,zero/)
      return
    end if

    shell_index=one+dble(p)*dlog(min(radius,Rsrs)/xprobmin1) &
                               /dlog(Rsrs/xprobmin1)
    ibase=floor(shell_index)
    if(ibase<=1) then
      radial_index=(/1,2,3,4/)
      t=shell_index-one
      radial_weight(1)=-(t-one)*(t-two)*(t-3.d0)/6.d0
      radial_weight(2)= t*(t-two)*(t-3.d0)/two
      radial_weight(3)=-t*(t-one)*(t-3.d0)/two
      radial_weight(4)= t*(t-one)*(t-two)/6.d0
    else if(ibase>=p-1) then
      radial_index=(/p-2,p-1,p,p+1/)
      t=shell_index-dble(p-2)
      radial_weight(1)=-(t-one)*(t-two)*(t-3.d0)/6.d0
      radial_weight(2)= t*(t-two)*(t-3.d0)/two
      radial_weight(3)=-t*(t-one)*(t-3.d0)/two
      radial_weight(4)= t*(t-one)*(t-two)/6.d0
    else
      radial_index=(/ibase-1,ibase,ibase+1,ibase+2/)
      t=shell_index-dble(ibase)
      t2=t*t
      t3=t2*t
      radial_weight(1)=-half*t+t2-half*t3
      radial_weight(2)= one-2.5d0*t2+1.5d0*t3
      radial_weight(3)= half*t+two*t2-1.5d0*t3
      radial_weight(4)=-half*t2+half*t3
    end if
  end subroutine pfss_radial_stencil

  double precision function pfss_avpp_at(ix2,ix3lo,ix3hi,radius) result(value)
    integer, intent(in)          :: ix2,ix3lo,ix3hi
    double precision, intent(in) :: radius
    integer                     :: ir,radial_index(4)
    double precision            :: radial_weight(4)

    call pfss_radial_stencil(radius,radial_index,radial_weight)
    value=zero
    do ir=1,4
      value=value+radial_weight(ir)*sum(Avpp(ix2,ix3lo:ix3hi,radial_index(ir)))
    end do
    value=value/dble(ix3hi-ix3lo+1)
  end function pfss_avpp_at

  double precision function pfss_avpt_at(ix2lo,ix2hi,ix3,radius) result(value)
    integer, intent(in)          :: ix2lo,ix2hi,ix3
    double precision, intent(in) :: radius
    integer                     :: ir,radial_index(4)
    double precision            :: radial_weight(4)

    call pfss_radial_stencil(radius,radial_index,radial_weight)
    value=zero
    do ir=1,4
      value=value+radial_weight(ir)*sum(Avpt(ix2lo:ix2hi,ix3,radial_index(ir)))
    end do
    value=value/dble(ix2hi-ix2lo+1)
  end function pfss_avpt_at

  subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    integer, intent(in)                :: ixI^L, ixC^L,idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)

    integer :: ix^D,ixg^D,dix
    double precision :: xlen^D
    double precision :: xCC(ixI^S,1:ndim)
    double precision :: A_vec(ixI^S,1:ndim)
    double precision :: A_rbsl(ixI^S,1:ndim)

    A=zero
    if (iprob==1) then 
      if (idir==3) then
        dix=2**(ilevel-1-block%level)
        if(dix>0) then
          {do ix^DB=ixCmin^DB,ixCmax^DB\}
            xlen^D=xC({ix^D,},^D)-xprobmin^D\
            if(.not.stretched_dim(1)) then
              ixg1=nint(xlen1/dx(1,ilevel))+1
            end if
            ixg2=nint(xlen2/dx(2,ilevel))+1
            ixg3=nint(xlen3/dx(3,ilevel))
            if(xC(ix^D,1)>Rsrs) then
              A(ix^D)=sum(Avpp(ixg2,ixg3-dix+1:ixg3+dix,p+1))/dble(2*dix)*Rsrs/xC(ix^D,1)
            else if(stretched_dim(1)) then
              A(ix^D)=pfss_avpp_at(ixg2,ixg3-dix+1,ixg3+dix,xC(ix^D,1))
            else
              A(ix^D)=sum(Avpp(ixg2,ixg3-dix+1:ixg3+dix,ixg1))/dble(2*dix)
            end if
          {end do\}
        else
          {do ix^DB=ixCmin^DB,ixCmax^DB\}
            xlen^D=xC({ix^D,},^D)-xprobmin^D\
            if(.not.stretched_dim(1)) then
              ixg1=nint(xlen1/dx(1,ilevel))+1
            end if
            ixg2=nint(xlen2/dx(2,ilevel))+1
            ixg3=ceiling(xlen3/dx(3,ilevel))
            if(xC(ix^D,1)>Rsrs) then
              A(ix^D)=Avpp(ixg2,ixg3,p+1)*Rsrs/xC(ix^D,1)
            else if(stretched_dim(1)) then
              A(ix^D)=pfss_avpp_at(ixg2,ixg3,ixg3,xC(ix^D,1))
            else
              A(ix^D)=Avpp(ixg2,ixg3,ixg1)
            end if
          {end do\}
        end if
      else if(idir==2) then
        dix=2**(ilevel-1-block%level)
        if(dix>0) then
          {do ix^DB=ixCmin^DB,ixCmax^DB\}
            xlen^D=xC({ix^D,},^D)-xprobmin^D\
            if(.not.stretched_dim(1)) then
              ixg1=nint(xlen1/dx(1,ilevel))+1
            end if
            ixg3=nint(xlen3/dx(3,ilevel))+1
            ixg2=nint(xlen2/dx(2,ilevel))
            if(xC(ix^D,1)>Rsrs) then
              A(ix^D)=sum(Avpt(ixg2-dix+1:ixg2+dix,ixg3,p+1))/dble(2*dix)*Rsrs/xC(ix^D,1)
            else if(stretched_dim(1)) then
              A(ix^D)=pfss_avpt_at(ixg2-dix+1,ixg2+dix,ixg3,xC(ix^D,1))
            else
              A(ix^D)=sum(Avpt(ixg2-dix+1:ixg2+dix,ixg3,ixg1))/dble(2*dix)
            end if
          {end do\}
        else
          {do ix^DB=ixCmin^DB,ixCmax^DB\}
            xlen^D=xC({ix^D,},^D)-xprobmin^D\
            if(.not.stretched_dim(1)) then
              ixg1=nint(xlen1/dx(1,ilevel))+1
            end if
            ixg3=nint(xlen3/dx(3,ilevel))+1
            ixg2=ceiling(xlen2/dx(2,ilevel))
            if(xC(ix^D,1)>Rsrs) then
              A(ix^D)=Avpt(ixg2,ixg3,p+1)*Rsrs/xC(ix^D,1)
            else if(stretched_dim(1)) then
              A(ix^D)=pfss_avpt_at(ixg2,ixg2,ixg3,xC(ix^D,1))
            else
              A(ix^D)=Avpt(ixg2,ixg3,ixg1)
            end if
          {end do\}
        end if
      else
        A(ixC^S)=0.d0
      end if
    end if

    if (iprob==2) then 
      xCC(ixI^S,1) = xC(ixI^S,1) * cos(half*dpi-xC(ixI^S,2)) * cos(xC(ixI^S,3))
      xCC(ixI^S,2) = xC(ixI^S,1) * cos(half*dpi-xC(ixI^S,2)) * sin(xC(ixI^S,3))
      xCC(ixI^S,3) = xC(ixI^S,1) * sin(half*dpi-xC(ixI^S,2))
      call RBSL_flux_rope(ixI^L,ixC^L,nnp,a0,F_flx,I_Helix,xCC,x_axis,A_vec)
      call Cart2SphereVector(ixI^L,ixC^L,xC,A_vec,A_rbsl)
      if (idir==3) then
        A(ixC^S)=A(ixC^S)+A_rbsl(ixC^S,3)
      else if(idir==2) then
        A(ixC^S)=A(ixC^S)+A_rbsl(ixC^S,2)
      else
        A(ixC^S)=A(ixC^S)+A_rbsl(ixC^S,1)
      end if
    end if

    if(weak_field) A(ixC^S)=A(ixC^S)*wf_scale ! weak-field process

  end subroutine initvecpot_usr

  !> This subroutine is designed for users to customize electric field.
  !> The custom-made electric field on photosphere (Eth, Eph and Er) 
  !> Note that any boundary except the photospheric boundary will be set
  !> as perfect conductor walls in order to ensure the normal magnetic
  !> field invariant at these boundaries.
  subroutine boundary_electric_field(ixI^L,ixO^L,qt,qdt,fE,s)
    ! specify tangential electric field at physical boundaries 
    ! to fix or drive normal magnetic field
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt,qdt
    type(state)                        :: s
    double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)

    double precision :: xC(ixI^S,1:ndim),xlen^D
    integer :: idir,ixC^L,ixA^L,ix^D,ixg^D,dix

    associate(w=>s%w,x=>s%x,ws=>s%ws)

    if(s%is_physical_boundary(1)) then
      ixCmin^D=ixOmin^D-kr(2,^D);
      ixCmax^D=ixOmax^D;
      fE(nghostcells^%1ixC^S,3)=0.d0
      ixCmin^D=ixOmin^D-kr(3,^D);
      ixCmax^D=ixOmax^D;
      fE(nghostcells^%1ixC^S,2)=0.d0
    end if

    ! perfect conductor theta boundaries
    !if(block%is_physical_boundary(3)) then
    !  ixCmin^D=ixOmin^D-kr(2,^D);
    !  ixCmax^D=ixOmax^D;
    !  ixAmin^D=ixCmin^D-kr(3,^D);
    !  ixAmax^D=ixCmax^D;
    !  fE(ixAmin2^%2ixA^S,1)=0.d0
    !  ixAmin^D=ixCmin^D-kr(1,^D);
    !  ixAmax^D=ixCmax^D;
    !  fE(ixAmin2^%2ixA^S,3)=0.d0
    !end if
    !if(block%is_physical_boundary(4)) then
    !  ixCmin^D=ixOmin^D-kr(2,^D);
    !  ixCmax^D=ixOmax^D;
    !  ixAmin^D=ixCmin^D-kr(3,^D);
    !  ixAmax^D=ixCmax^D;
    !  fE(ixAmax2^%2ixA^S,1)=0.d0
    !  ixAmin^D=ixCmin^D-kr(1,^D);
    !  ixAmax^D=ixCmax^D;
    !  fE(ixAmax2^%2ixA^S,3)=0.d0
    !end if
    !! perfect conductor phi boundaries
    !if(block%is_physical_boundary(5)) then
    !  ixCmin^D=ixOmin^D-kr(3,^D);
    !  ixCmax^D=ixOmax^D;
    !  ixAmin^D=ixCmin^D-kr(2,^D);
    !  ixAmax^D=ixCmax^D;
    !  fE(ixAmin3^%3ixA^S,1)=0.d0
    !  ixAmin^D=ixCmin^D-kr(1,^D);
    !  ixAmax^D=ixCmax^D;
    !  fE(ixAmin3^%3ixA^S,2)=0.d0
    !end if
    !if(block%is_physical_boundary(6)) then
    !  ixCmin^D=ixOmin^D-kr(3,^D);
    !  ixCmax^D=ixOmax^D;
    !  ixAmin^D=ixCmin^D-kr(2,^D);
    !  ixAmax^D=ixCmax^D;
    !  fE(ixAmax3^%3ixA^S,1)=0.d0
    !  ixAmin^D=ixCmin^D-kr(1,^D);
    !  ixAmax^D=ixCmax^D;
    !  fE(ixAmax3^%3ixA^S,2)=0.d0
    !end if

    end associate

  end subroutine boundary_electric_field

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
  ! Purpose: background heating and localized heating
  !==============================================================================
 subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: lQgrid(ixI^S),bQgrid(ixI^S)

    !! add global background heating bQ
    call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)
    if(liakh_mass_loading) then
      call add_liakh_filament_mass_source(qdt,ixI^L,ixO^L,qtC,w,x)
    end if
    !! add localized heating lQ
    !if(iprob==2) then
    !  call getlQ(lQgrid,ixI^L,ixO^L,qtC,wCT,x)
    !  w(ixO^S,e_)=w(ixO^S,e_)+qdt*lQgrid(ixO^S)
    !end if

  end subroutine specialsource

  subroutine my_special_resistivity(w,ixI^L,ixO^L,idirmin,x,current,eta)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L, idirmin
    double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision             :: current(ixI^S,7-2*ndir:3), eta(ixI^S)
    double precision             :: j_mag(ixI^S), B_mag(ixI^S),  zeta(ixI^S), eta_0 
    double precision             :: eta_un, eta_exp(ixI^S), eta_fp(ixI^S), delta_eta(ixI^S) 
    integer                      :: ix^D

    eta_0 = 200  !km^2/s
    eta_un = eta_0*1.0d10/(unit_length*unit_velocity) !normalization
    j_mag = dsqrt(sum(current(ixI^S,:)**2, ndim+1)) 
    B_mag = dsqrt(sum(w(ixI^S,mag(1):mag(3))**2, ndim+1))

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if(B_mag(ix^D) .ne. zero) then 
          delta_eta(ix^D) = minval(block%ds(ix^D,1:3))
          zeta(ix^D) = (j_mag(ix^D)*delta_eta(ix^D)/B_mag(ix^D))**2
          eta_fp(ix^D) = eta_un*10.0d0*zeta(ix^D)/&
                         (1+dexp(-(zeta(ix^D)-0.1)/0.01))
        else
          eta_fp(ix^D) = zero
        end if
        eta(ix^D) = eta_un + eta_fp(ix^D)
    {end do\}

  end subroutine my_special_resistivity


 ! =================== USER BUILT-IN SUBROUTINES ===================
 ! Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /Mag_list/ do_rbsl, weak_field, wf_scale, Bflux, I_Helix, Rmfr, npoints, flip_axis, &
                         load_filament_density, filament_density_factor, &
                         filament_core_fraction, filament_min_height_Mm, &
                         liakh_mass_loading, liakh_density_contrast, &
                         liakh_t_start_s, liakh_t_load_s, &
                         liakh_sigma_horizontal_Mm, liakh_sigma_radial_Mm, &
                         liakh_min_height_Mm, &
                         background_heating_cgs

    do_rbsl=.false.
    weak_field=.true.
    wf_scale=1.d-2
    Bflux=1.d20
    I_Helix=.true.
    Rmfr=4.5d1
    flip_axis=.true.
    load_filament_density=.false.
    filament_density_factor=20.d0
    filament_core_fraction=0.6d0
    filament_min_height_Mm=5.d0
    liakh_mass_loading=.false.
    liakh_density_contrast=30.d0
    liakh_t_start_s=0.d0
    liakh_t_load_s=100.d0
    liakh_sigma_horizontal_Mm=3.6d0
    liakh_sigma_radial_Mm=4.8d0
    liakh_min_height_Mm=5.d0
    background_heating_cgs=0.d0

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, Mag_list, end=111)
111    close(unitpar)
    end do

    !> For positive helicity MFR, I_cur > 0; For negative helicity MFR, I_cur < 0
    F_flx= Bflux/(unit_magneticfield*unit_length**2)
    a0   = Rmfr*1.d8/unit_length
    np   = npoints

    call read_atmosphere_params(files)

  end subroutine usr_params_read

  subroutine read_atmosphere_params(files)
    character(len=*), intent(in) :: files(:)
    integer :: n, atmosphere_unit, ios, ich, code
    character(len=2048) :: line
    character(len=512) :: io_message

    namelist /atmosphere_list/ atmosphere_model, atmosphere_normalization, &
      atmosphere_reference_height_Mm, atmosphere_reference_nH_cgs, &
      atmosphere_reference_pressure_cgs, atmosphere_chromosphere_K, atmosphere_tanh_top_K, &
      atmosphere_transition_height_Mm, atmosphere_transition_width_Mm, atmosphere_transition_K, &
      atmosphere_conductive_flux_cgs, atmosphere_floor_height_Mm, atmosphere_isothermal_K

    atmosphere_model='analytic_tr'
    atmosphere_normalization='number_density'
    atmosphere_reference_height_Mm=10.d0
    atmosphere_reference_nH_cgs=6.2d8
    atmosphere_reference_pressure_cgs=1.38365d2
    atmosphere_chromosphere_K=8.d3
    atmosphere_tanh_top_K=1.5d6
    atmosphere_transition_height_Mm=2.d0
    atmosphere_transition_width_Mm=0.25d0
    atmosphere_transition_K=1.6d5
    atmosphere_conductive_flux_cgs=2.d5
    atmosphere_floor_height_Mm=-0.2d0
    atmosphere_isothermal_K=1.d6

    ! Locate the group before reading it: EOF without a group is allowed,
    ! whereas EOF inside an explicitly supplied (truncated) group is an error.
    ! Later parameter files override only the values they supply.
    do n=1,size(files)
      open(newunit=atmosphere_unit,file=trim(files(n)),status='old',action='read', &
           iostat=ios,iomsg=io_message)
      if(ios/=0) then
        if(mype==0) write(*,*) 'Atmosphere parameter file: ',trim(files(n)),': ',trim(io_message)
        call mpistop('Cannot open atmosphere parameter file')
      end if
      do
        read(atmosphere_unit,'(a)',iostat=ios,iomsg=io_message) line
        if(ios<0) exit
        if(ios>0) then
          close(atmosphere_unit)
          if(mype==0) write(*,*) 'Atmosphere parameter file: ',trim(files(n)),': ',trim(io_message)
          call mpistop('Cannot scan atmosphere parameter file')
        end if
        do ich=1,len_trim(line)
          if(line(ich:ich)==achar(9)) line(ich:ich)=' '
          code=iachar(line(ich:ich))
          if(code>=iachar('A').and.code<=iachar('Z')) line(ich:ich)=achar(code+32)
        end do
        line=adjustl(line)
        if(line(1:1)/='&'.and.line(1:1)/='$') cycle
        line=adjustl(line(2:))
        if(line(1:15)/='atmosphere_list') cycle
        if(index(' /!',line(16:16))==0) cycle
        backspace(atmosphere_unit,iostat=ios)
        if(ios/=0) then
          close(atmosphere_unit)
          call mpistop('Cannot seek atmosphere namelist')
        end if
        read(atmosphere_unit,nml=atmosphere_list,iostat=ios,iomsg=io_message)
        if(ios/=0) then
          close(atmosphere_unit)
          if(mype==0) write(*,*) 'Invalid atmosphere_list in ',trim(files(n)),': ',trim(io_message)
          call mpistop('Invalid or incomplete atmosphere_list')
        end if
      end do
      close(atmosphere_unit)
    end do

    call validate_atmosphere_parameters()
  end subroutine read_atmosphere_params

  subroutine validate_atmosphere_parameters()
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    integer :: ich, code

    ! Model names are case-insensitive; no physics switches are changed here.
    atmosphere_model=adjustl(atmosphere_model)
    atmosphere_normalization=adjustl(atmosphere_normalization)
    do ich=1,len(atmosphere_model)
      code=iachar(atmosphere_model(ich:ich))
      if(code>=iachar('A').and.code<=iachar('Z')) atmosphere_model(ich:ich)=achar(code+32)
      code=iachar(atmosphere_normalization(ich:ich))
      if(code>=iachar('A').and.code<=iachar('Z')) atmosphere_normalization(ich:ich)=achar(code+32)
    end do
    if(atmosphere_model/='analytic_tr'.and.atmosphere_model/='isothermal') &
      call mpistop('atmosphere_model must be analytic_tr or isothermal')
    if(atmosphere_normalization/='number_density'.and.atmosphere_normalization/='pressure') &
      call mpistop('atmosphere_normalization must be number_density or pressure')
    if(.not.all(ieee_is_finite((/atmosphere_reference_height_Mm,atmosphere_reference_nH_cgs, &
        atmosphere_reference_pressure_cgs,atmosphere_chromosphere_K,atmosphere_tanh_top_K, &
        atmosphere_transition_height_Mm,atmosphere_transition_width_Mm,atmosphere_transition_K, &
        atmosphere_conductive_flux_cgs,atmosphere_floor_height_Mm,atmosphere_isothermal_K/)))) &
      call mpistop('Atmosphere parameters must all be finite')
    if(atmosphere_reference_nH_cgs<=0.d0.or.atmosphere_reference_pressure_cgs<=0.d0) &
      call mpistop('Atmosphere reference density and pressure must be positive')
    if(atmosphere_chromosphere_K<=0.d0.or.atmosphere_tanh_top_K<=0.d0.or. &
       atmosphere_transition_K<=0.d0.or.atmosphere_isothermal_K<=0.d0) &
      call mpistop('Atmosphere temperatures must be positive')
    if(atmosphere_transition_width_Mm<=0.d0) &
      call mpistop('Atmosphere transition width must be positive')
    if(atmosphere_conductive_flux_cgs<0.d0) &
      call mpistop('Atmosphere conductive flux must be nonnegative')
    if(atmosphere_model=='analytic_tr') then
      if(atmosphere_chromosphere_K>=atmosphere_transition_K.or. &
         atmosphere_transition_K>=atmosphere_tanh_top_K) &
        call mpistop('analytic_tr requires T_chromosphere < T_transition < T_tanh_top')
      if(atmosphere_floor_height_Mm>=atmosphere_transition_height_Mm) &
        call mpistop('Atmosphere floor height must be below the transition height')
    end if
  end subroutine validate_atmosphere_parameters

  subroutine enhance_filament_density(ixI^L,ixO^L,w,x)
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer          :: ix^D,iax,axis_stride
    double precision :: xcart(1:ndim),d2,d2min,xi,weight,rho_factor
    double precision :: height_Mm

    ! About 200 samples resolve the smooth observed axis to much better than a0.
    axis_stride=max(1,np/200)
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      height_Mm=(x(ix^D,1)-xprobmin1)*unit_length/1.d8
      if(height_Mm < filament_min_height_Mm) cycle

      xcart(1)=x(ix^D,1)*sin(x(ix^D,2))*cos(x(ix^D,3))
      xcart(2)=x(ix^D,1)*sin(x(ix^D,2))*sin(x(ix^D,3))
      xcart(3)=x(ix^D,1)*cos(x(ix^D,2))
      d2min=huge(one)
      do iax=1,np,axis_stride
        d2=sum((xcart(:)-x_axis(iax,:))**2)
        d2min=min(d2min,d2)
      end do
      d2=sum((xcart(:)-x_axis(np,:))**2)
      d2min=min(d2min,d2)
      xi=sqrt(d2min)/a0

      if(xi < one) then
        if(xi <= filament_core_fraction) then
          weight=one
        else
          weight=half*(one+cos(dpi*(xi-filament_core_fraction)/ &
                 (one-filament_core_fraction)))
        end if
        rho_factor=one+(filament_density_factor-one)*weight
        w(ix^D,rho_)=w(ix^D,rho_)*rho_factor
      end if
    {end do\}
  end subroutine enhance_filament_density

  ! Precompute the nearest observed RBSL-axis direction on the finest angular
  ! mesh. This makes the time-dependent source inexpensive during each RK stage.
  subroutine prepare_liakh_mass_source_geometry()
    integer          :: ith,iph,iax,axis_stride,best_axis
    double precision :: theta,phi,axis_theta,axis_phi,axis_r
    double precision :: theta_min,theta_max,phi_min,phi_max,margin
    double precision :: ray(3),axis_unit(3),dot_product,angle,best_angle

    if(liakh_t_load_s<=zero) call mpistop('Liakh t_load must be positive')
    if(liakh_sigma_horizontal_Mm<=zero.or.liakh_sigma_radial_Mm<=zero) &
      call mpistop('Liakh source widths must be positive')

    allocate(liakh_axis_radius(m,n),liakh_axis_angle(m,n))
    liakh_axis_radius=-one
    liakh_axis_angle=huge(one)

    theta_min=huge(one)
    theta_max=-huge(one)
    phi_min=huge(one)
    phi_max=-huge(one)
    do iax=1,np
      axis_r=sqrt(sum(x_axis(iax,:)**2))
      axis_theta=acos(max(-one,min(one,x_axis(iax,3)/axis_r)))
      axis_phi=atan2(x_axis(iax,2),x_axis(iax,1))
      if(axis_phi<zero) axis_phi=axis_phi+two*dpi
      theta_min=min(theta_min,axis_theta)
      theta_max=max(theta_max,axis_theta)
      phi_min=min(phi_min,axis_phi)
      phi_max=max(phi_max,axis_phi)
    end do
    margin=4.d0*liakh_sigma_horizontal_Mm*1.d8/unit_length
    axis_stride=max(1,np/600)

    do ith=1,m
      theta=xprobmin2+(dble(ith)-half)*dtheta
      if(theta<theta_min-margin.or.theta>theta_max+margin) cycle
      do iph=1,n
        phi=xprobmin3+(dble(iph)-half)*dphi
        if(phi<phi_min-margin.or.phi>phi_max+margin) cycle
        ray(1)=sin(theta)*cos(phi)
        ray(2)=sin(theta)*sin(phi)
        ray(3)=cos(theta)
        best_angle=huge(one)
        best_axis=1
        do iax=1,np,axis_stride
          axis_r=sqrt(sum(x_axis(iax,:)**2))
          axis_unit=x_axis(iax,:)/axis_r
          dot_product=max(-one,min(one,sum(ray*axis_unit)))
          angle=acos(dot_product)
          if(angle<best_angle) then
            best_angle=angle
            best_axis=iax
          end if
        end do
        axis_r=sqrt(sum(x_axis(np,:)**2))
        axis_unit=x_axis(np,:)/axis_r
        dot_product=max(-one,min(one,sum(ray*axis_unit)))
        angle=acos(dot_product)
        if(angle<best_angle) then
          best_angle=angle
          best_axis=np
        end if
        liakh_axis_radius(ith,iph)=sqrt(sum(x_axis(best_axis,:)**2))
        liakh_axis_angle(ith,iph)=best_angle
      end do
    end do

    if(mype==0) then
      write(*,*) 'Liakh mass source: chi, start/load [s] =', &
           liakh_density_contrast,liakh_t_start_s,liakh_t_load_s
      write(*,*) 'Liakh mass source: sigma_h/sigma_r [Mm] =', &
           liakh_sigma_horizontal_Mm,liakh_sigma_radial_Mm
    end if
  end subroutine prepare_liakh_mass_source_geometry

  ! Liakh et al. (2023), Eq. (8), generalized from a 2.5D straight tube to
  ! a curved 3D RBSL axis. Only density is added; momentum and energy are kept.
  subroutine add_liakh_filament_mass_source(qdt,ixI^L,ixO^L,qtC,w,x)
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qdt,qtC,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer          :: ix^D,ith,iph,na
    double precision :: t_start,t_end,t_load,height_Mm
    double precision :: axis_r,angle,d_horizontal_Mm,d_radial_Mm
    double precision :: exponent,weight,res,rho_background,rho_source

    t_start=liakh_t_start_s/unit_time
    t_load=liakh_t_load_s/unit_time
    t_end=t_start+t_load
    if(qtC<t_start.or.qtC>=t_end) return

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      height_Mm=(x(ix^D,1)-xprobmin1)*unit_length/1.d8
      if(height_Mm<liakh_min_height_Mm) cycle

      ith=max(1,min(m,int((x(ix^D,2)-xprobmin2)/dtheta)+1))
      iph=max(1,min(n,int((x(ix^D,3)-xprobmin3)/dphi)+1))
      axis_r=liakh_axis_radius(ith,iph)
      if(axis_r<zero) cycle
      angle=liakh_axis_angle(ith,iph)
      d_horizontal_Mm=angle*half*(x(ix^D,1)+axis_r)*unit_length/1.d8
      d_radial_Mm=(x(ix^D,1)-axis_r)*unit_length/1.d8
      exponent=(d_horizontal_Mm/liakh_sigma_horizontal_Mm)**4+ &
               (d_radial_Mm/liakh_sigma_radial_Mm)**4
      if(exponent>50.d0) cycle
      weight=exp(-exponent)

      na=floor((x(ix^D,1)-xprobmin1+gzone)/dr+half)
      na=max(1,min(jmax-1,na))
      res=x(ix^D,1)-xprobmin1+gzone-(dble(na)-half)*dr
      rho_background=ra(na)+(one-cos(dpi*res/dr))*half*(ra(na+1)-ra(na))
      rho_source=liakh_density_contrast*rho_background/t_load*weight
      w(ix^D,rho_)=w(ix^D,rho_)+qdt*rho_source
    {end do\}
  end subroutine add_liakh_filament_mass_source

  subroutine extrapolate_potential_field()
    integer          :: bcm,bcn
    double precision :: a,b,c,d,rsun,rssd
    double precision :: mflux
    double precision, allocatable :: scrb3d(:,:,:)
    character*80     :: brfile


    if(mype==0) then
      write(*,*) 'npe =',npe
      if(stretched_dim(1)) then
        rsun=-6.955d5 ! km
      else
        rsun=6.955d5 ! km
      end if
      rssd=ssd*abs(rsun)
      a=xprobmin2
      b=xprobmax2
      if(xprobmin3<0)then
        c=xprobmin3-xprobmin3
        d=xprobmax3-xprobmin3
      else if(xprobmax3>2.0d0*dpi) then
             c=xprobmin3-xprobmax3+2.0d0*dpi
             d=xprobmax3-xprobmax3+2.0d0*dpi
           else
             c=xprobmin3
             d=xprobmax3
      end if

      write(*,'(a)')      '========= PDFI_ss domain Info ========='
      write(*,'(a,f6.2)') 'PDFI_ss domain Info: a = ', a
      write(*,'(a,f6.2)') 'PDFI_ss domain Info: b = ', b
      write(*,'(a,f6.2)') 'PDFI_ss domain Info: c = ', c
      write(*,'(a,f6.2)') 'PDFI_ss domain Info: d = ', d
      write(*,'(a)')      '======================================='
      write(*,*) ' '

      write(*,'(a)') '============ Prepare Br fields at the initial time ============'
      write(*,'(a)') 'pf extrapolation: br read in from infile'
      write(*,'(a)') '======================== Reading Bfile ========================'
      call read_b(0,m,n,brll0)
      ! PDFI_SS needs magnetic field in Gauss
      brll0=brll0*unit_magneticfield
      write(*,'(a,f8.2,a,f8.2,a)') 'Input Brll range : ', &
                  minval(brll0), ' Gauss --> ', &
                  maxval(brll0), ' Gauss '
!     call magnetic_filter(brll0,m,n)
      if(do_rbsl) then
        write(*,'(a)') '====== brphot caculation in Regularized Biot_Savart Laws ======'
        call prep_rbsl(ixgI^L,ixgO^L,br_rbsl)
        write(*,'(a,f8.2,a,f8.2,a)') 'Br_rbsl range : ', &
                    minval(br_rbsl), ' Gauss --> ', &
                    maxval(br_rbsl), ' Gauss '
        brll0=brll0-br_rbsl
        write(*,'(a,f8.2,a,f8.2,a)') 'Final Brll range : ', &
                    minval(brll0), ' Gauss --> ', &
                    maxval(brll0), ' Gauss '
        write(*,'(a)') '====== Regularized Biot_Savart Laws caculation completed ======'
        write(*,*) ' '
      end if
      write(*,'(a)') '==============================================================='
      write(*,*) ' '

      allocate(scrb3d(m,n,p+1))

      write(*,'(a)') '========== Ahpot_ss pfss extrapolation ========='
      ! To use periodic BC in phi, set bcn=0
      !For magnetic field that does not penetrate sidewalls of spherical  wedge domain, set bcm=3 and bcn=3
      bcm=3
      bcn=3
      write(*,'(a,2I4)') 'Boundary condition set: bcm, bcn = ',bcm,bcn
      write(*,'(a)') 'Next to call scrbpot_ss ...'
      call scrbpot_ss(m,n,p,bcm,bcn,a,b,c,d,rsun,rssd,brll0,scrb3d)
      write(*,'(a)') 'After calling scrbpot_ss, there comes mflux_ss calling ...'
      call mflux_ss(m,n,a,b,c,d,rsun,brll0,mflux)
      write(*,'(a,ES16.5)') 'Just after mflux_ss call, mflux = ',mflux
      call ahpot_ss(m,n,p,bcm,bcn,a,b,c,d,rsun,rssd,scrb3d,mflux,Avpt,Avpp)

      ! nondimensionalize vector potential, Avpt and Avpp has unit of G km
      Avpt=Avpt*1.d5/(unit_length*unit_magneticfield) 
      Avpp=Avpp*1.d5/(unit_length*unit_magneticfield)

      write(*,'(a)') '========== ahpot_ss call completed ========='
      write(*,*) ' '
    end if

    if(npe>1) then
      call MPI_BCAST(Avpt,m*(n+1)*(p+1),MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(Avpp,n*(m+1)*(p+1),MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(brll0,m*n,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    end if

  end subroutine extrapolate_potential_field

  !==============================================================================
  ! Purpose: to read magnetic field from input fits files 
  !==============================================================================
  subroutine read_b(tc,mi,ni,b0rll,t_rec)
    use mod_global_parameters

    integer, intent(in) :: tc
    integer, intent(in) :: mi,ni
    double precision, intent(out) :: b0rll(mi,ni)
    double precision, intent(out), optional :: t_rec

    real              :: nullval=0.0
    real, allocatable :: btmp(:,:)
    integer           :: t_out,naxes(2),naxis,p0
    integer           :: status,readwrite,blocksize
    logical           :: anynull
    character*80      :: comment
    character*200     :: brll_magnetogram,blon_magnetogram,blat_magnetogram
   
    write (brll_magnetogram, '(a,i4.4,a)') 'Bphot/Brll',tc,'.fits' 
   
    b0rll = 0.d0

    ! read in photospheric array of Brll at cell centers (colat,lon) order.
    status = 0
    blocksize = 1
    call ftopen(126,brll_magnetogram,readwrite,blocksize,status)
    call ftgknj(126,'NAXIS',1,2,naxes,naxis,status)
    p0=naxes(1)*naxes(2)
    allocate(btmp(naxes(1),naxes(2)))
    if(present(t_rec)) then
      call ftgkyj(126, 'TIME', t_out, comment, status)
      t_rec = real(t_out,8)/unit_time
    end if

    call ftgpve(126,1,1,p0,nullval,btmp,anynull,status)
    call ftclos(126, status)
    call ftfiou(126, status)
    b0rll=dble(btmp)/unit_magneticfield
    deallocate(btmp)

  end subroutine read_b

  subroutine magnetic_filter(br_phot,ntheta,nphi)
    integer, intent(in)            :: ntheta,nphi
    double precision,intent(inout) :: br_phot(ntheta,nphi) 

    double precision :: br_max
    double precision :: br_thr=200.d0,br_top=1000.d0
    br_max = maxval(dabs(br_phot))
    where(dabs(br_phot)>br_thr)
      br_phot=br_phot/dabs(br_phot)*((dabs(br_phot)-br_thr)/(br_max-br_thr)*&
             (br_top-br_thr)+br_thr)
    end where
          
  end subroutine magnetic_filter

  subroutine flux_rope_insertion(x_axis)

  double precision, intent(out) :: x_axis(nnp,ndim)

  integer          :: iax, ixp, path_unit, path_iostat
  integer(kind=8)  :: path_bytes, expected_path_bytes
  logical          :: path_exists
  double precision :: lonax(np),latax(np),rax(np)
  double precision ::  lonm(np), latm(np), rm(np)

  double precision :: lon0,lat0,lon1,lat1
  double precision :: lonc,latc,rc
  double precision :: x0,y0,z0
  double precision :: x1,y1,z1
  double precision :: xc,yc,zc

  double precision ::  v1(nnp),  v2(nnp),  v3(nnp)
  double precision ::  xs(nnp),  ys(nnp),  zs(nnp)

    inquire(file='path.dat',exist=path_exists,size=path_bytes)
    if(.not.path_exists) call mpistop('Cannot find binary RBSL path.dat')
    expected_path_bytes=24_8*int(np,kind=8)
    if(path_bytes/=expected_path_bytes) then
      if(mype==0) then
        write(*,'(a,i0)') 'Binary path.dat bytes found: ',path_bytes
        write(*,'(a,i0)') 'Binary path.dat bytes expected: ',expected_path_bytes
      end if
      call mpistop('path.dat must contain npoints interleaved float64 triples')
    end if

    open(newunit=path_unit,file='path.dat',action='read',status='old',&
         access='stream',form='unformatted',convert='little_endian')
    do iax=1,np
      read(path_unit,iostat=path_iostat) lonax(iax),latax(iax),rax(iax)
      if(path_iostat/=0) call mpistop('Failed while reading binary path.dat')
    end do
    close(path_unit)

  if(mype==0) then
      write(6,'(a)') ' ===== The spherical coordinates of the filament from SCC Measure ====='
      write(6,'(a)') '            Longitude(deg) , Latitude(deg) , Radius( x Rsun)            '
    do iax=1,np
      write(6,'(a,3f16.5,a)') '       ',lonax(iax),latax(iax),rax(iax),'          '
    end do
      write(6,'(a)') ' ======================================================================'
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

    if(flip_axis) then
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
        xs(ixp)    = v1(ixp)    * cos(v2(ixp))    * cos(v3(ixp))
        ys(ixp)    = v1(ixp)    * cos(v2(ixp))    * sin(v3(ixp))
        zs(ixp)    = v1(ixp)    * sin(v2(ixp))
        xs(np+ixp) = v1(np+ixp) * cos(v2(np+ixp)) * cos(v3(np+ixp))
        ys(np+ixp) = v1(np+ixp) * cos(v2(np+ixp)) * sin(v3(np+ixp))
        zs(np+ixp) = v1(np+ixp) * sin(v2(np+ixp))
      end do
    end if

    x_axis(:,1) = xs
    x_axis(:,2) = ys
    x_axis(:,3) = zs

  end subroutine flux_rope_insertion

  subroutine get_lonlat(x_in,y_in,z_in,lon,lat)

  real*8 :: x_in,y_in,z_in
  real*8 :: x,y,z
  real*8 :: lon,lat
  real*8 :: Ebm

    Ebm=dsqrt(x_in**2+y_in**2+z_in**2)
    x=x_in/Ebm
    y=y_in/Ebm
    z=z_in/Ebm

    lat=asin(z)
    lon=atan2(y,x)

    return

  end subroutine get_lonlat

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

  subroutine prep_rbsl(ixI^L,ixO^L,br_mfr)

    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(out) :: br_mfr(m,n)

    integer :: ix^D
    double precision :: dx01,dx02,dx03,r_sun
    double precision :: sinth(m+1),sinth_hlf(m)
    double precision :: xgl(ixI^S,1:ndim)
    double precision :: xthC(ixI^S,1:ndim) , xphC(ixI^S,1:ndim)
    double precision :: xthCC(ixI^S,1:ndim), xphCC(ixI^S,1:ndim)
    double precision :: AthC(ixI^S,1:ndim),AthS(ixI^S,1:ndim)
    double precision :: AphC(ixI^S,1:ndim),AphS(ixI^S,1:ndim)
    double precision :: Ath_rbsl(m,n+1), Aph_rbsl(m+1,n)

    dx01=(xprobmax1-xprobmin1)/(domain_nx1*2**(refine_max_level-1))
    dx02=(xprobmax2-xprobmin2)/(domain_nx2*2**(refine_max_level-1))
    dx03=(xprobmax3-xprobmin3)/(domain_nx3*2**(refine_max_level-1))
    {do ix^DB=ixgImin^DB,ixgImax^DB\}
      xgl(ix^D,1)=(ix1-1)*dx01+xprobmin1
      xgl(ix^D,2)=(ix2-1)*dx02+xprobmin2
      xgl(ix^D,3)=(ix3-1)*dx03+xprobmin3
      xthC(ix^D,1)=xgl(ix^D,1)
      xthC(ix^D,2)=xgl(ix^D,2)+half*dx02
      xthC(ix^D,3)=xgl(ix^D,3)
      xphC(ix^D,1)=xgl(ix^D,1)
      xphC(ix^D,2)=xgl(ix^D,2)
      xphC(ix^D,3)=xgl(ix^D,3)+half*dx03
    {end do\}
    ! +++++++++ Special Regular :: 1=2, 2=3, 3=1 +++++++++
    xthCC(ixI^S,1) = xthC(ixI^S,1) * cos(0.50d0*dpi-xthC(ixI^S,2)) * cos(xthC(ixI^S,3))
    xthCC(ixI^S,2) = xthC(ixI^S,1) * cos(0.50d0*dpi-xthC(ixI^S,2)) * sin(xthC(ixI^S,3))
    xthCC(ixI^S,3) = xthC(ixI^S,1) * sin(0.50d0*dpi-xthC(ixI^S,2))
    xphCC(ixI^S,1) = xphC(ixI^S,1) * cos(0.50d0*dpi-xphC(ixI^S,2)) * cos(xphC(ixI^S,3))
    xphCC(ixI^S,2) = xphC(ixI^S,1) * cos(0.50d0*dpi-xphC(ixI^S,2)) * sin(xphC(ixI^S,3))
    xphCC(ixI^S,3) = xphC(ixI^S,1) * sin(0.50d0*dpi-xphC(ixI^S,2))

    print*,'Call prep_rbsl: Caculation Ath ...'
    call RBSL_flux_rope(ixI^L,ixO^L,nnp,a0,F_flx,I_Helix,xthCC,x_axis,AthC)
    call Cart2SphereVector(ixI^L,ixO^L,xthC,AthC,AthS)
    ix1 = 1
    do ix2 = ixOmin2, ixOmax2-1
      do ix3 = ixOmin3, ixOmax3
        Ath_rbsl(ix2,ix3)=AthS(ix1,ix2,ix3,2)*unit_magneticfield*unit_length/1.d5
      end do
    end do

    print*,'Call prep_rbsl: Caculation Aph ...'
    call RBSL_flux_rope(ixI^L,ixO^L,nnp,a0,F_flx,I_Helix,xphCC,x_axis,AphC)
    call Cart2SphereVector(ixI^L,ixO^L,xphC,AphC,AphS)
    ix1 = 1
    do ix2 = ixOmin2, ixOmax2
      do ix3 = ixOmin3, ixOmax3-1
        Aph_rbsl(ix2,ix3)=AphS(ix1,ix2,ix3,3)*unit_magneticfield*unit_length/1.d5
      end do
    end do

    print*,'Call prep_rbsl: Caculation br_rbsl ...'
    r_sun=6.955d5 ! km
    call sinthta_ss(xprobmin2,xprobmax2,m,sinth,sinth_hlf)
    call curlh_ce_ss(m,n,Ath_rbsl,Aph_rbsl,r_sun,sinth,sinth_hlf,dtheta,dphi,br_mfr)

  end subroutine prep_rbsl

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
        end if
        AIx(ix^D,:) = AIx(ix^D,:) + KIr*Rpl(:)
        AFx(ix^D,:) = AFx(ix^D,:) + KFr*Rcr(:)
      end do
      AIx(ix^D,:) = AIx(ix^D,:)*I_cur/a
      AFx(ix^D,:) = AFx(ix^D,:)*F_flx*0.25d0*re_pi/a**2
    {end do\}
    Atotal=AIx+AFx
    if(present(Bfr)) call curlvector(Atotal,ixI^L,ixO^L,Bfr,idirmin,1,ndir)

  end subroutine RBSL_flux_rope

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

  ! initialize solar atmosphere table in a vertical line through the global domain
  subroutine inithdstatic
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    integer          :: j,na,ibc,btlevel
    double precision :: Ta(jmax),gg(jmax)
    double precision :: ytc,rhohc,Tpho,Ttop,htra,ftra,wtra,Ttr,Fc,k_para,hgbot,phc
    double precision :: invT,dyt,rhot,pt,dytt,ratio,tr_argument
    double precision :: res,pb

    call validate_atmosphere_parameters()
    if(.not.all(ieee_is_finite((/unit_length,unit_numberdensity,unit_temperature, &
        unit_pressure,unit_velocity,unit_density,dr,gzone,SRadius,usr_grav/)))) &
      call mpistop('Atmosphere units and grid inputs must be finite')
    if(min(unit_length,unit_numberdensity,unit_temperature,unit_pressure, &
           unit_velocity,unit_density,dr,SRadius)<=0.d0) &
      call mpistop('Atmosphere units, spacing and solar radius must be positive')

     ytc = atmosphere_reference_height_Mm*1.d8/unit_length
   rhohc = atmosphere_reference_nH_cgs/unit_numberdensity
     phc = atmosphere_reference_pressure_cgs/unit_pressure
    Tpho = atmosphere_chromosphere_K/unit_temperature
    Ttop = atmosphere_tanh_top_K/unit_temperature
    htra = atmosphere_transition_height_Mm*1.d8/unit_length
    wtra = atmosphere_transition_width_Mm*1.d8/unit_length
     Ttr = atmosphere_transition_K/unit_temperature
      Fc = atmosphere_conductive_flux_cgs/unit_pressure/unit_velocity
   hgbot = atmosphere_floor_height_Mm*1.d8/unit_length

    ! Preserve the successful RBSL15 analytic temperature profile and
    ! trapezoidal hydrostatic integration. "isothermal" selects only this
    ! initial/boundary atmosphere; it does not disable the MHD energy equation.
    if(atmosphere_model=='analytic_tr') then
      k_para = 8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
      tr_argument=2.d0*(Ttr-Tpho)/(Ttop-Tpho)-1.d0
      if(.not.ieee_is_finite(tr_argument).or.abs(tr_argument)>=1.d0) &
        call mpistop('Atmosphere transition parameters are numerically degenerate')
      ftra = wtra*atanh(tr_argument)
    end if
    do j = 1,jmax
       ya(j) = (dble(j)-0.5d0)*dr-gzone
       if(atmosphere_model=='isothermal') then
         Ta(j) = atmosphere_isothermal_K/unit_temperature
       else if(ya(j)>htra) then
         Ta(j) = (3.5d0*Fc/k_para*(ya(j)-htra)+Ttr**3.5d0)**(2.d0/7.d0)
       else if(ya(j)>hgbot) then
         Ta(j) = Tpho+0.5d0*(Ttop-Tpho)*(tanh((ya(j)-htra+ftra)/wtra)+1.d0)
       else
         Ta(j) = Tpho+0.5d0*(Ttop-Tpho)*(tanh((hgbot-htra+ftra)/wtra)+1.d0)
       end if
       if(SRadius+ya(j)<=0.d0) call mpistop('Atmosphere table reaches a nonpositive radius')
       gg(j) = usr_grav*(SRadius/(SRadius+ya(j)))**2
    enddo
    if(.not.all(ieee_is_finite(Ta)).or.any(Ta<=0.d0)) &
      call mpistop('Atmosphere temperature table must be positive and finite')
    if(.not.all(ieee_is_finite(ya)).or..not.all(ieee_is_finite(gg))) &
      call mpistop('Atmosphere height/gravity table must be finite')
    if(.not.ieee_is_finite(ytc).or.ytc<ya(1).or.ytc>ya(jmax)) &
      call mpistop('Atmosphere reference height lies outside the atmosphere table')

    ra(1) = 1.d5
    pa(1) = ra(1)*Ta(1)
    invT = 0.d0
    do j = 2,jmax
         dyt = ya(j) - ya(j-1)
        invT = invT+dyt*(gg(j)/Ta(j)+gg(j-1)/Ta(j-1))*0.5d0
       pa(j) = pa(1)*dexp(invT)
       ra(j) = pa(j)/Ta(j)
    end do
    if(.not.all(ieee_is_finite(ra)).or..not.all(ieee_is_finite(pa)).or. &
       any(ra<=0.d0).or.any(pa<=0.d0)) &
      call mpistop('Unnormalized hydrostatic density/pressure must be positive and finite')

    ! Both normalizations use the same reference height and spatial
    ! interpolation. The denominator is a height difference, never a pressure difference.
    do j=2,jmax
      if(ytc<=ya(j)) exit
    end do
    dytt=ytc-ya(j-1)
    rhot=ra(j-1)+dytt*(ra(j)-ra(j-1))/(ya(j)-ya(j-1))
    pt=pa(j-1)+dytt*(pa(j)-pa(j-1))/(ya(j)-ya(j-1))
    if(.not.all(ieee_is_finite((/rhot,pt/))).or.rhot<=0.d0.or.pt<=0.d0) &
      call mpistop('Interpolated atmosphere reference values must be positive and finite')
    if(atmosphere_normalization=='number_density') then
      ratio=rhohc/rhot
    else
      ratio=phc/pt
    end if
    if(.not.ieee_is_finite(ratio).or.ratio<=0.d0) &
      call mpistop('Atmosphere normalization ratio must be positive and finite')
    ra=ra*ratio
    pa=pa*ratio
    if(.not.all(ieee_is_finite(ra)).or..not.all(ieee_is_finite(pa)).or. &
       any(ra<=0.d0).or.any(pa<=0.d0)) &
      call mpistop('Normalized hydrostatic density/pressure must be positive and finite')

    if(mype==0) then
      write(*,*) 'Atmosphere model: ',trim(atmosphere_model)
      write(*,*) 'Atmosphere normalization: ',trim(atmosphere_normalization)
      write(*,*) 'Atmosphere reference height [Mm]: ',atmosphere_reference_height_Mm
      write(*,*) 'Atmosphere requested nH [cm^-3], pressure [erg cm^-3]: ', &
                  atmosphere_reference_nH_cgs,atmosphere_reference_pressure_cgs
      write(*,*) 'Atmosphere actual reference nH, pressure: ', &
                  rhot*ratio*unit_numberdensity,pt*ratio*unit_pressure
    end if

    ! initialized rho and p in the fixed bottom boundary
    !! NOTE: just to find the value of gzone position
      na = floor(gzone/dr+0.5d0)      
     res = gzone-(dble(na)-0.5d0)*dr
    rhob = ra(na)+res/dr*(ra(na+1)-ra(na))
      pb = pa(na)+res/dr*(pa(na+1)-pa(na))
    allocate(rbc(nghostcells))
    allocate(pbc(nghostcells))

    !! Define rho,p of ghostcells
    btlevel=refine_max_level
    do ibc=nghostcells,1,-1
       na = floor((gzone-dx(1,btlevel)*(dble(nghostcells-ibc+1)-0.5d0))/dr+0.5d0)
      res = gzone-dx(1,btlevel)*(dble(nghostcells-ibc+1)-0.5d0)-(dble(na)-0.5d0)*dr
      rbc(ibc) = ra(na)+res/dr*(ra(na+1)-ra(na))
      pbc(ibc) = pa(na)+res/dr*(pa(na+1)-pa(na))
    end do
    
    if (mype==0) then
     print*,'Ratio, rhot:',ratio,rhot
     print*,'Range of rho solution:',minval(ra),maxval(ra)
     print*,'Range of pressure solution:',minval(pa),maxval(pa)
     print*,'Minimum of rho solution:',minval(ra)*unit_numberdensity
     print*,'Minimum of pressure solution:',minval(pa)*unit_pressure
     print*,'Maximum of Temperature solution:',Ta(jmax)*unit_temperature
     print*,'Rho at fixed bottom boundary:',minval(rbc)*unit_numberdensity, maxval(rbc)*unit_numberdensity
     print*,'P at fixed bottom boundary:',minval(pbc)*unit_pressure, maxval(pbc)*unit_pressure
    end if
  end subroutine inithdstatic
  
  subroutine rbsl_insert(ixIs^L, ixI^L, ixO^L, ws, x, s_in)                                

    integer, intent(in)                :: ixIs^L, ixI^L, ixO^L
    double precision, intent(inout)    :: ws(ixIs^S,1:nws)
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    type(state),intent(in)             :: s_in
                                                                            
    double precision                   :: Adummy(ixIs^S,1:3)                   
                                                                            
    Adummy=zero                                                                  
    call b_from_vector_potential_S(ixIs^L, ixI^L, ixO^L, ws, x, Adummy, s_in)
                                                                            
  end subroutine rbsl_insert

  !> calculate magnetic field from vector potential A at cell edges
  subroutine b_from_vector_potential_S(ixIs^L, ixI^L, ixO^L, ws, x, A, s_in)

    integer, intent(in)                :: ixIs^L, ixI^L, ixO^L
    double precision, intent(inout)    :: ws(ixIs^S,1:nws),A(ixIs^S,1:3)
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    type(state),intent(in)             :: s_in

    integer                            :: ixC^L, hxC^L, idim1, idim2, idir
    double precision                   :: xC(ixIs^S,1:ndim),xCC(ixIs^S,1:ndim)
    double precision                   :: circ(ixIs^S,1:ndim)

    A=zero
    ! extend one layer of cell center locations in xCC
    xC=0.d0
    xCC=0.d0
    xCC(ixI^S,1:ndim)=x(ixI^S,1:ndim)
    {
    xCC(ixIsmin^D^D%ixI^S,1:ndim)=x(ixImin^D^D%ixI^S,1:ndim)
    xCC(ixIsmin^D^D%ixIs^S,^D)=x({ixImin^DD,},^D)-s_in%dx({ixImin^DD,},^D)
    \}
    {^IFTHREED
    xCC(ixImin1:ixImax1,ixIsmin2,ixIsmin3,1)=x(ixImin1:ixImax1,ixImin2,ixImin3,1)
    xCC(ixIsmin1,ixImin2:ixImax2,ixIsmin3,2)=x(ixImin1,ixImin2:ixImax2,ixImin3,2)
    xCC(ixIsmin1,ixIsmin2,ixImin3:ixImax3,3)=x(ixImin1,ixImin2,ixImin3:ixImax3,3)
    }

    do idir=7-2*ndim,3
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-1+kr(idir,^D);
      do idim1=1,ndim
        ! Get edge coordinates
        if (idim1/=idir) then
          xC(ixC^S,idim1)=xCC(ixC^S,idim1)+half*s_in%dx(ixC^S,idim1)
        else
          xC(ixC^S,idim1)=xCC(ixC^S,idim1)
        end if
      end do
      ! Initialise vector potential at the edge
      call rbsl_vectorpotential(ixIs^L, ixC^L, xC, A(ixIs^S,idir), idir, s_in)
      A(ixC^S,idir)=A(ixC^S,idir)*s_in%dsC(ixC^S,idir)
    end do

    ! Take the curl of the vector potential 
    circ=zero
    ! Calculate circulation on each face
    do idim1=1,ndim ! Coordinate perpendicular to face 
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          if(lvc(idim1,idim2,idir)==0) cycle
          ! Assemble indices
          hxC^L=ixC^L-kr(idim2,^D);
          ! Add line integrals in direction idir
          circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                           +lvc(idim1,idim2,idir)* &
                            (A(ixC^S,idir)&
                            -A(hxC^S,idir))
        end do
      end do
      ! Divide by the area of the face to get B
      where(s_in%surfaceC(ixC^S,idim1)==0)
        circ(ixC^S,idim1)=zero
      elsewhere
        circ(ixC^S,idim1)=circ(ixC^S,idim1)/s_in%surfaceC(ixC^S,idim1)
      end where
      ws(ixC^S,idim1) = ws(ixC^S,idim1)+circ(ixC^S,idim1)
    end do

  end subroutine b_from_vector_potential_S

  subroutine rbsl_vectorpotential(ixI^L, ixC^L, xC, A, idir, s_in)
    ! calculate the rbsl vectorpotential on the edges
    ! used by b_from_vectorpotential_S()
    integer, intent(in)                :: ixI^L, ixC^L,idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)
    type(state),intent(in)             :: s_in

    integer :: ix^D,ixg^D,dix
    double precision :: xlen^D
    double precision :: xCC(ixI^S,1:ndim)
    double precision :: A_vec(ixI^S,1:ndim)
    double precision :: A_rbsl(ixI^S,1:ndim)

    A=zero
    xCC(ixI^S,1) = xC(ixI^S,1) * cos(half*dpi-xC(ixI^S,2)) * cos(xC(ixI^S,3))
    xCC(ixI^S,2) = xC(ixI^S,1) * cos(half*dpi-xC(ixI^S,2)) * sin(xC(ixI^S,3))
    xCC(ixI^S,3) = xC(ixI^S,1) * sin(half*dpi-xC(ixI^S,2))
    call RBSL_flux_rope(ixI^L,ixC^L,nnp,a0,F_flx,I_Helix,xCC,x_axis,A_vec)
    call Cart2SphereVector(ixI^L,ixC^L,xC,A_vec,A_rbsl)
    if (idir==3) then
      A(ixC^S)=A(ixC^S)+A_rbsl(ixC^S,3)
    else if(idir==2) then
      A(ixC^S)=A(ixC^S)+A_rbsl(ixC^S,2)
    else
      A(ixC^S)=A(ixC^S)+A_rbsl(ixC^S,1)
    end if

  end subroutine rbsl_vectorpotential

  !> calculate cell-center values from face-center values in fake process
  subroutine fake_face_to_center(ixIs^L,ixO^L,ws,s)
    use mod_global_parameters
    ! Non-staggered interpolation range
    integer, intent(in)                :: ixIs^L, ixO^L
    double precision, intent(in)       :: ws(ixIs^S,1:nws)
    type(state)                        :: s

    integer                            :: fxO^L, gxO^L, hxO^L, jxO^L, kxO^L, idim

    associate(w=>s%w)

    ! calculate cell-center values from face-center values in 2nd order
    do idim=1,ndim
      ! Displace index to the left
      ! Even if ixI^L is the full size of the w arrays, this is ok
      ! because the staggered arrays have an additional place to the left.
      hxO^L=ixO^L-kr(idim,^D);
      ! Interpolate to cell barycentre using arithmetic average
      ! This might be done better later, to make the method less diffusive.
      w(ixO^S,mag(idim))=half/s%surface(ixO^S,idim)*&
        (ws(ixO^S,idim)*s%surfaceC(ixO^S,idim)&
        +ws(hxO^S,idim)*s%surfaceC(hxO^S,idim))
    end do

    ! calculate cell-center values from face-center values in 4th order
    !do idim=1,ndim
    !  gxO^L=ixO^L-2*kr(idim,^D);
    !  hxO^L=ixO^L-kr(idim,^D);
    !  jxO^L=ixO^L+kr(idim,^D);

    !  ! Interpolate to cell barycentre using fourth order central formula
    !  w(ixO^S,mag(idim))=(0.0625d0/s%surface(ixO^S,idim))*&
    !         ( -ws(gxO^S,idim)*s%surfaceC(gxO^S,idim) &
    !     +9.0d0*ws(hxO^S,idim)*s%surfaceC(hxO^S,idim) &
    !     +9.0d0*ws(ixO^S,idim)*s%surfaceC(ixO^S,idim) &
    !           -ws(jxO^S,idim)*s%surfaceC(jxO^S,idim) )
    !end do

    ! calculate cell-center values from face-center values in 6th order
    !do idim=1,ndim
    !  fxO^L=ixO^L-3*kr(idim,^D);
    !  gxO^L=ixO^L-2*kr(idim,^D);
    !  hxO^L=ixO^L-kr(idim,^D);
    !  jxO^L=ixO^L+kr(idim,^D);
    !  kxO^L=ixO^L+2*kr(idim,^D);

    !  ! Interpolate to cell barycentre using sixth order central formula
    !  w(ixO^S,mag(idim))=(0.00390625d0/s%surface(ixO^S,idim))* &
    !     (  +3.0d0*ws(fxO^S,idim)*s%surfaceC(fxO^S,idim) &
    !       -25.0d0*ws(gxO^S,idim)*s%surfaceC(gxO^S,idim) &
    !      +150.0d0*ws(hxO^S,idim)*s%surfaceC(hxO^S,idim) &
    !      +150.0d0*ws(ixO^S,idim)*s%surfaceC(ixO^S,idim) &
    !       -25.0d0*ws(jxO^S,idim)*s%surfaceC(jxO^S,idim) &
    !        +3.0d0*ws(kxO^S,idim)*s%surfaceC(kxO^S,idim) )
    !end do

    end associate

  end subroutine fake_face_to_center

  !==============================================================================
  ! Purpose: to read Eth and Eph from input fits files 
  !==============================================================================
  subroutine read_electric(tc,mi,ni,Eth0,Eph0,Era0,t_rec)
    use mod_global_parameters

    integer, intent(in)           :: tc,mi,ni
    double precision, intent(out) :: Eth0(mi,ni+1), Eph0(mi+1,ni), Era0(mi+1,ni+1)
    double precision, intent(out) :: t_rec

    real :: Eth_tmp(mi,ni+1), Eph_tmp(mi+1,ni), Era_tmp(mi+1,ni+1)
    real :: nullval=0.0
    integer :: t_out,naxes(2),naxis,p0
    integer :: status,readwrite,blocksize
    logical :: anynull
    character*80 :: comment
    character*200 :: filename_Eth, filename_Eph, filename_Era

    write (filename_Eth,  '(a,i4.4,a)') 'Epdfi/Eth',  tc,'.fits'
    write (filename_Eph,  '(a,i4.4,a)') 'Epdfi/Eph',  tc,'.fits'
    write (filename_Era,  '(a,i4.4,a)') 'Epdfi/Era',  tc,'.fits'

    status = 0
    blocksize = 1

    ! read Eth
    call ftopen(127,filename_Eth,readwrite,blocksize,status)
    call ftgpve(127,1,1,p_eth,nullval,Eth_tmp,anynull,status)
    call ftgkyj(127, 'TIME', t_out, comment, status)
    call ftclos(127, status)
    call ftfiou(127, status)
    t_rec = real(t_out,8)/unit_time
    Eth0  = dble(Eth_tmp)/(unit_velocity*unit_magneticfield)

    ! read Eph
    call ftopen(128,filename_Eph,readwrite,blocksize,status)
    call ftgpve(128,1,1,p_eph,nullval,Eph_tmp,anynull,status)
    call ftclos(128, status)
    call ftfiou(128, status)
    Eph0  = dble(Eph_tmp)/(unit_velocity*unit_magneticfield)

    ! read Era
    call ftopen(129,filename_Era,readwrite,blocksize,status)
    call ftgpve(129,1,1,p_era,nullval,Era_tmp,anynull,status)
    call ftclos(129, status)
    call ftfiou(129, status)
    Era0  = dble(Era_tmp)/(unit_velocity*unit_magneticfield)

  end subroutine read_electric

  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    integer :: ix1
    do ix1=ixOmin1,ixOmax1
      ggrid(ix1^%1ixO^S)=usr_grav*(SRadius/x(ix1^%1ixO^S,1))**2
    end do
  end subroutine

  subroutine getbQ(bQgrid,ixI^L,ixO^L,qt,w,x)
    !!calculate background heating bQ, Mok 2016 ApJ
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)
    double precision, intent(out) :: bQgrid(ixI^S)

    double precision :: bQgrid0(ixI^S),bQgrid1(ixI^S)
    double precision :: Bmag(ixI^S),bvec(ixI^S,1:ndir),cvec(ixI^S,1:ndir),tmp(ixI^S)
    double precision :: curlo,curhi
    integer :: idims,idir,ix^D

    bQgrid(ixO^S)=bQ0*dexp(-(x(ixO^S,1)-xprobmin1)/heatscale)
    !bQgrid0(ixO^S)=bQ0*dexp(-(x(ixO^S,1)-xprobmin1)/heatscale)
    !if(B0field) then
    !  Bmag(ixI^S)=dsqrt(sum((w(ixI^S,mag(:))+block%B0(ixI^S,:,0))**2,dim=ndim+1))
    !  do idir=1,ndir
    !    bvec(ixI^S,idir)=(w(ixI^S,mag(idir))+block%B0(ixI^S,idir,0))/Bmag(ixI^S)
    !  end do
    !else
    !  Bmag(ixI^S)=dsqrt(sum(w(ixI^S,mag(:))**2,dim=ndim+1))
    !  do idir=1,ndir
    !    bvec(ixI^S,idir)=w(ixI^S,mag(idir))/Bmag(ixI^S)
    !  end do
    !end if
    !cvec=0.d0
    !! calculate local curvature of magnetic field
    !do idims=1,ndim
    !  call gradient(bvec(ixI^S,1),ixI^L,ixO^L,idims,tmp) 
    !  cvec(ixO^S,1)=cvec(ixO^S,1)+bvec(ixO^S,idims)*tmp(ixO^S)
    !  call gradient(bvec(ixI^S,2),ixI^L,ixO^L,idims,tmp) 
    !  cvec(ixO^S,2)=cvec(ixO^S,2)+bvec(ixO^S,idims)*tmp(ixO^S)
    !  call gradient(bvec(ixI^S,3),ixI^L,ixO^L,idims,tmp) 
    !  cvec(ixO^S,3)=cvec(ixO^S,3)+bvec(ixO^S,idims)*tmp(ixO^S)
    !end do
    !tmp(ixO^S)=dsqrt(sum(cvec(ixO^S,:)**2,dim=ndim+1))
    !! set lower and upper limit for curvature
    !curlo=0.1d0 ! 1/10 
    !curhi=2.d0 ! 1/0.5
    !where(tmp(ixO^S)<curlo)
    !  tmp(ixO^S)=curlo
    !elsewhere(tmp(ixO^S)>curhi)
    !  tmp(ixO^S)=curhi
    !end where
    !bQgrid1(ixO^S)=bQ0*Bmag(ixO^S)**1.75d0*w(ixO^S,rho_)**0.125d0*tmp(ixO^S)**0.75d0
    !
    !{do ix^DB=ixOmin^DB,ixOmax^DB\}
    !    if(bQgrid1(ix^D)<bQgrid0(ix^D)) then
    !       bQgrid(ix^D)=bQgrid0(ix^D)
    !    else
    !       bQgrid(ix^D)=bQgrid1(ix^D)
    !    end if
    !{end do\}

  end subroutine getbQ

end module mod_usr
