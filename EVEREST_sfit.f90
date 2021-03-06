module fit_module
!
#define debug_  1
!
use input
use timer
use accuracy
use caoh_param
!
implicit none
!
integer  :: verbose  = 3 ! Verbosity level    
!
public fitting_energies
!
private
!
integer,parameter :: Nestates = 3,ncouples = 1
!
integer(ik)  :: total_parameters = 0 
!
integer, parameter :: enermax=6000    ! maximal number of calculated energies for a given J and a given symmetry 
!integer, parameter :: obsmax=5000     ! maximal number of experimental energies
integer, parameter :: pointsmax=7000  ! maximal number of potential data points 
!
integer, parameter :: f_inp = 5, f_out = 6   ! read/write units 
!
integer,parameter          ::  findif = 3    ! Dinite difference differentiation with 2 or 3 points
!
real(rk),parameter :: fitfactordeltax=0.001   ! parameter for the finite differncies differention 
!
character(len=70) :: deriv_type = 'finite  '  ! hellmann or direct - how we calculate derivatives. 
                                             ! direct means the finite diferentiation of energies wrt parameters.
!
integer,parameter :: nofititer=6             ! the Jacobi matrix can be evaluated only each nofititer+1 step. For all 
                                             ! iterations in between the the same jacobi matrix is used 
integer,parameter :: f_en =206               ! All computed energies (not only those that have obs. counterparts) will be written into here.
integer,parameter :: f_pot=207               ! The current potential parameters will be written here. 
integer,parameter :: f_potpoints =208        ! Computed at the current iteration step potential energy points are printed here. 
character(len=200)::  char_job(10)           ! The derectives from the input file are stored in this array. Will be used to create all job-input files.
character(len=200)::  char_vib_job(100)      ! The vib. derectives from the input file are stored in this array. Will be used to create all job-input files.
character(len=200)::  char_rot_job(100,2)      ! The rot. derectives from the input file are stored in this array. Will be used to create all job-input files.
character(len=1)  :: mark
!
integer(ik) :: Nfields = 4
!
integer,parameter :: sym2parity(4) = (/0,1,0,1/) ! parity versus symmetry
!
! global variables, as they appear in the Everest input file.
!
integer ::   npnt2_,jrot_,neval_,nalf_,max2d_,max3d_,idia_,&
             kmin_,npnt1_,ipar_,max3d2_
integer ::   ibass_ = 0, nvib_, neval2_ = 0
integer      :: npnt2x(100),nevalx(100),jrot_card(100),icard_j(0:100)
!
!integer ::   jmax = 10, kmax =10    ! maximal value up to whicj all j = 0,1,2,..jmax will be evaluated. 
!
real(8) ::   ezero_ = -1                     ! zero point energy as computed at J=0,A1 
!
integer             :: ierr
character(len=wl)   :: line_buffer
integer,parameter   :: max_input_lines=500000  ! maximum length (in lines) of input. 500,000 lines is plenty..
!                                              ! to avoid feeding in GB of data by mistake.

  !
  ! type describing the parameter geneology
  !
  type linkT
    !
    integer(ik) :: ifield
    integer(ik) :: iparam
    !
  end type linkT
  !
  type fieldT
    !
    character(len=cl)    :: name='(unnamed)' ! Identifying name of the function (default to avoid garbled outputs)
    character(len=cl)    :: type='NONE'  ! Identifying type of the function
    character(len=cl)    :: class='NONE' ! Identifying class of the function (poten, spinorbit,dipole,abinitio etc)
    !
    ! variable used for GRID curves to specify interpolation and extrapolation  options
    ! 
    character(len=cl)    :: interpolation_type='CUBICSPLINES'
    !
    !character(len=cl)    :: interpolation_type='QUINTICSPLINES'
    !
    integer(ik)          :: iref         ! reference number of the term as given in input (bra in case of the coupling)
    integer(ik)          :: jref         ! reference number of the coupling term as given in input (ket in case of the coupling)
    integer(ik)          :: istate       ! the actual state number (bra in case of the coupling)
    integer(ik)          :: jstate       ! the actual state number (ket in case of the coupling)
    integer(ik)          :: Nterms       ! Number of terms or grid points
    real(rk)             :: factor=1.0_rk      ! defines if the term is imaginary or real
    real(rk)             :: r1,r2,alpha        ! equilibrium values
    real(rk)             :: alphamin           ! equilibrium values
    real(rk)             :: range(3,10)        ! ranges of values
    real(rk)             :: fit_factor=1.0     ! defines if the term is imaginary or real
    real(rk),pointer     :: value(:)=>null()   ! Expansion parameter or grid values from the input
    real(rk),pointer     :: weight(:)=>null()    ! fit (1) or no fit (0)
    integer(ik),pointer     :: ipower1(:)=>null()    ! power
    integer(ik),pointer     :: ipower2(:)=>null()    ! power
    integer(ik),pointer     :: ipower3(:)=>null()    ! power
    character(len=cl),pointer :: forcename(:)=>null() ! The parameter name
    integer(ik),pointer     :: iaddress(:)=>null()    ! iaddress in the potparam object 
    type(linkT),pointer     ::  link(:)=>null() ! linking to a value of a different field
    !
  end type fieldT
  !
  type quantaT
    real(rk)     :: Jrot       ! J - real
    integer(ik)  :: irot       ! index of the J value in J_list
    integer(ik)  :: Ka         ! index of the J value in J_list
    integer(ik)  :: istate     ! e-state
    integer(ik)  :: imulti     ! imulti = 1,..,multiplicity
    real(rk)     :: F          ! spin-projection = -spin,-spin+1,..,spin-1,spin
    real(rk)     :: omega      ! sigma+lambda
    real(rk)     :: spin       ! spin
    integer(ik)  :: ilambda    ! ilambda
    integer(ik)  :: v1  = 0     ! vibrational quantum
    integer(ik)  :: v2  = 0     ! vibrational quantum
    integer(ik)  :: v3  = 0     ! vibrational quantum
    integer(ik)  :: ivib = 1   ! vibrational quantum number counting all vib-contracted states
    integer(ik)  :: ilevel = 1  ! primitive quantum
    integer(ik)  :: iroot       ! running number
    integer(ik)  :: iJ_ID       ! running number within the same J
    integer(ik)  :: iparity = 0
    integer(ik)  :: igamma = 1
    integer(ik)  :: iomega = 1  ! countig number of omega
    character(len=cl) :: name   ! Identifying name of the  function
  end type quantaT
  !
  type obsT
     !
     real(rk)      :: Jrot
     real(rk)      :: Jrot_      ! J - real (lower)
     integer(ik)   :: irot       ! index of the J value in J_list
     integer(ik)   :: irot_
     integer(ik)   :: symmetry
     integer(ik)   :: N
     integer(ik)   :: N_         ! N lower
     integer(ik)   :: iparity
     integer(ik)   :: iref  ! reference group:  1 and 2
     integer(ik)   :: iparity_   ! parity +/- (0,1) of the lower state lower
     real(rk)      :: energy
     real(rk)      :: frequency
     real(rk)      :: weight
     type(quantaT) :: quanta
     type(quantaT) :: quanta_  ! quantum numbers of the lower state
     !
  end type obsT
  !
    type calcT
     !
     real(rk),pointer   :: energy(:)=>null()
     type(quantaT),pointer  :: quanta(:)=>null()
     !
  end type calcT
  !
  type fittingT
     !
     logical              :: run
     integer(ik)          :: nJ = 1        ! Number of J values processed 
     real(rk),pointer     :: j_list(:)=>null()     ! J-values processed in the fit
     integer(ik)          :: iparam(1:2) = (/1,100000/)
     integer(ik)          :: itermax = 500
     integer(ik)          :: Nenergies = 1
     integer(ik)          :: Nabinitio = 1
     integer(ik)          :: parmax =0          ! total number of all parameters used
     real(rk)             :: factor = 1.0_rk
     real(rk)             :: target_rms = 1e-8
     real(rk)             :: target_stability = 1e-12
     real(rk)             :: robust = 0
     character(len=cl)    :: geom_file = 'pot.fit'
     character(len=cl)    :: output_file = 'fitting'
     character(len=cl)    :: fit_type = 'LINUR'      ! to switch between fitting methods.
     real(rk)             :: threshold_coeff    = -1e-18
     real(rk)             :: threshold_lock     = -1e-18
     real(rk)             :: threshold_obs_calc  = -1e-16
     real(rk)             :: zpe=0
     logical              :: shift_to_zpe = .true.   ! 
     real(rk)             :: fit_scaling=1.0_rk         ! scaling the fitting correction with this factor >0 and <1
     integer(ik)          :: linear_search = 0 ! use linear scaling to achieve better convergence with the Armijo condition
     real(rk)             :: Vmin(1:4) = 0 
     type(obsT),pointer   :: obs(:)=>null()         ! experimental data
     real(rk),pointer    :: ai(:)=>null()           ! ab initio energies
     real(rk),pointer    :: r1(:)=>null()           ! ab initio geometry
     real(rk),pointer    :: r2(:)=>null()           ! ab initio geometry
     real(rk),pointer    :: r3(:)=>null()           ! ab initio geometry
     real(rk),pointer    :: weight(:)=>null()       ! ab initio geometry
     integer(ik),pointer :: iPES(:)=>null()         ! ab initio geometry
     !
     !type(paramT),pointer :: param(:)         ! fitting parameters
     !
  end type fittingT
  !
  type(fieldT),pointer :: poten(:)=>null(),spinorbit(:)=>null()
  !
  type(fittingT)            :: fitting

!
contains
 !
subroutine fitting_energies
 
  use input 
  use caoh_param 

  real(rk)          :: jrot                    ! current value of the rotational quantum number 
  !
  character(len=30) :: evib_exe =   "evvib.e"  ! Everest exe-file. 
  character(len=30) :: erot_exe =   "evrot.e"  ! Everest exe-file. 
  character(len=30) :: dervib_exe = "evvder.e" ! Everest exe-file: vibrational derivatives 
  character(len=30) :: derrot_exe = "evrder.e"          ! Everest exe-file: vibrational derivatives 

  integer           :: alloc, ierror      ! error state variables 
  !
  integer,parameter  :: jlist_max = 500
  integer(ik) :: maxiter_as = 3     ! maximal number of iterations to find a match for assignement  !
  !
  real(rk)              :: j_list_(1:jlist_max)=-1.0_rk,jmin,jmax,kmax,jrot2,f_t
  !
  character(len=3)   :: kmax_ch,jmax_ch
  character(len=2)   :: fititer_ch_
  character(len=2)   :: fititer_ch
  !
  real(rk),allocatable  :: j_list(:)
  !
  character(len=cl) :: w
  !
  logical :: eof
  !
  integer(ik) :: i,nJ,ipot,iso,iref,istate,jref,Nparam,Nparam_check,iparam,ifield,ifield_,iaddr,Total_size=0
  integer(ik) :: iOMP = 1
  !
  type(fieldT),pointer      :: field
  !
  integer(ik)       :: iut !  iut is a unit number. 
  integer(ik)  :: iobs,itau,iai,ic,iparams,irange,iparity,ipar,irot,Nsize,iaddress,iterm,k0,MXRoot
  integer(ik)  :: itau_,iter_th,iener_,jobs
  logical :: matchfound
  character(len=wl) :: large_fmt
  character(len=cl) :: ioname
  !
  integer(ik)  :: nmax,nrot
  integer(ik),allocatable  :: Nstates(:,:,:)
  type(calcT),allocatable  :: calc(:,:,:)
  character(len=1),allocatable  :: mark(:)


  integer           :: npts,pot_npts,nused      ! number of data points: current, all obs. energies, all pot. points, and actuall used. 
  integer           :: nrow,ncol 
  integer           :: numpar                  ! number of varying parameters 
  integer           :: ndigits                 ! number of rounded digits in the parameters output. 
  !
  ! objects to store the pot. parameters, geometries, values of poten. function.  
  !
  real(rk),allocatable :: parold(:),potparam(:)
  !
  !integer,allocatable:: ipower(:,:)          ! powers of the potential function expansion. Can be used 
  !                                           ! to form the potential function os just as a parameter name.
  character(len=8),allocatable :: nampar(:)   ! parameter names 
  
  integer,allocatable :: ivar(:)              ! switch for the parameters fit: ivar(i) = 0 for the parameters that do not vary.
  integer,allocatable ::   ifitparam(:),iparamfit(:)   ! address of the refiedn parameters
  !
  !integer,allocatable :: J_obs(:),sym_obs(:) ! Arrays where we 
  !                                           ! store the information about obs. energies: 
  !                                           ! J, symmetry, number in the block, vib. quantum numbers
  !
  ! objects to store the energies 
  real(rk),allocatable :: enercalc(:),ener_obs(:),eps(:)
  !
  ! objects to needed for the fitting procedure: jacobi matrix, derivatives, matrices to solve the 
  ! the ststem of linear equations, standard errors and so on. 
  real(rk),allocatable :: rjacob(:,:),al(:,:),bl(:),dx(:),ai(:,:),sterr(:),derj0(:)
  real(rk),allocatable :: sigma(:)

  real(rk),allocatable :: wt_bit(:),wtall(:) ! weight factors - only for the energies and total.
  !
  !
  real(rk) :: wtsum,fit_factor,r1,r2,theta,ssq,rms,sum_sterr,conf_int
  real(rk) :: stadev_old,stability,stadev,tempx,deltax,v,potright,potleft
  real(rk) :: ssq1,ssq2,rms1,rms2
  logical          :: still_run
  integer          :: j,l,info_
  integer          :: iener,irow,icolumn,ipar_
  !                   and their weights will we reduced. 
  !
  ! different character objects 
  !
  character(len=68)  :: fmt1 ! to have nice parameters output   
  character(len=2)   :: fmt0
  character(len=7)   :: out_name    ! file name for outputs
  !
  ! parameters needed for the dgelss lapack routine: 
  !
  real(rk), allocatable :: Tsing(:,:)
  integer                       :: rank0,lwork
  real(rk),allocatable  :: wspace(:)
  !
  integer                     :: itime0,itime2,itime,irate2,imax2 ! time control variables 
  !
  integer                     :: fititer   ! iteration variable
  !
  logical                     :: ifopen,isys
  !
  logical                     ::  SYSTEMQQ  ! system function for  calling extyernal programs, 
                                            ! platform dependent.  
   !
   ! START PROGRAM HERE 
   !
   if (verbose>=2) write(f_out,"(/'The fitting procedure starts...'/)")
   !
   ! read the beginning time. 
   !
   call accuracyInitialize
   !
   call SYSTEM_CLOCK(itime0,irate2,imax2)
   !
   !dir = 'C:\sergei\programs\fiiting_dvr\myprogs\fit\'
   !
   ! Read the input file, count all energies, parameters, data points, so that we can 
   ! allocate all needed arrays and only then actuly read the input again. 
   !
   allocate(poten(Nfields),stat=alloc)
   call ArrayStart('poten',alloc,1,4)
   !allocate(spinorbit(ncouples),stat=alloc)
   !call ArrayStart('spinorbit',alloc,1,4)
   !
   ! input.f90
   !
   write(ioname, '(a, i4)') 'write to a temporary (scratch) file.'
   call IOstart(trim(ioname), iut)
   !
   open(unit=iut, status='scratch', action='readwrite')
   write(large_fmt, '(A,i0,A)') '(A', wl, ')'
   trans_loop: do i=1, max_input_lines
     read(unit=*,fmt=large_fmt,iostat=ierr) line_buffer
     if(ierr /=0) exit trans_loop
     write(iut, '(a)') trim(line_buffer)

     ! This is a hack; I need to know if to echo the input or not before processing it
     ! The option 'do_not_echo_input' is dealt with here as a special case
     line_buffer = adjustl(line_buffer) ! remove leading spaces

     do j=1, len(trim(line_buffer)) ! convert to uppercase
      ic = ichar( line_buffer(j:j))
      if( ic >= 97) line_buffer(j:j) = achar(ic-32)
     enddo
   enddo trans_loop
   rewind(iut)
   !
   call input_options(echo_lines=.true.,error_flag=1)
   !
   ipot = 0
   iso = 0
   total_parameters = 0
   irange = 0
   nmax = 0
   jmax = 0.5
   kmax = 0.5
   !
   do
      call read_line(eof,iut) ; if (eof) exit
      call readu(w)
      select case(w)
      !
      case("STOP","FINISH","END")
        exit
        !
      case("")
       print "(1x)"    !  Echo blank lines
        !
      case('VERBOSE')
        !
        call readi(verbose)
        !
      case('J_LIST','JLIST','JROT','J')
        !
        jmin =  1e6
        jmax = -1.0
        i = 0
        !
        do while (item<Nitems.and.i<jlist_max)
           !
           i = i + 1
           !
           call readu(w)
           !
           if (trim(w)/="-") then
             !
             read(w,*) jrot
             !
             j_list_(i) = jrot
             !
           else
             !
             call readf(jrot2)
             !
             do while (item<=Nitems.and.nint(2.0*jrot)<nint(2.0*jrot2))
               !
               jrot = jrot + 1.0
               j_list_(i) = jrot
               i = i + 1
               !
             enddo
             i = i - 1
             !
           endif
           !
           jmin = min(jmin,j_list_(i))
           jmax = max(jmax,j_list_(i))
           !
        enddo
        !
        nJ = i
        !
        allocate(j_list(i),stat=alloc)
        !
        call ArrayStart('j_list',alloc,size(j_list),kind(j_list))
        !
        J_list(1:i) = J_list_(1:i)
        !
      case('KMAX')
        !
        call readf(kmax)
        !
      case('FIT_FACTOR')
        !
        call readf(fit_factor)
        !
      case('OUTNAME','OUT_NAME')
        !
        call reada(out_name)
        !
      case('EVIB')
        !
        call reada(evib_exe)
        !
      case('EROT')
        !
        call reada(erot_exe)
        !
      case('DERVIB','DER_VIB')
        !
        call reada(dervib_exe)
        !
      case('DERROT','DER_ROT')
        !
        call reada(derrot_exe)
        !
      case('OMP')
        !
        call readi(iOMP)
        !
      case ("MEM","MEMORY")
        !
        call readf(memory_limit)
        !
        call readu(w)
        !
        select case(w)
            !
          case default 
            !
            call report("Unexpected argument in MEMORY",.true.)
            !
          case("TB","T")
            !
            memory_limit = memory_limit*1024_rk
            !
          case("GB","G")
            !
            memory_limit = memory_limit
            !
          case("MB","M")
            !
            memory_limit = memory_limit/1024_rk
            !
          case("KB","K")
            !
            memory_limit = memory_limit/1024_rk**2
            !
          case("B")
            !
            memory_limit = memory_limit/1024_rk**3
            !
        end select
        !
      case("SPIN-ORBIT","POTEN","POTENTIAL") 
          !
          select case (w)
             !
          case("POTEN","POTENTIAL")
             !
             ipot = ipot + 1
             !
             field => poten(ipot)
             field%istate = ipot
             !
             call readi(field%iref)
             field%jref = field%iref
             field%class = "POTEN"
             !
             ! Check if it was defined before 
             do istate=1,ipot-1
                if (field%iref==poten(istate)%iref) then
                  call report ("poten object is repeated",.true.)
                endif
             enddo
             !
          case("SPIN-ORBIT")
             !
             iso = iso + 1
             !
             call readi(iref)
             call readi(jref)
             !
             field => spinorbit(iso)
             !
             field%iref = iref
             field%jref = jref
             field%istate = iref
             field%jstate = jref
             !
             field%class = "SPINORBIT"
             !
          case default
             call report ("Unrecognized keyword (error 02): "//trim(w),.true.)
          end select
          !
          call read_line(eof,iut) ; if (eof) exit
          call readu(w)
          !
          irange = 0 
          !
          do while (trim(w)/="".and.trim(w)/="END") !we read until the END the object input block
            !
            select case(w)
              !
            case("TYPE")
              !
              call readu(w)
              !
              field%type = trim(w)
              !
            case("NAME")
              !
              call reada(w)
              !
              field%name = trim(w)
              !
            case("R1")
              !
              call readf(field%r1)
              !
            case("R2")
              !
              call readf(field%r2)
              !
            case("ALPHA")
              !
              call readf(field%alpha)
              !
            case("ALPHAMIN")
              !
              call readf(field%alphamin)
              !
            case("RANGE")
              !
              irange = irange + 1
              !
              call readf(field%range(1,irange))
              call readf(field%range(2,irange))
              call readf(field%range(3,irange))
              !
            case("FIT_FACTOR")
              !
              call readf(field%fit_factor)
              !
            case("NPARAM","N","NPOINTS")
              !
              call readi(Nparam)
              !
              if (Nparam<0) then
                  call report ("The size of the potential is illegar (<1)",.true.)
              endif
              !
              field%Nterms = Nparam
              !
            case("VALUES")
              !
              Nparam_check = 0
              !
              call input_options(echo_lines=.true.,error_flag=1)
              !
              do while (trim(w)/="END")
                 !
                 call read_line(eof,iut) ; if (eof) exit
                 !
                 call readu(w)
                 !
                 Nparam_check = Nparam_check+1
                 !
              enddo
              !
              Nparam_check = Nparam_check-1
              !
              if (trim(w) /= "END") then
                  call report ("ERROR: Cannot find `END' statement)",.true.)
              endif
              !
              ! go back to beginning of VALUES block and reset `w' to original value
              do i=1, Nparam_check+1
                backspace(unit=iut)
              enddo
              !
              w = "VALUES"
              !
              Nparam = Nparam_check
              !
              if (Nparam <= 0) then
                  call report ("ERROR: Number of points or parameters <= 0 )",.true.)
              endif
              !
              field%Nterms = Nparam
              !
              ! Allocation of the pot. parameters
              !
              allocate(field%value(Nparam),field%forcename(Nparam),field%weight(Nparam),stat=alloc)
              call ArrayStart('field%value',alloc,size(field%value),kind(field%value))
              call ArrayStart('field%forcename',alloc,size(field%forcename),kind(field%forcename))
              call ArrayStart('field%weight',alloc,size(field%weight),kind(field%weight))
              !
              allocate(field%ipower1(Nparam),field%ipower2(Nparam),field%ipower3(Nparam),stat=alloc)
              call ArrayStart('field%ipower',alloc,size(field%ipower1),kind(field%ipower1))
              call ArrayStart('field%ipower',alloc,size(field%ipower2),kind(field%ipower2))
              call ArrayStart('field%ipower',alloc,size(field%ipower3),kind(field%ipower3))
              !
              allocate(field%iaddress(Nparam),stat=alloc)
              call ArrayStart('field%iaddress',alloc,size(field%ipower1),kind(field%ipower1))
              !
              allocate(field%link(Nparam),stat=alloc)
              call ArrayStart('field%link',alloc,size(field%link),8_ik)
              field%link(:)%ifield = 0
              field%link(:)%iparam = 1
              !
              field%value = 0
              field%forcename = 'dummy'
              field%weight = 0
              !
              iparam = 0
              !
              do while (trim(w)/="".and.iparam<Nparam.and.trim(w)/="END")
                 !
                 call read_line(eof,iut) ; if (eof) exit
                 !
                 iparam = iparam+1
                 !
                 if (nitems<2) then
                    !
                    write(out,"(a,i4)") "wrong number of records for an analytical field-type," // &
                                    "must be two at least (name value)",nitems
                    call report ("illegal number of records (<2) in the current field-line "//trim(w),.true.)
                    !
                 endif
                 !
                 !
                 call readu(field%forcename(iparam))
                 !
                 call readi(field%ipower1(iparam))
                 call readi(field%ipower2(iparam))
                 call readi(field%ipower3(iparam))
                 !
                 call readf(f_t)
                 !
                 field%value(iparam) = f_t
                 field%weight(iparam) = 0
                 !
                 if(nitems>5) then
                   !
                   call readu(w)
                   !
                   if (trim(w(1:1))=="F") then 
                     !
                     field%weight(iparam) = 1.0_rk
                     !
                   endif
                   !
                 endif 
                 !
                 if(trim(w(1:1))=="L".and.nitems>7) then
                   !
                   call readi(field%link(iparam)%ifield)
                   call readi(field%link(iparam)%iparam)
                   !
                   ! set the weight of the linked parameter to zero
                   !
                   field%weight(iparam) = 0
                   !
                 endif
                 !
                 total_parameters = total_parameters + 1
                 !
              enddo
              !
              call readu(w)
              !
              if (iparam/=Nparam.or.(trim(w)/="".and.trim(w)/="END")) then
                 !
                 print "(2a,2i6)","wrong number of rows in section: ",trim(field%name),iparam,Nparam
                 call report("illegal number of rows in a section: "//trim(w),.true.)
                 !
              endif
                 !
            case default
                 !
                 call report ("Unrecognized keyword (error 03): "//trim(w),.true.)
                 !
            end select
            !
            call read_line(eof,iut) ; if (eof) exit
            call readu(w)
          enddo
          !
       case("FITTING")
         !
         call read_line(eof,iut) ; if (eof) exit
         call readu(w)
         !
         do while (trim(w)/="".and.trim(w)/="END")
           !
           select case(w)
              !
           case('J_LIST','JLIST','J','JROT')
             !
             i = 0
             !
             do while (item<Nitems.and.i<jlist_max)
                !
                i = i + 1
                !
                call readu(w)
                !
                if (trim(w)/="-") then
                  !
                  read(w,*) jrot
                  !
                  j_list_(i) = jrot
                  !
                else
                  !
                  call readf(jrot2)
                  !
                  do while (item<=Nitems.and.nint(2.0*jrot)<nint(2.0*jrot2))
                    !
                    jrot = jrot + 1.0
                    j_list_(i) = jrot
                    i = i + 1
                    !
                  enddo
                  !
                  i = i - 1
                  !
                endif
                !
             enddo
             !
             jmax = maxval(J_list_)
             Nmax = nint(jmax-0.5)
             !
             nJ = i
             allocate(fitting%j_list(nJ),stat=alloc)
             !
             fitting%nJ = nJ
             fitting%J_list(1:i) = J_list_(1:i)
             !
           case('ITMAX','ITERMAX','ITER')
             !
             call readi(fitting%itermax)
             !
           case('ROBUST')
             !
             call readf(fitting%robust)
             !
           case('TARGET_RMS')
             !
             call readf(fitting%target_rms)
             !
           case('TARGET_STABILITY')
             !
             call readf(fitting%target_stability)
             !
           case('FIT_TYPE')
             !
             call readu(fitting%fit_type)
             !
           case('DERIVE_TYPE')
             !
             call readu(deriv_type)
             !
             select case (deriv_type)
             case ('HELLMANN-FEYNMAN','HELLMAN')
              deriv_type = 'HELLMANN'
             end select
             !
           case('THRESH_ASSIGN','THRESH_REASSIGN','THRESH_LOCK','LOCK','LOCK_QUANTA')
             !
             call readf(fitting%threshold_lock)
             !
           case('THRESH_OBS-CALC') 
             ! switch off weights for residuals larger than THRESH_OBS-CALC
             !
             call readf(fitting%threshold_obs_calc)
             !
           case("IPARAM")
             !
             call readi(fitting%iparam(1))
             call readi(fitting%iparam(2))
             !
           case('GEOMETRIES')
             !
             call readl(fitting%geom_file)
             !
           case('OUTPUT')
             !
             call reada(fitting%output_file)
             !
           case('ZPE')
             !
             call readf(fitting%zpe)
             fitting%shift_to_zpe = .false.
             !
           case('FIT_FACTOR')
             !
             call readf(fitting%factor)
             !
           case('FIT_SCALE')
             !
             call readf(fitting%fit_scaling)
             !
           case('LINEAR_SEARCH','LINEAR-SEARCH')
             !
             call readi(fitting%linear_search)
             !
           case('ABINIIO')
             !
             ! ignore the experiment and fit to the ab initio curves only
             !
             fitting%factor = small_
             !
           case('ENER','ENERGIES')
             !
             Nparam_check = 0
             !
             call input_options(echo_lines=.true.,error_flag=1)
             !
             do while (trim(w)/="END")
                !
                call read_line(eof,iut) ; if (eof) exit
                !
                call readu(w)
                !
                Nparam_check = Nparam_check+1
                !
             enddo
             !
             Nparam_check = Nparam_check-1
             !
             if (trim(w) /= "END") then
                 call report ("ERROR: Cannot find `END' statement)",.true.)
             endif
             !
             ! go back to beginning of VALUES block and reset `w' to original value
             do i=1, Nparam_check+1
               backspace(unit=iut)
             enddo
             !
             fitting%Nenergies = Nparam_check
             !
             !call readi(fitting%Nenergies)
             !
             allocate (fitting%obs(1:fitting%Nenergies),stat=alloc)
             !
             if (alloc/=0) then
               write (out,"(' Error ',i0,' initializing obs. energy related arrays')") alloc
               stop 'obs. energy arrays - alloc'
             end if
             !
             iobs = 0
             MXRoot = 1
             !
             call read_line(eof,iut) ; if (eof) exit
             call readu(w)
             !
             do while (trim(w)/="END".and.iobs<fitting%Nenergies)
                !
                iobs = iobs + 1
                !
                call reread(-1)
                !
                call readf(fitting%obs(iobs)%Jrot)
                !
                i = 0
                matchfound = .false.
                do while( i<fitting%nJ.and..not.matchfound )
                  !
                  i = i + 1
                  if (fitting%J_list(i)/=fitting%obs(iobs)%Jrot) cycle
                  matchfound = .true.
                  !
                enddo
                !
                ! skip current line if these Jrot-s are not processed
                if (.not.matchfound) then
                  iobs = iobs-1
                  call read_line(eof,iut) ; if (eof) exit
                  call readu(w)
                  cycle
                endif
                !
                fitting%obs(iobs)%irot = i
                !
                ! parity:
                !
                call readi(itau)
                !
                fitting%obs(iobs)%iparity = itau
                !
                call readi(fitting%obs(iobs)%N)
                call readf(fitting%obs(iobs)%energy)
                !
                call readi(fitting%obs(iobs)%quanta%istate)
                !
                ! skip current line if this state is outside the range 
                if (fitting%obs(iobs)%quanta%istate>Nestates) then 
                  iobs = iobs-1
                  call read_line(eof,iut) ; if (eof) exit
                  call readu(w)
                  cycle
                endif
                !
                ! iref = 1 (state=1) and iref=2 (state=2,3)
                !
                fitting%obs(iobs)%iref = min(fitting%obs(iobs)%quanta%istate,2)
                !
                call readi(fitting%obs(iobs)%quanta%v1)
                call readi(fitting%obs(iobs)%quanta%v2)
                call readi(fitting%obs(iobs)%quanta%v3)
                call readi(fitting%obs(iobs)%quanta%ilambda)
                call readf(fitting%obs(iobs)%quanta%omega)
                !
                call readf(fitting%obs(iobs)%weight)
                !
                MXRoot = max(MXRoot,fitting%obs(iobs)%N)
                !
                ! make sure the first line is the lowest state
                !
                if ( iobs==1.and.( nint(fitting%obs(iobs)%Jrot-0.5)/=0.or.itau/=1.or.fitting%obs(iobs)%N/=1.or.fitting%obs(iobs)%quanta%istate/=1) ) then
                  write(out,"('Illegal fitting input: 1st line must be the g.s. not J,tau,state,n = ',f8.1,3i,f12.3)") &
                            fitting%obs(iobs)%Jrot,itau,fitting%obs(iobs)%quanta%istate,fitting%obs(iobs)%N,fitting%obs(iobs)%energy
                  stop 'Illegal fitting input: 1st line must be the g.s.'
                endif
                !
                call read_line(eof,iut) ; if (eof) exit
                call readu(w)
                !
                ! make sure that the first entry corresponds to J=Jmin, iparity=1
                !
                if (iobs==1.and.(itau/=1.or.nint(fitting%obs(iobs)%Jrot-0.5_rk)/=0)) then 
                   write(f_out,"('Input error: first observed fitting entry must be J=0.5, parity=1',f8.1,i2)") & 
                                 fitting%obs(iobs)%Jrot,itau
                   stop 'Input error: first observed fitting entry must be J=0.5, parity=1'
                endif
                !
             enddo
             !
           case default
             !
             call report ("Unrecognized keyword name (error 01): "//trim(w),.true.)
             !
           end select
           !
           call read_line(eof,iut) ; if (eof) exit
           call readu(w)
           !
         enddo
         !
         fitting%Nenergies = iobs
         !
         if (trim(w)/="".and.trim(w)/="END") then
            call report ('wrong last line in FITTING ',.false.)
         endif

    case ('ABINITIO')
         !
         call read_line(eof,iut) ; if (eof) exit
         call readu(w)
         !
         do while (trim(w)/="".and.trim(w)/="END")
           !
           select case(w)
             !
           case('VMIN')
             !
             call readf(fitting%Vmin(1))
             call readf(fitting%Vmin(2))
             call readf(fitting%Vmin(3))
             call readf(fitting%Vmin(4))
             !
           case('VALUES')
             !
             Nparam_check = 0
             !
             call input_options(echo_lines=.true.,error_flag=1)
             !
             do while (trim(w)/="END")
                !
                call read_line(eof,iut) ; if (eof) exit
                !
                call readu(w)
                !
                Nparam_check = Nparam_check+1
                !
             enddo
             !
             Nparam_check = Nparam_check-1
             !
             if (trim(w) /= "END") then
                 call report ("ERROR: Cannot find `END' statement)",.true.)
             endif
             !
             ! go back to beginning of VALUES block and reset `w' to original value
             do i=1, Nparam_check+1
               backspace(unit=iut)
             enddo
             !
             fitting%Nabinitio = Nparam_check
             !
             !call readi(fitting%Nenergies)
             !
             allocate (fitting%ai(1:fitting%Nabinitio),fitting%weight(1:fitting%Nabinitio),stat=alloc)
             call ArrayStart('fitting%ai',alloc,size(fitting%ai),kind(fitting%ai))
             call ArrayStart('fitting%weight',alloc,size(fitting%weight),kind(fitting%weight))
             !
             allocate (fitting%r1(1:fitting%Nabinitio),fitting%r2(1:fitting%Nabinitio),fitting%r3(1:fitting%Nabinitio),stat=alloc)
             call ArrayStart('fitting%r1-2-3',alloc,size(fitting%r1),kind(fitting%r1))
             call ArrayStart('fitting%r1-2-3',alloc,size(fitting%r2),kind(fitting%r2))
             call ArrayStart('fitting%r1-2-3',alloc,size(fitting%r3),kind(fitting%r3))
             !
             allocate (fitting%iPES(1:fitting%Nabinitio),stat=alloc)
             call ArrayStart('fitting%iPES',alloc,size(fitting%iPES),kind(fitting%iPES))
             !
             iai = 0
             !
             call read_line(eof,iut) ; if (eof) exit
             call readu(w)
             !
             do while (trim(w)/="END".and.iai<fitting%Nabinitio)
                !
                iai = iai + 1
                !
                call reread(-1)
                call readf(fitting%r1(iai))
                call readf(fitting%r2(iai))
                call readf(fitting%r3(iai))
                call readf(fitting%ai(iai))
                call readf(fitting%weight(iai))
                call readi(fitting%iPES(iai))
                !
                ! subtract the corresponding minim values 
                !
                fitting%ai(iai) = fitting%ai(iai) - fitting%Vmin(fitting%iPES(iai))
                !
                call read_line(eof,iut) ; if (eof) exit
                call readu(w)
                !
             enddo
             !
           case default
             !
             call report ("Unrecognized keyword name (error 01): "//trim(w),.true.)
             !
           end select
           !
           call read_line(eof,iut) ; if (eof) exit
           call readu(w)
           !
         enddo
         !
         !fitting%Nenergies = iobs
         !
         if (trim(w)/="".and.trim(w)/="END") then
            call report ('wrong last line in ABINIIO ',.false.)
         endif
         !
    case default
        !
        call report ("Unrecognized unit name "//trim(w),.true.)
        !
    end select
    !
   end do
   !
   !
   write(f_out,"('Number of obs. data points: ',i9)") fitting%Nenergies
   !
   ! Allocate the obs. energy related arrays
   !
   !allocate (ener_obs(en_npts),enercalc(en_npts),J_obs(en_npts),&
   !          sym_obs(en_npts),gamma_obs(en_npts),N_obs(en_npts),&
   !          quanta_obs(3,en_npts),stat=alloc)
   !if (alloc/=0) then
   !  write (f_out,"(' Error ',i,' initializing obs. energy related arrays')") alloc
   !  stop 'obs. energy arrays - alloc'
   !end if
   !
   ! Potential energy points 
   !
   ! Count them:
   !
   pot_npts = fitting%Nabinitio
   !
   write(f_out,"('Number of potential energy data points: ',i9)") pot_npts
   !
   npts = pot_npts + fitting%Nenergies
   !
   allocate (wtall(npts),wt_bit(npts),stat=alloc)
   !
   call ArrayStart('wtall',alloc,size(wtall),kind(wtall))
   call ArrayStart('wt_bit',alloc,size(wt_bit),kind(wt_bit))
   !
   allocate(Nstates(2,0:Nmax,2),stat=alloc)
   call ArrayStart('Nstates',alloc,size(Nstates),kind(Nstates))
   !
   allocate(calc(2,0:Nmax,2),stat=alloc)
   call ArrayStart('calc',alloc,1,4)
   !
   if (fitting%robust>0) then 
     allocate (sigma(npts),stat=alloc)
     call ArrayStart('sigma',alloc,size(sigma),kind(sigma))
   endif
   !
   wtall(1:fitting%Nenergies) = fitting%obs(1:fitting%Nenergies)%weight
   wtall(fitting%Nenergies+1:npts) = fitting%weight(1:pot_npts)
   !
   wtsum = sum(wtall(1:fitting%Nenergies))
   wtall(1:fitting%Nenergies) = wtall(1:fitting%Nenergies)/wtsum
   !
   wtsum = sum(wtall(fitting%Nenergies+1:npts))
   wtall(fitting%Nenergies+1:npts) = wtall(fitting%Nenergies+1:npts)/wtsum
   !
   ! 2. Factorizing the obs. weights by the factor "fit_factor":
   !
   wtall(1:fitting%Nenergies) = wtall(1:fitting%Nenergies)*fitting%factor
   !
   ! 3. And normilizing the all weight factors. 
   !
   wtsum = sum(wtall(1:npts))
   wtall(1:npts) = wtall(1:npts)/wtsum
   ! 
   ! Count how many data points actually will be fitted. 
   !
   nused=0
   wtsum=0
   !
   wt_bit = 0 
   !
   do i=1,npts
     if (wtall(i) .gt. 0.0_rk) then 
       nused=nused+1
       !
       wt_bit(i) = 1.0d0
       !
     endif
   enddo
   !
   !
   !     For the details on the robust fit see J.K.G. Watson, JMS, 219, 326 (2003).
   !     The idea is optimize the weigth factors iteratevely depending on the fitting results, as  follows:
   !     w_i = 1 / sigma_i^2+alpha*eps_i^2 
   !     where sigma_i represent the experimental data precision, alpha is the Watson parameter (can be also adjsuted) and 
   !     eps_i = Obs - Calc deviatons. 
   !
   !                   robust_fit = 0 corresponds to the standard fit, the robust fittting procedure is switched off. 
   !                   robust_fit <> =0, the robust fittting procedure is on.
   !                   robust_fit also defines the target accuracy of the fit. For example for robust_fit  = 0.001 
   !                   all data with (obs-calc) < 0.001 will be considered as relatively bad (i.e. outliers ) 
   ! sigma = exp. data precision for robust fit 
   !
   if (fitting%robust>0) then
     !
     sigma = 1.0_rk
     do i=1,npts
       if (wtall(i)>0) sigma(i) = sigma(i)/sqrt(wtall(i))*fitting%robust
     enddo
     !
     !wtsum = 1.0d0 ! sqrt(sum(sigma(1:fitting%Nenergies)**2))
     !sigma(:) = sigma(:)*fitting%robust/wtsum
     !
   endif 
   !
   write(f_out,"('Number of data points used in the fit: ',i9)") nused
   !
   ! Allocate objects, that will be used for the fitting procedure:
   !
   allocate (rjacob(npts,total_parameters),eps(npts),stat=alloc)
   call ArrayStart('rjacob',alloc,size(rjacob),kind(rjacob))
   call ArrayStart('eps',alloc,size(eps),kind(eps))
   !
   allocate (potparam(total_parameters),stat=alloc)
   call ArrayStart('potparam',alloc,size(potparam),kind(potparam))
   allocate (parold(total_parameters),stat=alloc)
   call ArrayStart('parold',alloc,size(parold),kind(parold))
   allocate (ivar(total_parameters),stat=alloc)
   call ArrayStart('ivar',alloc,size(ivar),kind(ivar))
   allocate (nampar(total_parameters),stat=alloc)
   call ArrayStart('nampar',alloc,size(nampar),kind(nampar))
   allocate (ifitparam(total_parameters),stat=alloc)
   call ArrayStart('ifitparam',alloc,size(ifitparam),kind(ifitparam))
   allocate (iparamfit(total_parameters),stat=alloc)
   call ArrayStart('iparamfit',alloc,size(iparamfit),kind(iparamfit))
   !
   allocate (ener_obs(fitting%Nenergies),enercalc(fitting%Nenergies),stat=alloc)
   call ArrayStart('ener_obs',alloc,size(ener_obs),kind(ener_obs))
   call ArrayStart('enercalc',alloc,size(enercalc),kind(enercalc))
   !
   allocate (mark(fitting%Nenergies),stat=alloc)
   if (alloc/=0) stop 'error allocating mark' 
   !
   allocate (al(total_parameters,total_parameters),ai(total_parameters,total_parameters),bl(total_parameters),&
             dx(total_parameters),sterr(total_parameters),Tsing(total_parameters,total_parameters),stat=alloc)
   call ArrayStart('al',alloc,size(al),kind(al))
   call ArrayStart('ai',alloc,size(ai),kind(ai))
   call ArrayStart('bl',alloc,size(bl),kind(bl))
   call ArrayStart('dx',alloc,size(dx),kind(dx))
   call ArrayStart('sterr',alloc,size(sterr),kind(sterr))
   call ArrayStart('Tsing',alloc,size(Tsing),kind(Tsing))
   !
   ! build the mapping betweeb potparam and fields 
   !
   iaddress = 0 
   !
   do ifield =1,Nfields
     !
     do iterm = 1,poten(ifield)%Nterms
       !
       iaddress = iaddress + 1
       !
       poten(ifield)%iaddress(iterm) = iaddress
       !
     enddo
   enddo
   !
   call map_parameters(dir=.true.)
   !   
   parold=potparam
   !
   mark(:) = ' '
   !
   ! The last object to allocate - the lapack related work array
   !
   lwork = 50*total_parameters
   !
   allocate (wspace(lwork),stat=alloc)
   call ArrayStart('wspace',alloc,size(wspace),kind(wspace))
   !
   ! prepare the file to write all computed energies 
   !
   open(f_en,file=trim(out_name)//'.en',status='replace')
   !
   ! fititer will count the iterations. 
   !
   fititer = 0
   !
   ! The initial parameters to be remembered. 
   !
   ! Create the input file for the expectaion values jobs, 
   ! to be used in connaction with the program 'xpext3.x'
   !
   !call create_the_xpect_job_file(total_parameters,ivar)
   !
   ! The outer fitting loop - allows to restart the fitting in case 
   ! we decide to remove some of the varying parameters. 
   ! This option is working together with fit_type ='linur'
   !
   still_run = .true.
   outer_loop: do while (still_run)  
      !
      ! Initial values for the standard error and  stability.
      !
      stadev_old = 1e10
      stability = 1e10
      stadev    = 1e10
      !
      ! Count the actually varying number of parameters:
      !
      numpar  = 0
      ifitparam = 1
      iparamfit = 0
      do i=1,total_parameters
        if (ivar(i) > 0) then 
          numpar=numpar+1
          ifitparam(numpar) = i
          iparamfit(i) = numpar
        endif
      enddo 
      !
      rjacob = 0 
      !
      Total_size = 0
      !
      ! Prepaper the first run:
      !
      !!!open(1,file='evib.out',status='replace') ; write(1,"('Start fitting...')") ; close(1)
      !!!open(1,file='erot.out',status='replace') ; write(1,"('Start fitting...')") ; close(1)
      !!!open(1,file='xpect.out',status='replace')  ; write(1,"('Start fitting...')") ; close(1)
      !
      !isys = systemqq('rm evib.out erot.out xpect.out')
      !
      ! The loop starts here. 
      !
      do while( fititer.le.fitting%itermax .and. stadev.ge.fitting%target_rms.and. stability.ge.fitting%target_stability)
        !
        fititer = fititer + 1
        !
        Nstates = 0 
        !
        ! If fitting%factor is set zero - there would be no fit to the experimental energies. 
        !
        if (fitting%factor<1e-12) rjacob(1:fitting%Nenergies,:) = 0
        !
        write(f_out,"(/'Iteration = ',i8)") fititer
        !
        ! create files with potential parameters 
        !
        call map_parameters(dir=.false.)
        !
        call write_input_potential_parameters(Nfields,poten)
        !
        if (verbose>=3) write(f_out,"(/'calling evvib Everest program')")
        !
        write(kmax_ch,"(i3)") nint(kmax+0.5)
        write(jmax_ch,"(i3)") nint(2.0_rk*jmax)
        write(fititer_ch_,"(i2)") fititer
        write(fititer_ch,"(a2)") adjustl(fititer_ch_)
        !
        ! generate the vib-input and rot-input files using the tempates provided
        !
#if (debug_ == 0)
         !
         if (verbose>=4) then 
             write(f_out,"('Make sure that the folliwng files are present in the folder:')")
             write(f_out,"('  - vib_inp.sh')")
             write(f_out,"('  - rot_inp.sh')")
             write(f_out,"('  - ',a)") trim(evib_exe)
             write(f_out,"('  - ',a)") trim(erot_exe)
             write(f_out,"('  - ',a)") trim(dervib_exe)
             write(f_out,"('  - ',a)") trim(derrot_exe)
             write(f_out,"(a,a,' and ',a)") '  - user defined potential function executables in ',trim(evib_exe),trim(erot_exe)
           endif
           !
           if (verbose>=5.and.fititer>1) then 
               isys = systemqq('cp evib.inp evib.inp'//fititer_ch)
               isys = systemqq('cp erot.inp erot.inp'//fititer_ch)
               isys = systemqq('cp erot.out erot.out'//fititer_ch)
           endif
           !
           isys = systemqq('vib_inp.sh '//kmax_ch)
           !
           isys = systemqq('rot_inp.sh '//kmax_ch//' '//jmax_ch)
           !
#endif
        !
        if (verbose>=3) write(f_out,"('calling vibrational ',a,' program, kmax = ',f8.1)") trim(evib_exe),kmax
        !
#if (debug_ == 0)
          isys = systemqq(evib_exe//'< evib.inp > evib.out')
#endif
        !
        if (verbose>=3) write(f_out,"('calling  rovibronic ',a,' program, kmax, jmax= ',2f8.1)") trim(erot_exe),kmax,jmax
#if (debug_ == 0)
          isys = systemqq(erot_exe//'<erot.inp >erot.out')
#endif
        !
        if (verbose>=3) write(f_out,"('collecting the rovibronic energies ...')") 
        !
#if (debug_ == 0)
           isys = systemqq('./collect_energies_rot.sh erot.out > erot.log')
#endif
        !
        ! initial counting of states
        !
        if (verbose>=4) write(f_out,"('   Counting states...')")
        !
        call get_vibronic_energies(Nmax,Nstates,calc)
        !
        ! check if the number of the states changed and the energy arrays have wrong  size and need to be reallocated 
        !
        if (associated(calc(1,0,1)%energy).and.sum(Nstates)/=Total_size) then 
           !
           do iref  = 1,2
             do irot = 0,Nmax
               do ipar = 1,2
                 !
                 deallocate(calc(iref,irot,ipar)%energy,calc(iref,irot,ipar)%quanta)
                 !
               enddo
             enddo 
           enddo
           call ArrayStop('calc%energy')
           call ArrayStop('calc%quanta')
        endif
        !
        Total_size = sum(Nstates)
        !
        ! allocate energy arrays before reading the states if necessary 
        !
        if (.not.associated(calc(1,0,1)%energy)) then
          !
          do iref  = 1,2
            do irot = 0,Nmax
              do ipar = 1,2
                !
                Nsize = max(1,Nstates(iref,irot,ipar))
                !
                allocate(calc(iref,irot,ipar)%energy(Nsize),calc(iref,irot,ipar)%quanta(Nsize),stat=alloc)
                call ArrayStart('calc%energy',alloc,size(calc(iref,irot,ipar)%energy),kind(calc(iref,irot,ipar)%energy))
                call ArrayStart('calc%quanta',alloc,size(calc(iref,irot,ipar)%quanta),80_ik)
                !
              enddo
            enddo 
          enddo
        endif
        !
        if (verbose>=4) write(f_out,"('   Collecting energies...')")
        !
        call get_vibronic_energies(Nmax,Nstates,calc)
        !
        ! Zero point energy:
        !
        ezero_ = calc(1,0,1)%energy(1)
        !
        if (verbose>=3) write(f_out,"('   ZPE = ',f15.6)") ezero_
        !
        ! if threshold_lock>0, correct the addresses of the energies in case of accidential swaps
        ! by comparing with energies within the "threshold_lock"-range and with "exp" quantum numbers.
        ! if threshold_lock<0, find closets matches with experimental energies 
        ! if threshold_lock=0, do nothing
        !
        mark = " "
        !
        if (abs(fitting%threshold_lock)>0) then 
          !
          k0 = 1
          !
          if (fitting%threshold_lock<0)  maxiter_as = 1
          !
          jrot_ = fitting%obs(1)%jrot
          itau_ = 1 ; if ( fitting%obs(1)%iparity==-1 ) itau_ = 2
          !
          do iobs = 1,fitting%Nenergies
            !
            Jrot = fitting%obs(iobs)%Jrot
            iparity = fitting%obs(iobs)%iparity
            iref = fitting%obs(iobs)%iref
            irot = nint(Jrot-0.5_rk)
            !
            itau = 1 ; if (iparity==-1) itau = 2
            !
            if (Jrot>jmax) cycle
            !
            iener = fitting%obs(iobs)%N
            !
            mark(iobs) = '*'
            !
            loop_thresh : do iter_th = 1,maxiter_as
               !
               loop_iener : do  iener_= k0,Nstates(iref,irot,itau)
                    !
                    do jobs = max(iobs-iener_,1),iobs-1
                      !
                      if (nint(jrot-fitting%obs(jobs)%jrot)/=0.or.itau-1/=fitting%obs(jobs)%iparity) cycle

                      if (nint(Jrot-fitting%obs(jobs)%Jrot)/=0.or.iref/=fitting%obs(jobs)%iref.or.&
                         iparity/=fitting%obs(jobs)%iparity ) cycle
                      !
                      if ( fitting%obs(jobs)%N==iener ) cycle loop_iener
                      !
                    enddo
                    !
                    if (iener_>Nstates(iref,irot,itau)) then
                      write(out,"('fitting error: the size of calc energy is too small ',i8,';')") Nstates(iref,irot,itau)
                      write(out,"('               J= ',f9.2,' and state  ',i8)") jrot,iener_
                      stop 'error: the size of calc energy is too small'
                    endif
                    !
                    ! Count MaxRoots for the Hellmann-Feynman derivatives in Everest 
                    !
                    MXRoot = max(MXRoot,iener_)
                    !
                    if ( abs( fitting%obs(iobs)%energy-( calc(iref,irot,itau)%energy(iener_)-ezero_ ) )<= &
                             real(iter_th,rk)*abs(fitting%threshold_lock).and.&
                          ! threshold_lock>0, QN match
                        ( ( calc(iref,irot,itau)%quanta(iener_)%v1==fitting%obs(iobs)%quanta%v1.and.&
                            calc(iref,irot,itau)%quanta(iener_)%v2==fitting%obs(iobs)%quanta%v2.and.&
                            calc(iref,irot,itau)%quanta(iener_)%v3==fitting%obs(iobs)%quanta%v3.and.&
                            calc(iref,irot,itau)%quanta(iener_)%ilambda==fitting%obs(iobs)%quanta%ilambda.and.&
                            nint(calc(iref,irot,itau)%quanta(iener_)%omega-fitting%obs(iobs)%quanta%omega)==0.and.&
                            calc(iref,irot,itau)%quanta(iener_)%istate ==fitting%obs(iobs)%quanta%istate ).or.& 
                           ! energy match + state + vib QN match 
                           ( fitting%threshold_lock<0.and. &
                            calc(iref,irot,itau)%quanta(iener_)%istate ==fitting%obs(iobs)%quanta%istate ) ) ) then 
                        !
                        fitting%obs(iobs)%N = iener_
                        mark(iobs) = ' '
                        exit loop_thresh
                        !
                    endif 
                    !
                 enddo loop_iener
                 !
            enddo loop_thresh
            !
            if (mark(iobs) == ' ')  k0 = max(fitting%obs(iobs)%N-iener_,1)
            !
          enddo
          !
        endif
        !
        !
        ! The derivatives of the energy wrt parameters will be evalueted 
        ! only if 1) it is a fitting job, not just energy calculations (j/=0) and
        !         2) it is a energy fitting job, not just a fit to the ab initio points 
        !            (fitting%factor is not zero). 
        !
        if (fitting%itermax.ge.1.and.fitting%factor>1e-12.and.mod(fititer+2,3)==0) then 
          !
          ! calculating the derivatives 
          ! differencies. It is essentially slower and we use it only for 
          ! the testing of the xpect3 derivativies.
          !
          if (trim(deriv_type)/='HELLMANN'.and.fitting%itermax.ge.1.and.fitting%factor>1e-12) then
            !
            if (verbose>=3) write(f_out,"('Using finete differences for the derivatives')")
            !
            call finite_diff_of_energy(rjacob,Nstates,potparam)
            !
          else
            !
#if (debug_ == 0)
            !
            ! Everst derivatives using Helmann-Feynman method
            !
            if (verbose>=3) write(f_out,"('   Calling the vib. derivs ',a,' program ...')") trim(dervib_exe)
            !
            isys = systemqq(dervib_exe//'< vde.inp > vde.out')
            !
            if (verbose>=3) write(f_out,"('   Calling the rovib. derivs Helmann-Feynman ',a,' program ...')") trim(derrot_exe)
            !
            isys = systemqq(derrot_exe//'< rde.inp > rde.out')
            !
#endif
            !
            ! Create the input files for the Hellmann-Feyman derivatives
            !
            if (trim(deriv_type)=='HELLMANN') then 
              !
              call create_the_deriv_job_file(MXRoot,iOMP)
              !
            endif
            !
            if (verbose>=4) write(f_out,"('   Collecting derivatives...')")
            !
            call get_vibronic_derivatives(Nmax,total_parameters,calc,iparamfit,rjacob)
            !
            ! Jacobi matrix stores derivatives only for energies that have obs. counterparts 
            ! in input with weight /= 0, i.e. participating in the fit. 
            !
          endif 
          !
        endif 
        !
        ! Print out the calc. and obs.-calc., i.e. result of the fit. 
        ! Only fitted energies are printed. 
        !
        write(f_out,"(/1x,100('-'))")
        write(f_out,"('| ## |  N |      J  |par |    Obs.    |   Calc.    | Obs.-Calc. |   Weight |    St v1 v2 v3 La Ome  |')")
        write(f_out,"(1x,100('-'))")
        !
        eps = 0
        !
        do nrow = 1,fitting%Nenergies
          !
          Jrot = fitting%obs(nrow)%Jrot
          iparity = fitting%obs(nrow)%iparity
          iref = fitting%obs(nrow)%iref
          Nrot = nint(Jrot-0.5_rk)
          !
          ipar_ = 1 ; if (iparity==-1) ipar_ = 2
          !
          if (Jrot>jmax) cycle
          !
          iener = fitting%obs(nrow)%N
          !
          if (iener>Nstates(iref,Nrot,ipar_)) cycle 
          !enercalc(nrow) = 0
          !
          enercalc(nrow) = calc(iref,Nrot,ipar_)%energy(iener)-ezero_
          !
          eps(nrow) = fitting%obs(nrow)%energy-enercalc(nrow)
          !
          write (f_out,"(2i5,f8.1,i5,2x,3f13.4,2x,e9.2,4x,5(i3),f5.1,2x,5(i3),f5.1,a)") &
             nrow,iener,Jrot,iparity,enercalc(nrow)+eps(nrow),enercalc(nrow),-eps(nrow),&
             wtall(nrow),&
             calc(iref,Nrot,ipar_)%quanta(iener)%istate,&
             calc(iref,Nrot,ipar_)%quanta(iener)%v1,&
             calc(iref,Nrot,ipar_)%quanta(iener)%v2,&
             calc(iref,Nrot,ipar_)%quanta(iener)%v3,&
             calc(iref,Nrot,ipar_)%quanta(iener)%ilambda,&             
             calc(iref,Nrot,ipar_)%quanta(iener)%omega,&
             fitting%obs(nrow)%quanta%istate,&
             fitting%obs(nrow)%quanta%v1,&
             fitting%obs(nrow)%quanta%v2,&
             fitting%obs(nrow)%quanta%v3,&
             fitting%obs(nrow)%quanta%ilambda,&
             fitting%obs(nrow)%quanta%omega,&
             mark(nrow)
             !
        enddo ! --- nrow
        !
        ! Printing all calculated term values. If the obs. counterpats exist, 
        ! the obs.-calc. are printed as well. 
        ! This list can be used to identify the obs-calc pairs, i.e. the number of 
        ! a term value as it appear in the list of calculated energies. 
        !
        write(f_en,"(/'Iteration = ',i4)") fititer
        !
        write(f_en,"(/1x,100('-'))")
        write(f_en,"('| ## |  N |      J  |par |    Obs.    |   Calc.    | Obs.-Calc. |   Weight |    St v1 v2 v3 La Ome  |')")
        write(f_en,"(1x,100('-'))")
        !
        do iref = 1,2
           !
           do Nrot = 0,Nmax
              !
              Jrot = real(Nrot,rk)+0.5_rk
              !
              do ipar_ = 1,2
                !
                !calc(Nrot,ipar)%energy(iener)
                !
                iparity = 1 ; if (ipar_==2) iparity = -1
                !
                do irow = 1,Nstates(iref,Nrot,ipar_)
                  !
                  iobs = 0
                  !
                  loop_iobs : do 
                    !
                    iobs = iobs+1
                    if (irow==fitting%obs(iobs)%N.and.nint(Jrot-fitting%obs(iobs)%Jrot)==0 &
                        .and.iref==fitting%obs(iobs)%iref &
                        .and.iparity==fitting%obs(iobs)%iparity ) then
                      !
                      exit loop_iobs
                      !
                    endif
                    !
                    if (iobs==fitting%Nenergies) exit loop_iobs 
                    !
                  enddo loop_iobs
                  !
                  if (iobs<fitting%Nenergies) then
                     !
                     mark = " "
                     !
                     if (abs(eps(iobs))>1.0) mark = "!"
                     !
                     write(f_en,"(2i5,f8.1,i5,2x,3f13.4,2x,e9.2,4x,5(i3),f5.1,2x,5(i3),f5.1,a)") &
                                     irow,fitting%obs(iobs)%n,Jrot,iparity,&
                                     !fitting%obs(iobs)%iref,&
                                     enercalc(iobs)+eps(iobs),enercalc(iobs),-eps(iobs),&
                                     wtall(iobs),&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%istate,&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%v1,&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%v2,&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%v3,&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%ilambda,&                                     
                                     calc(iref,Nrot,ipar_)%quanta(irow)%omega,&
                                     fitting%obs(iobs)%quanta%istate,&
                                     fitting%obs(iobs)%quanta%v1,&
                                     fitting%obs(iobs)%quanta%v2,&
                                     fitting%obs(iobs)%quanta%v3,&
                                     fitting%obs(iobs)%quanta%ilambda,&
                                     fitting%obs(iobs)%quanta%omega,&
                                     mark(iobs)

                     !
                  else
                     !
                     write(f_en,"(2i5,f8.1,i5,2x,3f13.4,2x,e9.2,4x,5(i3),f5.1)") irow,0,Jrot,iparity,0.0_rk,&
                                     calc(iref,Nrot,ipar_)%energy(irow)-ezero_,0.0_rk,0.0_rk,&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%istate,&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%v1,&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%v2,&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%v3,&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%ilambda,&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%omega

                                     
                     !
                  endif
                  !
                enddo
                !
              enddo
              !
           enddo
           !
        enddo 
        !
        !
        ! switch off energies with large obs-calc assuming unwanted swapping 
        if (fitting%threshold_obs_calc>small_) then
           do iener = 1,fitting%Nenergies
              if (abs(eps(iener))>fitting%threshold_obs_calc.and.wtall(iener)>small_) then
                wtall(iener) = 0
                if (fitting%robust>small_) sigma(iener) = 0
                eps(iener) = 0
                nused = nused - 1
              endif
           enddo
           nused = max(nused,1)
           wtsum = sum(wtall(1:npts))
           wtall(1:npts) = wtall(1:npts)/wtsum
        endif 
        !
        ! Here the potential energy section starts. 
        !
        do nrow=1,pot_npts
          !
          !  we assume here that the angles are written in degrees and bond length in in Angstrom
          !     
          r1 =fitting%r1(nrow)
          r2 =fitting%r2(nrow)
          theta=fitting%r3(nrow)*pi/180.0_rk
          !
          ifield = fitting%iPES(nrow)
          !
          if (ifield>Nfields) then 
            stop 'Potential section: istate is outside the range'
          endif
          !
          ! Call the potential function. It has to be included into the "poten" subroutine
          !
          call poten_func(poten(ifield),r1,r2,theta,v)
          !
          ! eps = ab initio energies - calculated pot. energies,
          ! where we comntinue counting the fitting data points starting with fitting%Nenergies - 
          ! obs. data. 
          !
          eps(nrow+fitting%Nenergies) = fitting%ai(nrow)-v
          !
          ! Calculate derivatives with respect to parameters 
          ! using the finite diff. method. 
          !
          if (fitting%itermax.ge.1.and.(mod(fititer-1,nofititer+1).eq.0)) then
            !
            do  ncol=1,numpar
              !
              i = ifitparam(ncol)
              !
              ! parameter from wrong states can be skipped 
              !
              call which_field_address(i,ifield_,iaddr)
              !
              if (ifield/=ifield_) cycle
              tempx=potparam(i)
              deltax=fitfactordeltax*abs(tempx)
              if (deltax .le. 1e-15) deltax=1e-5
              !
              potparam(i)=tempx+deltax
              poten(ifield)%value(iaddr) = potparam(i)
              !
              ! update linked values if necessary 
              call update_linked_parameters
              !
              call poten_func(poten(ifield),r1,r2,theta,potright)
              !
              potparam(i)=tempx-deltax
              !
              ! update linked values if necessary 
              call update_linked_parameters
              !
              poten(ifield)%value(iaddr) = potparam(i)
              call poten_func(poten(ifield),r1,r2,theta,potleft)
              !
              potparam(i)=tempx
              poten(ifield)%value(iaddr) = potparam(i)
              !
              ! update linked values if necessary 
              call update_linked_parameters
              !
              rjacob(nrow+fitting%Nenergies,ncol)=(potright-potleft)/(2.0_rk*deltax)
            enddo ! --- ncol
            !
          endif     
          !
        enddo  ! ---  nrow
        !
        ! ssq  - weighted rms**2, rms  - root mean square deviation. 
        !
        ssq=sum(eps(1:npts)*eps(1:npts)*wtall(1:npts))
        rms=sqrt(sum(eps(1:npts)*eps(1:npts))/npts)
        !
        ! Prepare the linear system a x = b as in the Newton fitting approach.  
        !
        if (fitting%itermax.ge.1) then
           !----- form the a and b matrix ------c
           ! form A matrix 
           do irow=1,numpar       
             do icolumn=1,irow    
               al(irow,icolumn)=sum(rjacob(1:npts,icolumn)*rjacob(1:npts,irow)*wtall(1:npts))
               al(icolumn,irow)=al(irow,icolumn)
             enddo
           enddo
           !
           ! form B matrix 
           do irow=1,numpar      
             bl(irow)=sum(eps(1:npts)*rjacob(1:npts,irow)*wtall(1:npts))
           enddo   
           !
           ! Two types of the linear solver are availible: 
           ! 1. linur (integrated into the program, from Ulenikov Oleg)
           ! 2. dgelss - Lapack routine (recommended).
           !
           select case (trim(fitting%fit_type)) 
           !
           case default
             !
             write (f_out,"('fit_type ',a,' unknown')") trim(fitting%fit_type)
             stop 'fit_type unknown'
             !
           case('LINUR') 
             !
             call linur(numpar,numpar,al(1:numpar,1:numpar),bl(1:numpar),dx(1:numpar),ierror)
             !
             ! In case of dependent parameters  "linur" reports an error = ierror, 
             ! which is a number of the dependent parameter. We can remove this paramter 
             ! from the fit and set its value to zero. And start the iteration again. 
             !
             if (ierror.ne.0) then 
               ncol=0
               do i=1,total_parameters
                  if (ivar(i) .ne. 0) then
                    ncol=ncol+1
                    if  ( ncol.eq.ierror ) then 
                        ivar(i) = 0
                        potparam(i) = parold(i)
                        write(f_out,"(i,'-th is out - ',a8)") i,nampar(i)
                    endif 
                  endif 
               enddo 
               cycle outer_loop    
              endif 
              !
           case ('DGELSS')
             !
             ai = al 
             call dgelss(numpar,numpar,1,ai(1:numpar,1:numpar),numpar,bl(1:numpar),numpar,Tsing,1.D-12,RANK0,wspace,lwork,info_)
             !
             if (info_/=0) then
               write(f_out,"('dgelss:error',i)") info_
               stop 'dgelss'
             endif
             !
             dx = bl
             !
           end select 
           !
           !----- update the parameter values ------!
           !
           do  ncol=1,numpar
              i = ifitparam(ncol)
              potparam(i)=potparam(i)+dx(ncol)
           enddo
           !
           if (fitting%robust>0) then
             !
             call robust_fitting(sigma(1:npts),eps(1:npts),wtall(1:npts))
             !
             ssq=sum(eps(1:npts)*eps(1:npts)*wtall(1:npts))
             !
           endif 
          !
           ! write the potential parameters into the file, where  
           ! we store the current values 
           !
           call map_parameters(dir=.false.)
           call write_input_potential_parameters(Nfields,poten)
           !
           ! Estimate standard deviation error. 
           !
           if ( nused.ne.numpar ) then 
             stadev=dsqrt(ssq/float(nused-numpar))
           else 
             stadev=dsqrt(ssq/nused)
           endif
           !
           ! Estimate the standard errors for each parameter using 
           ! the inverse matrix of a. 
           !
           call invmat(al,ai,numpar,total_parameters)
           !
           sum_sterr=0.d0
           do  ncol=1,numpar
              i = ifitparam(ncol)
              if (nused.eq.numpar) then  
                 sterr(ncol)=0
              else
                 sterr(ncol)=sqrt(abs(ai(ncol,ncol)))*stadev
                 sum_sterr=sum_sterr+abs(sterr(ncol)/potparam(i))
              endif
           enddo    
           !
           sum_sterr=sum_sterr/numpar 
           !
           ! This is how we define stability of the fit:
           ! as a relative change of stadev comparing with the step before. 
           !
           stability=abs( (stadev-stadev_old)/stadev )
           stadev_old=stadev
           !
        else
           !
           stadev=dsqrt(ssq/nused)
           !
        endif
        !
        ! Print the updated parameters. 
        !
        write(f_out,"(/'Potential paramters:')")
        !
        call print_parameters
        !
        !do i=1,total_parameters
        !  write (f_out,"(a8,4x,i2,e22.14)") nampar(i),ivar(i),potparam(i)
        !enddo
        !
        ! Print the potential energy points into a separate unit. 
        !
        inquire(f_potpoints,opened=ifopen)
        if ( ifopen ) then
           rewind(f_potpoints)
        else
           open (f_potpoints,file=trim(out_name)//'.pot',status='replace' )
        endif
        !
        write(f_potpoints,"(1h1,5x,12('*'),' ab initio points ',  &
             12('*')// &
             4x,'r1',5x,'r2',5x,'theta',7x,'ab initio PES',3x, &
             'cal.PES',3x,'a-c',3x,'weight'/)")
        !
        do nrow=1,pot_npts
          !
          r1 =fitting%r1(nrow)
          r2 =fitting%r2(nrow)
          theta=fitting%r3(nrow)*pi/180.0_rk
          !
          istate = fitting%iPES(nrow)
          !
          v = fitting%ai(nrow)-eps(nrow+fitting%Nenergies)
          write (f_potpoints,"(2f11.4,1x,f12.5,1x,i3,f12.4,2f12.4,e10.2)") r1,r2,theta, &
                 istate,fitting%ai(nrow),v, &
                 eps(nrow+fitting%Nenergies),wtall(nrow+fitting%Nenergies) 
          !
        enddo
        !
        ! Output some staqtistics and results 
        !
        !  only if we are fitting:  
        !
        if (fitting%itermax.ne.0) then
          !
          !still_run = .false.
          ! 
          write (f_out,"(//'Potential parameters rounded in accord. with their standart errors'/)")
          l = 0 
          do i=1,total_parameters
            if (ivar(i) .ne. 0) then
               l=l+1
               ndigits = 0
               conf_int = sterr(l)
               do while (conf_int.le.10.0)
                 ndigits = ndigits +1 
                 conf_int = conf_int*10
               enddo
               !
               !if (ndigits<1) ndigits = 12
               !
               if (conf_int>1e8) conf_int = 0 
               !
               write(fmt0,"(i2)") ndigits
               fmt0 = adjustl(fmt0)
               fmt1 = '(a8,i4,2x,f22.'//fmt0//'    ''    ('',i14,'')'')'
               write (f_out,FMT=fmt1) nampar(i),ivar(i),potparam(i),nint(conf_int)
            else 
               ndigits =2
               if (potparam(i).ne.0.0) ndigits = 8

               write(fmt0,"(i2)") ndigits
               fmt0 = adjustl(fmt0)
               fmt1 = '(a8,I4,2x,f22.'//fmt0//')'
               write (f_out,FMT=fmt1) nampar(i),ivar(i),potparam(i)
            endif
          enddo  ! --- i
          !
        endif 
        !
        still_run = .false.
        !
        ! Print out the ssq for the rovib. energies and pot. data points separetely:
        !
        ssq1 = 0 ; ssq2 = 0 
        !
        wtsum = sum(wt_bit(1:fitting%Nenergies))
        !
        if (wtsum/=0) ssq1 = sqrt( sum(eps(1:fitting%Nenergies)**2*dble(wt_bit(1:fitting%Nenergies)))/wtsum )
        !
        wtsum = sum(wt_bit(1+fitting%Nenergies:npts))
        !
        if (wtsum/=0) ssq2 = sqrt( sum(eps(1+fitting%Nenergies:npts)**2*dble(wt_bit(1+fitting%Nenergies:npts)))/wtsum )

        rms1=sqrt(sum(eps(1:fitting%Nenergies)**2)/fitting%Nenergies)
        rms2=sqrt(sum(eps(1+fitting%Nenergies:npts)**2)/pot_npts)

        !
        write (f_out,6552) fititer,nused,numpar,stadev,ssq1,ssq2,stability
        !
      enddo  ! --- fititer
      !
      call MemoryReport
      !
   enddo outer_loop
   !
   ssq1 = 0 ; ssq2 = 0 
   !
   wtsum = sum(wt_bit(1:fitting%Nenergies))
   !
   if (wtsum/=0) ssq1 = sqrt( sum(eps(1:fitting%Nenergies)**2*dble(wt_bit(1:fitting%Nenergies)))/wtsum )
   !
   wtsum = sum(wt_bit(1+fitting%Nenergies:npts))
   !
   if (wtsum/=0) ssq2 = sqrt( sum(eps(1+fitting%Nenergies:npts)**2*dble(wt_bit(1+fitting%Nenergies:npts)))/wtsum )
   !
   call memoryreport 
   !
  


6552   format(/3X,86('-')/'   |  Iter  | Points | Params |    Deviat     |',&
       '     ssq_ener  |    ssq_pot  | Stability |'/&
       3X,86('-')/,&
       '   | ',I6,' | ',I6,' | ',I6,' |  ',E12.5,' | ',E12.5,'  |  ',&
            E10.3,' |',E10.3,' |',/3X,86('-')/)

6553   format(/3X,86('-')/'   |  Iter  | Points | Params |    Deviat     |',&
       '     rms_ener  |    rms_pot  | Stability |'/&
       3X,80('-')/,&
       '   | ',I6,' | ',I6,' | ',I6,' |  ',E12.5,' | ',E12.5,'  |  ',&
            E10.3,' |',E10.3,' |',/3X,86('-')/)

       !
       call SYSTEM_CLOCK(itime2,irate2,imax2)
       !
       itime=(itime2-itime0)/irate2
       write(f_out,"(/i10,' secs CPU time used'/)") itime
       !

  contains 


    subroutine finite_diff_of_energy(rjacob,Nstates,potparam)
      !
      real(rk) :: rjacob(:,:),potparam(:)
      integer(ik),intent(inout)  :: Nstates(2,0:Nmax,2)
      real(rk),allocatable :: enerleft(:),enerright(:)
      type(calcT) :: calc_(2,0:Nmax,2)
      integer(ik) :: N,Nsize,irot
      !
      ! calculate derivatives with respect to parameters
      !
      allocate (enerright(fitting%Nenergies),enerleft(fitting%Nenergies),stat=alloc)
      call ArrayStart('enerright',alloc,size(enerright),kind(enerright))
      call ArrayStart('enerleft',alloc,size(enerleft),kind(enerleft))
      !
      do iref = 1,2
        do irot = 0,Nmax
          do ipar_ = 1,2
            !
            Nsize = max(1,Nstates(iref,irot,ipar_))
            !
            allocate(calc_(iref,irot,ipar_)%energy(Nsize),calc_(iref,irot,ipar_)%quanta(Nsize),stat=alloc)
            call ArrayStart('calc_%energy',alloc,size(calc_(iref,irot,ipar_)%energy),kind(calc_(iref,irot,ipar_)%energy))
            call ArrayStart('calc_%quanta',alloc,size(calc_(iref,irot,ipar_)%quanta),80_ik)
            !
          enddo
          !
        enddo 
      enddo
      !
      do  ncol=1,numpar
          !
          i = ifitparam(ncol)
          !
          if (verbose>=5) write(f_out,"('Derivative wrt to parameter ',i4,2x,a9)") i,trim(nampar(i))
          !
          tempx=potparam(i)
          deltax=fitfactordeltax*abs(tempx)
          if (deltax .le. 1e-15) deltax=1e-6
          !
          potparam(i)=tempx+deltax
          !
          call run_everest
          !
          call get_vibronic_energies(Nmax,Nstates,calc_)
          ! 
          do iobs = 1,fitting%Nenergies
             !
             Jrot = fitting%obs(iobs)%Jrot
             iparity = fitting%obs(iobs)%iparity
             iref = fitting%obs(iobs)%iref
             nrot = int(jrot)
             !
             ipar_ = 1 ; if (iparity==-1) ipar_ = 2
             !
             if (Jrot>jmax) cycle
             !
             N = fitting%obs(iobs)%N
             !
             if (N>Nstates(iref,nrot,ipar_)) then
               write(out,"('finite_diff_of_energy: illegal N in the observed list:, J,par,N = ',f9.1,i3,i7)") Jrot,iparity,N
               stop 'finite_diff_of_energy: illegal N in the observed list' 
             endif
             !
             enerright(iobs) = calc_(iref,nrot,ipar_)%energy(N)
             ! 
          enddo 
          !
          potparam(i)=tempx-deltax
          !
          call run_everest 
          !
          call get_vibronic_energies(Nmax,Nstates,calc_)
          ! 
          do iobs = 1,fitting%Nenergies
             !
             Jrot = fitting%obs(iobs)%Jrot
             iparity = fitting%obs(iobs)%iparity
             iref = fitting%obs(iobs)%iref
             nrot = int(jrot)
             !
             ipar_ = 1 ; if (iparity==-1) ipar_ = 2
             !
             if (Jrot>jmax) cycle
             !
             N = fitting%obs(iobs)%N
             !
             if (N>Nstates(iref,nrot,ipar_)) then
               write(out,"('finite_diff_of_energy: illegal N in the observed list:, J,par,N = ',f9.1,i3,i7)") Jrot,iparity,N
               stop 'finite_diff_of_energy: illegal N in the observed list' 
             endif
             !
             enerleft(iobs) = calc_(iref,nrot,ipar_)%energy(N)
             ! 
          enddo
          !
          if (verbose>=5) write(f_out,"('      rjacob...'/)")
          !
          rjacob(1:fitting%Nenergies,ncol)=(enerright(1:fitting%Nenergies)-enerleft(1:fitting%Nenergies))/(2.0_rk*deltax)
          !
          rjacob(1:fitting%Nenergies,ncol) = rjacob(1:fitting%Nenergies,ncol) - rjacob(1,ncol)
          !
          potparam(i)=tempx
          !
          if (verbose>=5) write(f_out,"('      ...done'/)")
          !
          if (verbose>=5) write(f_out,"('      Mapping parameters...'/)")
          !
          call map_parameters(dir=.false.)
          !
          if (verbose>=5) write(f_out,"('      ...done'/)")
          !
      enddo ! --- ncol
      !
      ! printing out derivatives for debugging purposes
      !
      if (verbose>=6) then 
         !
         write(f_out,"(/'   Derivative of energies wrt to pot. parameters:')")
         write(f_out,"('    iobs,Jrot,iparity,iref,istate,npar,Energy(obs),DE_i/Df_n')")
         !
         do  ncol=1,numpar
             !
             i = ifitparam(ncol)
             ! 
             do iobs = 1,fitting%Nenergies
                !
                Jrot = fitting%obs(iobs)%Jrot
                iparity = fitting%obs(iobs)%iparity
                iref = fitting%obs(iobs)%iref
                nrot = int(jrot)
                !
                write(out,"(i4,1x,f8.1,3(1x,i2),i5,1x,f15.4,1x,g18.8)") &
                     iobs,Jrot,iparity,iref,fitting%obs(iobs)%quanta%istate,&
                     ncol,fitting%obs(iobs)%Energy,rjacob(iobs,ncol)
                !
             enddo
         enddo
         !
      endif 
      !
      deallocate (enerright,enerleft)
      !
      call ArrayStop('enerright')
      call ArrayStop('enerleft')
      !
      do iref = 1,2
        do irot = 0,Nmax
          do ipar_ = 1,2
            !
            if (associated(calc_(iref,irot,ipar_)%energy)) deallocate(calc_(iref,irot,ipar_)%energy)
            if (associated(calc_(iref,irot,ipar_)%quanta)) deallocate(calc_(iref,irot,ipar_)%quanta)
            !
            calc_(iref,irot,ipar_)%energy=>null()
            calc_(iref,irot,ipar_)%quanta=>null()
            !                                                                
          enddo                                                              
        enddo                                                                
      enddo
      !
      call ArrayStop('calc_%energy') 
      call ArrayStop('calc_%quanta')                                                               
      !
      if (verbose>=5) write(f_out,"('... done!'/)")
      !                                                                      
    end subroutine finite_diff_of_energy                                     
    !
    subroutine map_parameters(dir)
     !
     logical,intent(in) :: dir
     !
     type(fieldT),pointer      :: field
     integer(ik)  :: ifield,iterm,iaddress,iaddress_
     type(linkT),pointer :: flink
     !
     ! mapping 
     !
     if (dir) then 
       !
       do ifield =1,Nfields
         !
         field => poten(ifield)
         !
         do iterm = 1,poten(ifield)%Nterms
           !
           iaddress = poten(ifield)%iaddress(iterm)
           !
           potparam(iaddress) = field%value(iterm)
           nampar(iaddress)   = field%forcename(iterm)
           ivar(iaddress)     = int(field%weight(iterm),ik)
           !
           flink => poten(ifield)%link(iterm)
           !
           if (flink%ifield/=0) then
             !
             iaddress_ = poten(flink%ifield)%iaddress(flink%iparam)
             !
             potparam(iaddress) = potparam(iaddress_)
             !
           endif 
           !
         enddo
         !
       enddo
       !
     else
       !
       do ifield =1,Nfields
         !
         field => poten(ifield)
         !
         do iterm = 1,poten(ifield)%Nterms
           !
           iaddress = poten(ifield)%iaddress(iterm)
           !
           field%value(iterm) = potparam(iaddress)
           field%forcename(iterm) = nampar(iaddress)
           field%weight(iterm) = real(ivar(iaddress),rk)
           !
           flink => poten(ifield)%link(iterm)
           !
           if (flink%ifield/=0) then
             !
             poten(ifield)%value(iterm) = poten(flink%ifield)%value(flink%iparam)
             !
           endif 
           !
         enddo
         !
       enddo
       !
     endif 
     !
   end subroutine map_parameters
   !

    subroutine print_parameters
     !
     type(fieldT),pointer      :: field
     integer(ik)  :: ifield,iterm
     type(linkT),pointer :: flink
     character(len=3) :: fit
     character(len=130) :: fmt_str
     !
     fmt_str = "(a8,1x,3(i3),1x,e22.14,2x,a3,2x,a4,1x,2i4)"
     !
     do ifield =1,Nfields
       !
       field => poten(ifield)
       !
       write (f_out,"(/'Poten',2x,i4)") ifield
       write (f_out,"(a)") trim(field%name)
       write (f_out,"(a)") "values"
       !
       do iterm = 1,poten(ifield)%Nterms
         !
         flink => poten(ifield)%link(iterm)
         !
         if (flink%ifield/=0) then
           !
           write (f_out,fmt_str)  poten(ifield)%forcename(iterm),&
                                            poten(ifield)%ipower1(iterm),&
                                            poten(ifield)%ipower2(iterm),&
                                            poten(ifield)%ipower3(iterm),&
                                            poten(ifield)%value(iterm),"   ",&
                                            "Link",flink%ifield,flink%iparam
                                            !
         elseif (poten(ifield)%weight(iterm)>0) then
           !
           write (f_out,fmt_str)  poten(ifield)%forcename(iterm),&
                                            poten(ifield)%ipower1(iterm),&
                                            poten(ifield)%ipower2(iterm),&
                                            poten(ifield)%ipower3(iterm),&
                                            poten(ifield)%value(iterm),"fit"
         else
           !
           write (f_out,fmt_str)  poten(ifield)%forcename(iterm),&
                                            poten(ifield)%ipower1(iterm),&
                                            poten(ifield)%ipower2(iterm),&
                                            poten(ifield)%ipower3(iterm),&
                                            poten(ifield)%value(iterm)
           !
         endif 
         !
         !
       enddo
       !
       write (f_out,"(a)") "end"
       !
     enddo
     !
   end subroutine print_parameters
   !



    ! check which field and the field address of the given potparameter 
    !
    subroutine which_field_address(iparamater,ifield,iaddr)
     !
     integer(ik),intent(in)  :: iparamater
     integer(ik),intent(out) :: ifield,iaddr
     !
     type(fieldT),pointer      :: field
     !
     if (iparamater<=poten(1)%Nterms) then 
         ifield = 1
         iaddr = iparamater
     elseif (iparamater<=poten(1)%Nterms+poten(2)%Nterms) then
         ifield = 2
         iaddr = iparamater-poten(1)%Nterms
     elseif (iparamater<=poten(1)%Nterms+poten(2)%Nterms+poten(3)%Nterms) then
         ifield = 3
         iaddr = iparamater-(poten(1)%Nterms+poten(2)%Nterms)
     elseif (iparamater<=poten(1)%Nterms+poten(2)%Nterms+poten(3)%Nterms+poten(4)%Nterms) then
         ifield = 4
         iaddr = iparamater-(poten(1)%Nterms+poten(2)%Nterms+poten(3)%Nterms)
     else
       stop 'which_field_address: iparameter is out of range'
     endif
     !
   end subroutine which_field_address
   !   
   !
   subroutine run_everest
      !
      logical  ::  SYSTEMQQ  ! system function for  calling extyernal programs, 
                             ! platform dependent.  
      !
      if (verbose>=5) write(f_out,"('   Everesting...')")
      !
      if (verbose>=5) write(f_out,"('   Mapping parameters...')")
      call map_parameters(dir=.false.)
      !
      if (verbose>=5) write(f_out,"('   Prepating input files with pot. parameters...')")
      !
      call write_input_potential_parameters(Nfields,poten)
      !
      if (verbose>=3) write(f_out,"(/'   Calling Everest program with kmax =  ',i5,'...')") nint(kmax+0.5)
      !
#if (debug_ == 0)
         isys = systemqq(evib_exe//'< evib.inp > evib.out')
#endif
      !
      if (verbose>=3) write(f_out,"('   Calling evrot Everest program for J=0.5...',f8.1)") Jmax
      !
#if (debug_ == 0)
         isys = systemqq(erot_exe//'<erot.inp >erot.out')
         !
#endif
         
#if (debug_ == 0)
         isys = systemqq('./collect_energies_rot.sh erot.out > erot.log')
#endif
      !
      if (verbose>=3) write(f_out,"('   ... done!')")
      !
   end subroutine run_everest
   !

    subroutine write_input_potential_parameters(Nfields,poten)
    !
    integer(ik),intent(in) :: Nfields
    type(fieldT),target,intent(in)  :: poten(Nfields) ! field1.fit
    character(len=10) :: filename
    character(len=1)  :: fit_str = " "
    type(fieldT),pointer      :: field
    integer(ik) :: f_p = 21,ifield,i
    type(linkT),pointer :: flink
      !
      do ifield = 1,Nfields
        !
        field => poten(ifield)
        !
        if (ifield<10) then 
           write(filename,"('field',i1,'.fit')") ifield
        else
           write(filename,"('field',i2,'.fit')") ifield
        endif
        !
        open(f_p,file=trim(filename),status='replace')
        !
        write(f_p,*) field%Nterms
        write(f_p,"('r1',4x,f19.10)") field%r1
        write(f_p,"('r2',4x,f19.10)") field%r2
        write(f_p,"('alpha',4x,f19.10)") field%alpha
        !
        do i=1,field%Nterms
           !
           fit_str = " " ; if (field%weight(i)>0) fit_str = "*"
           !
           ! check if the paramter is linked to a varied paramter; it needs to be marked as varied in EVEREST  with *
           !
           flink => poten(ifield)%link(i)
           !
           if (flink%ifield/=0) then
             !
             if (poten(flink%ifield)%weight(flink%iparam)>0) fit_str = "*"
             !
           endif
           !
           ! write the record into the input file
           !
           write(f_p,"(a7,a1,1x,3i4,1x,e24.12)") adjustl((field%forcename(i))),fit_str,&
                                                field%ipower1(i),field%ipower2(i),field%ipower3(i),&
                                                field%value(i)
        end do
        !
        write(f_p,*) field%alphamin
        write(f_p,*) field%range(1,1),field%range(2,1),field%range(3,1)
        write(f_p,*) field%range(1,2),field%range(2,2),field%range(3,2)
        write(f_p,*) field%range(1,3),field%range(2,3),field%range(3,3)
        !
        close(f_p)
        !
      enddo
      !
    end subroutine write_input_potential_parameters


   !
  end subroutine fitting_energies
  !
  !
  !
  subroutine linur(dimen,npar,coeff,constant,solution,error)

  integer,intent(in)  :: dimen,npar
  integer,intent(out) :: error 
  real(rk),intent(in)  :: coeff(npar,npar),constant(npar)
  real(rk),intent(out) :: solution(npar)
  real(rk)          :: a0(npar,npar)
  real(rk)          :: c
  integer                   :: i1,i2,i,k8,k,l,k9

  !----- begin ----!
  
    do i1=1,dimen
    do i2=1,dimen 
       a0(i1,i2)=coeff(i1,i2)
    enddo
    enddo

    do i=1,dimen
      solution(i)=constant(i)
    enddo
    error=0
    do i=1,dimen
      c=0
      k8=i-1
      do k=1,k8
        c=c+a0(k,i)*a0(k,i)
      enddo

      if (c.ge.a0(i,i)) then
      !      write(f_out,*) '(',i,'-th adj. parameter is wrong )'
       error=i
       return
      endif

      a0(i,i)=sqrt(a0(i,i)-c)
      if (a0(i,i).eq.0) then
      !      write(f_out,*) '(',i,'-th adj. parameter is wrong )'
         error=i
         return
      endif
      k8=i+1
      do l=k8,dimen
         k9=i-1
         c=0.0
         do k=1,k9 
            c=c+a0(k,i)*a0(k,l)
         enddo
         a0(i,l)=(a0(l,i)-c)/a0(i,i)
      enddo
    enddo
    do i=1,dimen
      k8=i-1
      c=0.0
      do k=1,k8
         c=c+a0(k,i)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
    do i1=1,dimen
      i=1+dimen-i1
      k8=i+1
      c=0.0
      do k=k8,dimen
          c=c+a0(i,k)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
  return
  end subroutine linur


!------------------------------------------!
  subroutine invmat(al,ai,dimen,npar)
  integer,intent(in)           :: npar,dimen
  real(rk),intent(in)  :: al(npar,npar)
  real(rk),intent(out) :: ai(npar,npar)
  real(rk)             :: h(npar),p,q
  integer                      :: i1,k,i,j,k8,k9
      

    ai(1:dimen,1:dimen)=al(1:dimen,1:dimen)
 
    do i1=1,dimen
      k=dimen-i1+1
      p=ai(1,1)
      do i=2,dimen
        q=ai(i,1)
        h(i)=q/p
        if(i.le.k) h(i)=-q/p
        do j=2,i
          k8=i-1
          k9=j-1
          ai(k8,k9)=ai(i,j)+q*h(j)
        enddo 
      enddo 
      ai(dimen,dimen)=1.0/p
      do i=2,dimen
        k8=i-1
        ai(dimen,k8)=h(i)
      enddo 
   end do 
   do i=1,dimen
     k8=i-1
     do j=1,k8
       ai(j,i)=ai(i,j)
     enddo 
   enddo 
   return
 end subroutine invmat
!------------------------------------------!


    subroutine get_vibronic_energies(Nmax,Nstates,calc)
      !
      integer(ik),intent(in) :: Nmax
      integer(ik),intent(inout)  :: Nstates(2,0:Nmax,2)
      type(calcT),intent(inout) :: calc(2,0:Nmax,2)
      !
      integer(ik)              :: tunit,tau,Kc,Pot,Kma,ISR,J2,ipar,Nrot,Ntot,Nsize,nn
      integer(ik)    :: i,iref,v1,v2,v3,ka
      real(rk)  :: Energy,De,wei,Ome,Wome,Wma,Jrot
      character(len=5)   :: Jch
      integer(ik)  :: Nstates_(2,0:Nmax,2)
      !
      tunit = 22
      !
      if (verbose>=5) write(f_out,"('   Extracting rovibronic energies ...')")
      !
      open(tunit,file='erot.log',action='read',status='old')
      !
      !      Jch,Par,j1,        Energy            De     V  Ka Kc    v1 v2 v3     Wei    Ome    Wome Pot Kma    Wma   ISR  N
      !
      !
      Nstates_ = 0
      !
      Ntot = sum(Nstates)
      !
      ! initial count of calculated state 
      !
      if (Ntot==0) then 
         !
         do 
           read(tunit,"(a5,2x,i2,1x,i5,2(1x,f12.4),1x,i3)",end=119) Jch,tau,nn,Energy,De,iref !Ka,Kc,v1,v2,v3,Wei,Ome,Wome,Pot,Kma,Wma,ISR,  N
           !
           do i=1,len_trim(Jch)
               if (Jch(i:i) == "/") Jch(i:i) = " "
           end do
           !
           read(Jch,*) J2
           Jrot = J2*0.5_rk
           Nrot = nint(Jrot-0.5_rk)
           ipar  = 1 ; if (tau == -1) ipar = 2
           !
           if (Nrot>Nmax) cycle
           !
           ! assing ref to 1 (V=1) or 2 (V=2,3)
           iref = min(iref,2)
           !
           Nstates_(iref,Nrot,ipar) = Nstates_(iref,Nrot,ipar) + 1
           !
           cycle
           !         
         119 continue
             exit
     
         end do
         !
         Ntot = sum(Nstates_)
         !
         if (Ntot==0) then  
            write(out,"('get_vibronic_energies: no states to read ')")
            stop 'get_vibronic_energies: no states to read '
         endif
         !
         Nstates = Nstates_
         !
         rewind(tunit)
         !
         return
         !
      endif 
      !
      Nstates_ = 0
      !
      do 
        read(tunit,"(a5,2x,i2,1x,i5,2(1x,f12.4),1x,i3,2(1x,i3),3(1x,i2),1x,f5.3,f4.1,1x,f5.3,1x,i2)",end=121) &
             Jch,tau,nn,Energy,De,iref,Ka,Kc,v1,v2,v3,Wei,Ome,Wome,Pot !,Kma,Wma,ISR,  N
        !
        do i=1,len_trim(Jch)
            if (Jch(i:i) == "/") Jch(i:i) = " "
        end do
        !
        read(Jch,*) J2
        Jrot = J2*0.5_rk
        Nrot = nint(Jrot-0.5_rk)
        ipar  = 1 ; if (tau == -1) ipar = 2
        !
        if (Nrot>Nmax) cycle
        !
        iref = min(iref,2)
        !
        Nstates_(iref,Nrot,ipar) = Nstates_(iref,Nrot,ipar) + 1
        !
        nn = Nstates_(iref,Nrot,ipar)
        !
        if (nn>size(calc(iref,Nrot,ipar)%energy)) cycle
        !
        calc(iref,Nrot,ipar)%energy(nn) = Energy
        calc(iref,Nrot,ipar)%quanta(nn)%v1 = v1
        calc(iref,Nrot,ipar)%quanta(nn)%v2 = v2
        calc(iref,Nrot,ipar)%quanta(nn)%v3 = v3
        calc(iref,Nrot,ipar)%quanta(nn)%istate = Pot
        calc(iref,Nrot,ipar)%quanta(nn)%Omega = Ome
        calc(iref,Nrot,ipar)%quanta(nn)%ilambda = Ka
        !
        cycle
        !         
      121 continue
          exit

      end do
      !
      if (verbose>=5) write(f_out,"('   ... done!')")
      !
      close(tunit,status='keep')
      !
    end subroutine get_vibronic_energies
    !

    subroutine get_vibronic_derivatives(Nmax,Nparams,calc,iparamfit,rjacob)
      !
      integer(ik),intent(in) :: Nmax,Nparams
      integer(ik),intent(in)  :: iparamfit(Nparams)
      real(rk) :: rjacob(:,:)
      type(calcT),intent(in) :: calc(2,0:Nmax,2)
      !
      integer(ik)              :: tunit,tau,Kc,Pot,Kma,ISR,J2,ipar,Nrot,Ntot,Nsize,nn,n1,n2,n3
      integer(ik)    :: i,iref,v1,v2,v3,ka,iab2,iobs,iparam,iState,iaddress,ifit,ifit_,iparam_
      real(rk)  :: Energy,De,wei,Ome,Wome,Wma,Jrot,Deriv
      character(len=5)   :: Jch
      type(linkT),pointer :: flink
      !
      tunit = 32
      !
      if (verbose>=5) write(f_out,"('   Extracting rovibronic energies ...')")
      !
      open(tunit,file='evrder.dat',action='read',status='old')
      !
      ! skip first line
      !
      read(tunit,*) Jch
      !
      ! set all derivatives to zero
      !
      rjacob = 0 
      !
      do 
        read(tunit,*,end=129) iref,J2,tau,iab2,nn,iState,iparam_,iparam,n2,n3,Energy,Deriv
        !
        Jrot = real(J2,rk)*0.5_rk
        Nrot = nint(Jrot-0.5_rk)
        ipar  = 1 ; if (tau == -1) ipar = 2
        !
        if (Nrot>Nmax) cycle
        !
        ! check just in case if the energy agtrees for the state from calc
        !
        if (abs(Energy-calc(iref,Nrot,ipar)%energy(nn))>1e-2) then 
          write(out,"('get_vibronic_deriv er: internal/external energies dont agree for 2J,tau,state,nn,E = ',4i5,2(1x,f12.3))") &
                    J2,tau,istate,nn,calc(iref,Nrot,ipar)%energy(nn),Energy
          stop 'get_vibronic_deriv er: internal/external energies dont agree'
        endif
        !
        iobs = 0
        !
        ! Find a match in the observed set of states
        !
        loop_iobs : do 
          !
          iobs = iobs+1
          if (nn==fitting%obs(iobs)%N.and.nint(Jrot-fitting%obs(iobs)%Jrot)==0 &
              .and.iref==fitting%obs(iobs)%iref &
              .and.tau==fitting%obs(iobs)%iparity ) then
            !
            exit loop_iobs
            !
          endif
          !
          if (iobs==fitting%Nenergies) exit loop_iobs 
          !
        enddo loop_iobs
        !
        if (iobs>=fitting%Nenergies) cycle 
        !
        ! check if the parameter is linked to a different object 
        !
        flink => poten(istate)%link(iparam)
        !
        ! obtain the address of the potentila parameters in the global register 
        !
        if (flink%ifield==0) then 
           !
           iaddress = poten(istate)%iaddress(iparam)
           ifit = iparamfit(iaddress)
        else
           iaddress = poten(flink%ifield)%iaddress(flink%iparam)
           ifit = iparamfit(iaddress)
        endif
        !
        if (ifit==0) then 
          write(out,"(/'Error in get_vibronic_deriv: illegal parameter in EV derivs, 2J,tau,state,nn,E,ipar = ',4i5,1x,f12.3)") &
                    J2,tau,istate,nn,calc(iref,Nrot,ipar)%energy(nn),iparam
          stop 'Error in get_vibronic_deriv: illegal parameter in EVEREST derivatives'
          !
        endif
        !
        ! Storing the derivatives by adding to the previosuly defined values in case of the linked parameters:
        ! which are defined as d E / d C =  d E / d C1 + d E / d C2 
        ! for the constraint C1=C2=C
        !
        rjacob(iobs,ifit) = rjacob(iobs,ifit) + deriv
        !
        continue
        !
        cycle
        !         
      129 continue
          exit
     
      end do
      !
      ! subtract the ZPE derivatives
      !
      do iobs = 2,fitting%Nenergies
        rjacob(iobs,:) = rjacob(iobs,:) - rjacob(1,:)
      enddo
      rjacob(1,:) = 0 
      !
      if (verbose>=5) write(f_out,"('   ... done!')")
      !
      close(tunit,status='keep')
      !
    end subroutine get_vibronic_derivatives

    !
    subroutine update_linked_parameters
      !
      !integer(ik),intent(in) :: Nparams
      !real(rk),intent(inout) :: potparam(Nparams)
      !
      integer(ik) :: ifield,iterm,iaddress
      type(linkT),pointer :: flink
      !
      !
      ! update linked parameters 
      !
      do ifield =1,Nfields
        !
        do iterm = 1,poten(ifield)%Nterms
          !
          flink => poten(ifield)%link(iterm)
          !
          if (flink%ifield/=0) then
            !
            poten(ifield)%value(iterm) = poten(flink%ifield)%value(flink%iparam)
            !
            !iaddress = poten(ifield)%iaddress(iterm)
            !
            !potparam(iaddress) = poten(ifield)%value(iterm)
            !
          endif 
          !
        enddo
        !
      enddo
      !
    end subroutine update_linked_parameters
    !
    !
    ! Defining potential energy function 
    !
    subroutine poten_func(field,r12,r32,theta,f)
      !
      implicit none
      !
      type(fieldT),intent(in)  :: field ! field1.fit
      integer,parameter           :: N = 99
      !
      double precision,intent(in) :: r12,r32,theta
      double precision,intent(out):: f
      double precision            :: f0,alphae
      !
      alphae = field%alpha*pi/180.0_rk
      !
      f0 = poten_xyz(field%Nterms,field%r1,field%r2,alphae,field%ipower1,field%ipower2,field%ipower3,&
                     field%value,field%r1,field%r2,alphae)
      !
      f = poten_xyz(field%Nterms,field%r1,field%r2,alphae,field%ipower1,field%ipower2,field%ipower3,field%value,r12,r32,theta)
      !
      f = f-f0
      !
     end subroutine poten_func


    !
    subroutine create_the_deriv_job_file(MXRoot,iOMP)
	  !
      integer,intent(in)  :: MXRoot,iOMP
      integer :: Mem = 1000,iprint = 1, tunit
      !
      tunit = 42
      !
      open(tunit,file='vde.inp',action='write',status='replace')
      !
      write(tunit,"(' &EVEREST')")
      write(tunit,"('Memo')")
      write(tunit,"(i8)") Mem 
      write(tunit,"('print')") 
      write(tunit,"(i8)") iprint
      write(tunit,"('OMP')") 
      write(tunit,"(i4)") iOMP
      write(tunit,"('End')") 
      !
      close(tunit,status='keep')
      !
      tunit = 43
      !
      open(tunit,file='rde.inp',action='write',status='replace')
      !
      write(tunit,"(' &EVEREST')")
      write(tunit,"('Memo')")
      write(tunit,"(i8)") Mem 
      write(tunit,"('MXRoot')")
      write(tunit,"(i8)") MXRoot 
      write(tunit,"('print')") 
      write(tunit,"(i8)") iprint
      write(tunit,"('End')") 
      !
      close(tunit,status='keep')
      !
    end subroutine create_the_deriv_job_file
    !
    !
    subroutine read_jacobi_matrix(enermax,parmax,derj)
	  !
      integer,intent(in)  :: enermax,parmax
      real(8),intent(out) :: derj(enermax,parmax)
      integer             :: ie,ipar
      !
      open(1, file='fort.12',status='old')
      !
      !derj(1,:) = 0 
      !
      do ie = 1,enermax 
        !
        read(1,*)
        !
        do ipar = 1,parmax 
          !
          read(1,*) derj(ie,ipar)
          !    
        enddo
        !
      enddo
      !
      close(1,status='keep')
      !
    end subroutine read_jacobi_matrix

    subroutine robust_fitting(sigma,eps,wt)

      !real(8),intent(inout) :: a_wats
      real(rk),intent(in)    :: sigma(:),eps(:)
      real(rk),intent(inout) :: wt(:)
      !
      integer(ik)            :: npts,i,nrow,nused
      real(rk)              :: da1,da2,wtsum,da,a_wats
      !
      npts = size(sigma)
      !
      nused = 0 
      do i=1,npts
        if (wt(i)>0) nused=nused+1
      enddo
      !
      !Watson alpha-parameter
      ! 
      a_wats = 0.001_rk
      !
      if (verbose>=4) write(out,"('Watson parameter =',f18.8)") a_wats
      !
      !Watson alpha-parameter
      ! 
      !do i = 1,-1
      !  !
      !  da1 = 0
      !  da2 = 0
      !  !
      !  do nrow=1,npts
      !    if (wt(nrow)>small_) then 
      !      da1 = da1+eps(nrow)**2/( sigma(nrow)**2+a_wats*eps(nrow)**2 )
      !      da2 = da2+eps(nrow)**4/( sigma(nrow)**2+a_wats*eps(nrow)**2 )**2
      !    endif 
      !  enddo 
      !  !
      !  da =( da1 -  real(nused-numpar,rk) )/da2
      !  a_wats = a_wats + da
      !  !
      !  if (a_wats<sqrt(small_)) a_wats = 1e-3+real(i,rk)*1e-2
      !  !
      !enddo
      ! 
      !  adjusting the weights  
      ! 
      do nrow=1,npts
         if (wt(nrow)>0) wt(nrow) = 1.0d0/( sigma(nrow)**2 + a_wats*eps(nrow)**2 )
      enddo 
      ! 
      ! "re-normalizing" the weight factors
      !
      wtsum = sum(wt(1:npts))
      wt(1:npts) = wt(1:npts)/wtsum
      !
    end subroutine robust_fitting
 

  end module fit_module
  !
  program driver
    use fit_module

    call fitting_energies()

  end program driver


