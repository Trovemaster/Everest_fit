module fit_module
!
#define debug_  1
!
use input
use timer
use accuracy
!
implicit none
!
integer, parameter    :: verbose  = 3 ! Verbosity level    
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
real(rk),parameter :: stadev_best=1e-04,stab_best=1e-12  ! best standard error and srability 
                                             ! the fit will finish when these reached. 
!
real(rk),parameter :: fitfactordeltax=0.001   ! parameter for the finite differncies differention 
!
character(len=70) :: deriv_type = 'finite  '  ! hellman or direct - how we calculate derivatives. 
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
! type to deal with the calculated Everest energies 
!
type  FTEverestT
  real(8),pointer   ::  energy(:)   => null()   ! to store the calculated energies for  given (J,symmetry)
  real(8),pointer   ::  derj(:,:)   => null() ! to store the calculated derivatives 
end type FTEverestT
!
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
    !
  end type fieldT
  !
  type quantaT
    real(rk) :: I1 ! nuclear spin 1
    real(rk) :: F1  ! \hat{F1} = \hat{I1} + \hat{J}
    INTEGER(ik) :: index_F1 ! index of the F1 in F1_list
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
  character(len=30) :: evib_exe               ! Everest exe-file. 
  character(len=30) :: erot_exe              ! Everest exe-file. 
  character(len=30) :: xpect_exe               ! Everest exe-file. 

  integer           :: alloc, ierror      ! error state variables 
  !
  integer,parameter  :: jlist_max = 500
  !
  real(rk)              :: j_list_(1:jlist_max)=-1.0_rk,jmin,jmax,kmax,jrot2,f_t
  !
  character(len=3)   :: kmax_ch,jmax_ch
  !
  real(rk),allocatable  :: j_list(:)
  !
  character(len=cl) :: w
  !
  logical :: eof
  !
  integer(ik) :: i,nJ,ipot,iso,iref,istate,jref,Nparam,Nparam_check,iparam
  !
  type(fieldT),pointer      :: field
  !
  integer(ik)       :: iut !  iut is a unit number. 
  integer(ik)  :: iobs,itau,iai,ic,iparams,irange,iparity
  logical :: matchfound
  character(len=wl) :: large_fmt
  character(len=cl) :: ioname
  !
  integer(ik)  :: nmax,nrot
  integer(ik),allocatable  :: Nstates(:,:,:)
  type(calcT),allocatable  :: calc(:,:,:)


  integer           :: npts,pot_npts,nused      ! number of data points: current, all obs. energies, all pot. points, and actuall used. 
  !
  integer           :: nrow,ncol               ! loop-variables
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
  !
  !integer,allocatable :: J_obs(:),sym_obs(:) ! Arrays where we 
  !                                           ! store the information about obs. energies: 
  !                                           ! J, symmetry, number in the block, vib. quantum numbers
  !
  ! objects to store the energies 
  real(rk),allocatable :: enercalc(:),ener_obs(:),eps(:)
  !
  type(FTEverestT)                    :: Everest(2)   ! caclualted term Everest values;
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
  real(rk) :: ssq1,ssq2,rms1,rms2,info_
  logical          :: still_run
  integer          :: j,l
  integer          :: iener,irow,icolumn,ipar_
  !
  !     For the details on the robust fit see J.K.G. Watson, JMS, 219, 326 (2003).
  !     The idea is optimize the weigth factors iteratevely depending on the fitting results, as  follows:
  !     w_i = 1 / sigma_i^2+alpha*eps_i^2 
  !     where sigma_i represent the experimental data precision, alpha is the Watson parameter (can be also adjsuted) and 
  !     eps_i = Obs - Calc deviatons. 
  !
  real(rk) :: robust_fit = 0 
  !                   robust_fit = 0 corresponds to the standard fit, the robust fittting procedure is switched off. 
  !                   robust_fit <> =0, the robust fittting procedure is on.
  !                   robust_fit also defines the target accuracy of the fit. For example for robust_fit  = 0.001 
  !                   all data with (obs-calc) < 0.001 will be considered as relatively bad (i.e. outliers ) 
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
   allocate(poten(nestates),stat=alloc)
   call ArrayStart('poten',alloc,1,4)
   allocate(spinorbit(ncouples),stat=alloc)
   call ArrayStart('spinorbit',alloc,1,4)
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
      case('ROBUST')
        !
        call readf(robust_fit)
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
        !
      case('DERIV','XPECT')
        !
        call reada(xpect_exe)
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
              call ArrayStart('field%forcename',alloc,size(field%value),kind(field%value))
              call ArrayStart('field%forcename',alloc,size(field%forcename),kind(field%forcename))
              call ArrayStart('field%forcename',alloc,size(field%weight),kind(field%weight))
              !
              allocate(field%ipower1(Nparam),field%ipower2(Nparam),field%ipower3(Nparam),stat=alloc)
              call ArrayStart('field%forcename',alloc,size(field%ipower1),kind(field%ipower1))
              call ArrayStart('field%forcename',alloc,size(field%ipower2),kind(field%ipower2))
              call ArrayStart('field%forcename',alloc,size(field%ipower3),kind(field%ipower3))
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
           case('FIT_TYPE')
             !
             call readu(fitting%fit_type)
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
                ! iref = 1 (state=1) and iref=2 (state=2,3)
                !
                fitting%obs(iobs)%iref = min(fitting%obs(iobs)%quanta%istate,2)
                !
                call readi(fitting%obs(iobs)%quanta%v1)
                call readi(fitting%obs(iobs)%quanta%v2)
                call readi(fitting%obs(iobs)%quanta%v3)
                call readf(fitting%obs(iobs)%quanta%F)
                call readi(fitting%obs(iobs)%quanta%ilambda)
                !
                call readf(fitting%obs(iobs)%weight)
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
   ! sigma = exp. data precision for robust fit 
   !
   if (fitting%robust>0) then
     !
     sigma = 1.0d0
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
   !
   allocate (ener_obs(fitting%Nenergies),enercalc(fitting%Nenergies),stat=alloc)
   call ArrayStart('ener_obs',alloc,size(ener_obs),kind(ener_obs))
   call ArrayStart('enercalc',alloc,size(enercalc),kind(enercalc))
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
   call map_parameters(dir=.true.)
   !   
   parold=potparam
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
      numpar  = 0
      !
      ! Count the actually varying number of parameters:
      !
      numpar  = 0
      do i=1,total_parameters
        if (ivar(i) .gt. 0) numpar=numpar+1
      enddo 
      !
      rjacob = 0 
      !
      ! Prepaper the first run:
      !
      !!!open(1,file='evib.out',status='replace') ; write(1,"('Start fitting...')") ; close(1)
      !!!open(1,file='erot.out',status='replace') ; write(1,"('Start fitting...')") ; close(1)
      !!!open(1,file='xpect.out',status='replace')  ; write(1,"('Start fitting...')") ; close(1)
      !
      !isys = systemqq('rm evib.out erot.out xpect.out')
      !
      Nstates = 0
      !
      ! The loop starts here. 
      !
      do while( fititer.le.fitting%itermax .and. stadev.ge.stadev_best.and. stability.ge.stab_best)
        !
        fititer = fititer + 1
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
        call write_potential_input_parameters(poten(1),'field1.fit')
        call write_potential_input_parameters(poten(2),'field2.fit')
        call write_potential_input_parameters(poten(3),'field3.fit')
        call write_potential_input_parameters(spinorbit(1),'field4.fit')
        !
        if (verbose>=3) write(f_out,"(/'calling Everest program')")
        !
        write(kmax_ch,"(i3)") nint(kmax+0.5)
        write(jmax_ch,"(i3)") nint(2.0_rk*jmax)
        !
        ! generate the vib-input and rot-input files using the tempates provided
        !
#if (debug_ == 0)
           isys = systemqq('vib_inp.sh '//kmax_ch)
           !
           isys = systemqq('rot_inp.sh '//kmax_ch//' '//jmax_ch)
#endif
        !
        if (verbose>=3) write(f_out,"('calling the vibrational evib.x program, kmax = ',f8.1)") kmax
        !
#if (debug_ == 0)
          isys = systemqq(evib_exe//'< evib.inp > evib.out')
#endif
        !
        if (verbose>=3) write(f_out,"('calling the rovibronic erot.x program, kmax, jmax= ',2f8.1)") kmax,jmax
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
        call get_vibronic_energies(Nmax,Nstates,calc)
        !
        ! Zero point energy:
        !
        ezero_ = calc(1,0,1)%energy(1)
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
          if (trim(deriv_type)/='hellman'.and.fitting%itermax.ge.1.and.fitting%factor>1e-12) then
            !
            call finite_diff_of_energy(rjacob,Nstates,potparam)
            !
          else
            !
            write(out,"('hellmanm is not imeplemented yet')")
            !
            stop 'hellmanm is not imeplemented yet'
            !
            !
            ! Jacobi matrix stores derivatives only for energies that have obs. counterparts 
            ! in input with weight /= 0, i.e. participating in the fit. 
            !
            do nrow = 1,fitting%Nenergies
              !
              !iJ = J_obs(nrow) ; jsym = sym_obs(nrow) !; ipar = mod(isym,2)
              !
              !if (iJ==Jrot.and.(isym==jsym.or.iasym==jsym)) then 
              !  !
              !  iener  = n_obs(nrow) ;  if (iener>meval.or.isym>4.or.iasym>4) cycle 
              !  !
              !  ! Use only the derivatives for parameters with ivar/=0, i.e. varying paramaters. 
              !  !
              !  ncol=0
              !  do  i=1,total_parameters
              !    if (ivar(i) .ne. 0) then
              !       !
              !       ncol=ncol+1
              !       !
              !       rjacob(nrow,ncol) = Everest(jsym)%derj(iener,i)-derj0(i)
              !       !
              !     endif
              !     ! 
              !  enddo 
              !  !
              !endif 
              !
            enddo
            !
          endif 
          !
        endif 
        !
        ! Print out the calc. and obs.-calc., i.e. result of the fit. 
        ! Only fitted energies are printed. 
        !
        write(f_out,"(/1x,100('-'))")
        write(f_out,"('| ## |  N |      J  |par |    Obs.    |   Calc.    | Obs.-Calc. |   Weight |    v1 v2 v3 St    Ome  |')")
        write(f_out,"(1x,100('-'))")
        !
        eps = 0
        !
        do nrow = 1,fitting%Nenergies
          !
          Jrot = fitting%obs(nrow)%Jrot
          iparity = fitting%obs(nrow)%iparity
          iref = fitting%obs(nrow)%iref
          nrot = int(jrot)
          !
          ipar_ = 1 ; if (iparity==-1) ipar_ = 2
          !
          if (Jrot>jmax) cycle
          !
          iener = fitting%obs(nrow)%N
          !
          if (iener>Nstates(iref,Nrot,ipar_)) then 
             enercalc(nrow) = 0
          else
             enercalc(nrow) = calc(iref,Nrot,ipar_)%energy(iener)-ezero_
          endif
          !
          eps(nrow) = fitting%obs(nrow)%energy-enercalc(nrow)
          !
          write (f_out,"(2i5,f8.1,i5,2x,3f13.4,2x,e9.2,4x,4(i3),f8.1,2x,3(i3))") &
             nrow,iener,Jrot,iparity,enercalc(nrow)+eps(nrow),enercalc(nrow),-eps(nrow),&
             wtall(nrow),&
             calc(iref,Nrot,ipar_)%quanta(iener)%v1,&
             calc(iref,Nrot,ipar_)%quanta(iener)%v2,&
             calc(iref,Nrot,ipar_)%quanta(iener)%v3,&
             calc(iref,Nrot,ipar_)%quanta(iener)%istate,&
             calc(iref,Nrot,ipar_)%quanta(iener)%omega,&
             fitting%obs(nrow)%quanta%v1,&
             fitting%obs(nrow)%quanta%v2,&
             fitting%obs(nrow)%quanta%v3
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
        write(f_en,"('| ## |  N |      J  |par |    Obs.    |   Calc.    | Obs.-Calc. |   Weight |    v1 v2 v3 St    Ome  |')")
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
                  nrow = 0
                  !
                  loop_nrow : do 
                    !
                    nrow = nrow+1
                    if (irow==fitting%obs(nrow)%N.and.nint(Jrot-fitting%obs(nrow)%Jrot)==0 &
                        .and.iref==fitting%obs(nrow)%iref &
                        .and.iparity==fitting%obs(nrow)%iparity ) then
                      !
                      exit loop_nrow
                      !
                    endif
                    !
                    if (nrow==fitting%Nenergies) exit loop_nrow 
                    !
                  enddo loop_nrow
                  !
                  if (nrow<fitting%Nenergies) then
                     !
                     mark = " "
                     !
                     if (abs(eps(nrow))>1.0) mark = "!"
                     !
                     write(f_en,"(2i5,f8.1,i5,' ',3f13.4,2x,e9.2,2x,a1,2x,4(i3),f8.1,2x,3(i3))") &
                                     irow,fitting%obs(nrow)%n,Jrot,fitting%obs(nrow)%iparity,&
                                     !fitting%obs(nrow)%iref,&
                                     enercalc(nrow)+eps(nrow),enercalc(nrow),-eps(nrow),&
                                     wtall(nrow),mark,&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%v1,&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%v2,&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%v3,&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%istate,&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%omega,&
                                     fitting%obs(nrow)%quanta%v1,&
                                     fitting%obs(nrow)%quanta%v2,&
                                     fitting%obs(nrow)%quanta%v3

                     !
                  else
                     !
                     write(f_en,"(2i5,f8.1,i5,' ',3f13.4,2x,e9.2,5x,4(i3),f8.1)") irow,0,Jrot,ipar_,0.0_rk,&
                                     calc(iref,Nrot,ipar_)%energy(irow)-ezero_,0.0_rk,0.0_rk,&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%v1,&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%v2,&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%v3,&
                                     calc(iref,Nrot,ipar_)%quanta(irow)%istate,&
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
        ! Here the potential energy section starts. 
        !
        do nrow=1,pot_npts
          !
          !  we assume here that the angles are written in degrees, 
          !  not in radians. I.e. we have to convert the degrees to in radians.
          !     
          r1 =fitting%r1(nrow)
          r2 =fitting%r2(nrow)
          theta=fitting%r3(nrow)*pi/180.0_rk
          !
          istate = fitting%iPES(nrow)
          !
          ! Call the potential function. It has to be included into the "poten" subroutine
          !
          call poten_func(potparam,v,istate,r1,r2,theta)
          !
          ! eps - epsilon = ab initio energies - calculated pot. energies,
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
            ncol=0
            numpar = 0
            do  i=1,total_parameters
              if (ivar(i) .ne. 0) then
                 !
                 ncol=ncol+1
                 tempx=potparam(i)
                 deltax=fitfactordeltax*abs(tempx)
                 if (deltax .le. 1e-15) deltax=1e-5
                 !
                 potparam(i)=tempx+deltax
                 call poten_func(potparam,potright,istate,r1,r2,theta)
                 !
                 potparam(i)=tempx-deltax
                 !
                 call poten_func(potparam,potleft,istate,r1,r2,theta)
                 !
                 potparam(i)=tempx
                 rjacob(nrow+fitting%Nenergies,ncol)=(potright-potleft)/(2.0_rk*deltax)
                 !
              endif
            enddo ! --- ncol
            !
            numpar = ncol
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
           ncol=0
           do i=1,total_parameters
            if (ivar(i) > 0) then
                 ncol=ncol+1
                 potparam(i)=potparam(i)+dx(ncol)
              endif
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
           call write_potential_input_parameters(poten(1),'field1.fit')
           call write_potential_input_parameters(poten(2),'field2.fit')
           call write_potential_input_parameters(poten(3),'field3.fit')
           call write_potential_input_parameters(spinorbit(1),'field4.fit')
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
           ncol = 0 
           do i=1,total_parameters
              if (ivar(i) > 0) then
                  ncol=ncol+1
                 if (nused.eq.numpar) then  
                    sterr(ncol)=0
                 else
                    sterr(ncol)=sqrt(abs(ai(ncol,ncol)))*stadev
                    sum_sterr=sum_sterr+abs(sterr(ncol)/potparam(i))
                 endif
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
        do i=1,total_parameters
          write (f_out,"(a8,4x,i2,e22.14)") nampar(i),ivar(i),potparam(i)
        enddo
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
  

6550   format(/3X,67('-')/'   |  Iter  | Points | Params |   Deviat    |',&
       '    rms      | Stability |'/&
       3X,67('-')/,&
       '   | ',I6,' | ',I6,' | ',I6,' |  ',E12.5,' | ',E12.5,'  |  ',&
            E10.3,' |',/3X,67('-')/)


6551   format(/3X,67('-')/'   |  Iter  | Points | Params |   Deviat    |',&
       '    ssq1     |   ssq2    |'/&
       3X,67('-')/,&
       '   | ',I6,' | ',I6,' | ',I6,' |  ',E12.5,' | ',E12.5,'  |  ',&
            E10.3,' |',/3X,67('-')/)

6552   format(/3X,80('-')/'   |  Iter  | Points | Params |   Deviat    |',&
       '    ssq_ener |   ssq_pot | Stability |'/&
       3X,80('-')/,&
       '   | ',I6,' | ',I6,' | ',I6,' |  ',E12.5,' | ',E12.5,'  |  ',&
            E10.3,' |',E10.3,' |',/3X,80('-')/)

6553   format(/3X,80('-')/'   |  Iter  | Points | Params |   Deviat    |',&
       '    rms_ener |   rms_pot | Stability |'/&
       3X,80('-')/,&
       '   | ',I6,' | ',I6,' | ',I6,' |  ',E12.5,' | ',E12.5,'  |  ',&
            E10.3,' |',E10.3,' |',/3X,80('-')/)

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
            allocate(calc_(iref,irot,ipar_)%energy(Nsize),stat=alloc)
            call ArrayStart('calc_%energy',alloc,size(calc_(iref,irot,ipar_)%energy),kind(calc_(iref,irot,ipar_)%energy))
            !
          enddo
          !
        enddo 
      enddo
      !
      ncol=0
      !
      do  i=1,total_parameters
        if (ivar(i) > 0) then
          !
          ncol=ncol+1
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
          rjacob(1:fitting%Nenergies,ncol)=(enerright(1:fitting%Nenergies)-enerleft(1:fitting%Nenergies))/(2.0_rk*deltax)
          !
          rjacob(2:fitting%Nenergies,ncol) = rjacob(2:fitting%Nenergies,ncol) - rjacob(1,ncol)
          !
          potparam(i)=tempx
          !
          call map_parameters(dir=.false.)
          !
        endif
        !
      enddo ! --- ncol
      !
      deallocate (enerright,enerleft)
      !
      call ArrayStop('enerright')
      call ArrayStop('enerleft')
      !
      do iref = 1,2
        do irot = 1,Nrot
          do ipar_ = 1,2
            !
            if (associated(calc_(iref,irot,ipar_)%energy)) deallocate(calc_(iref,irot,ipar_)%energy)
            !                                                                
          enddo                                                              
        enddo                                                                
      enddo
      !
      call ArrayStop('calc_%energy')                                                               
      !                                                                      
    end subroutine finite_diff_of_energy                                     
    !
    subroutine map_parameters(dir)
     !
     logical,intent(in) :: dir
     !
     type(fieldT),pointer      :: field
     !
     ! mapping 
     !
     if (dir) then 
       !
       field => poten(1) 
       potparam(1:field%Nterms) = field%value(1:field%Nterms)
       nampar(1:field%Nterms)   = field%forcename(1:field%Nterms)
       ivar(1:field%Nterms)     = int(field%weight(1:field%Nterms))
       iparams = field%Nterms
       !
       field => poten(2) 
       potparam(iparams+1:iparams+field%Nterms) = field%value(1:field%Nterms)
       nampar(iparams+1:iparams+field%Nterms)   = field%forcename(1:field%Nterms)
       ivar(iparams+1:iparams+field%Nterms)     = int(field%weight(1:field%Nterms))
       iparams = iparams+field%Nterms
       !
       field => poten(3) 
       potparam(iparams+1:iparams+field%Nterms) = field%value(1:field%Nterms)
       nampar(iparams+1:iparams+field%Nterms)   = field%forcename(1:field%Nterms)
       ivar(iparams+1:iparams+field%Nterms)     = int(field%weight(1:field%Nterms))
       iparams = iparams+field%Nterms
       !
       field => spinorbit(1) 
       potparam(iparams+1:iparams+field%Nterms) = field%value(1:field%Nterms)
       nampar(iparams+1:iparams+field%Nterms)   = field%forcename(1:field%Nterms)
       ivar(iparams+1:iparams+field%Nterms)     = int(field%weight(1:field%Nterms))
       iparams = iparams+field%Nterms
       !
     else
       !
       field => poten(1) 
       field%value(1:field%Nterms) = potparam(1:field%Nterms)
       field%weight(1:field%Nterms) = real(ivar(1:field%Nterms),rk)
       iparams = field%Nterms
       !
       field => poten(2) 
       field%value(1:field%Nterms) = potparam(iparams+1:iparams+field%Nterms)
       field%weight(1:field%Nterms) = real(ivar(iparams+1:iparams+field%Nterms),rk)
       iparams = iparams+field%Nterms
       !
       field => poten(3) 
       field%value(1:field%Nterms) = potparam(iparams+1:iparams+field%Nterms)
       field%weight(1:field%Nterms) = real(ivar(iparams+1:iparams+field%Nterms),rk)
       iparams = iparams+field%Nterms
       !
       field => spinorbit(1) 
       potparam(iparams+1:iparams+field%Nterms) = field%value(1:field%Nterms)
       field%weight(1:field%Nterms) = real(ivar(iparams+1:iparams+field%Nterms),rk)
       !
     endif 
     !
   end subroutine map_parameters
   !
   
   subroutine run_everest
      !
      logical  ::  SYSTEMQQ  ! system function for  calling extyernal programs, 
                             ! platform dependent.  
      !
      call map_parameters(dir=.false.)
      call write_potential_input_parameters(poten(1),'field1.fit')
      call write_potential_input_parameters(poten(2),'field2.fit')
      call write_potential_input_parameters(poten(3),'field3.fit')
      call write_potential_input_parameters(spinorbit(1),'field4.fit')
      !
      if (verbose>=3) write(f_out,"(/'calling Everest program for J = ',f9.1,', iparity =  ',i2'...')") jrot,iparity
      !
#if (debug_ == 0)
         isys = systemqq(evib_exe//'< evib.inp >> evib.out')
#endif
      !
      if (verbose>=3) write(f_out,"('callling rotlev program ...')") 
      !
#if (debug_ == 0)
         isys = systemqq(erot_exe//'<erot.inp >>erot.out')
         isys = systemqq('./collect_energies_rot.sh erot.out > erot.log')
#endif
      !
   end subroutine run_everest
   !

  end subroutine fitting_energies
  !

  subroutine write_potential_input_parameters(field,filename)
    !
    type(fieldT),intent(in)  :: field ! field1.fit
    character(len=10) :: filename
    !
    integer(ik) :: f_p = 21,i
      !
      open(f_p,file=trim(filename),status='replace')
      !
      write(f_p,*) field%Nterms
      write(f_p,"('r1',4x,f19.10)") field%r1
      write(f_p,"('r2',4x,f19.10)") field%r2
      write(f_p,"('alpha',4x,f19.10)") field%alpha
      !
      do i=1,field%Nterms
         write(f_p,"(a8,1x,3i4,1x,e24.12)") adjustl((field%forcename(i))),field%ipower1(i),field%ipower2(i),field%ipower3(i),&
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
    end subroutine write_potential_input_parameters
  

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
      integer(ik)    :: alloc,i,irot,iref,v1,v2,v3,ka
      real(rk)  :: Energy,De,wei,Ome,Wome,Wma,Jrot
      character(len=5)   :: Jch
      integer(ik)  :: Nstates_(2,0:Nmax,2)
      !
      tunit = 22
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
         if (associated(calc(1,0,1)%energy)) then 
           !
           do iref  = 1,2
             do irot = 0,Nmax
               do ipar = 1,2
                   deallocate(calc(iref,irot,ipar)%energy)
                   deallocate(calc(iref,irot,ipar)%quanta)
               enddo
             enddo
           enddo
           !
           call ArrayStop('calc%energy')
           call ArrayStop('calc%quanta')
           !
         endif
         !
         do iref  = 1,2
           do irot = 0,Nmax
             do ipar = 1,2
               !
               Nsize = max(1,Nstates(iref,irot,ipar))
               !
               allocate(calc(iref,irot,ipar)%energy(Nsize),calc(iref,irot,ipar)%quanta(Nsize),stat=alloc)
               call ArrayStart('calc%energy',alloc,size(calc(iref,irot,ipar)%energy),kind(calc(iref,irot,ipar)%energy))
               call ArrayStart('calc%quanta',alloc,size(calc(iref,irot,ipar)%quanta),kind(calc(iref,irot,ipar)%quanta))
               !
             enddo
           enddo 
         enddo
         !
         rewind(tunit)
         !
      endif 
      !
      Nstates_ = 0
      !
      do 
        read(tunit,"(a5,2x,i2,1x,i5,2(1x,f12.4),1x,i3,1x,2(1x,i3),3(1x,i2),1x,f5.3,1x,f3.1,1x,f5.3,1x,i2)",end=121) &
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
        calc(iref,Nrot,ipar)%energy(nn) = Energy
        calc(iref,Nrot,ipar)%quanta(nn)%v1 = v1
        calc(iref,Nrot,ipar)%quanta(nn)%v2 = v2
        calc(iref,Nrot,ipar)%quanta(nn)%v3 = v3
        calc(iref,Nrot,ipar)%quanta(nn)%istate = Pot
        calc(iref,Nrot,ipar)%quanta(nn)%Omega = Ome
        calc(iref,Nrot,ipar)%quanta(nn)%Ka = Ka
        !
        cycle
        !         
      121 continue
          exit

      end do      
      !
      close(tunit,status='keep')
      !
    end subroutine get_vibronic_energies
    !

    !
    subroutine create_the_xpect_job_file(nparams,ivar)
	  !
      integer,intent(in)  :: nparams,ivar(nparams)
      integer             :: lpot,npropin,nprt,nv1,i
      !
      open(1,file='xpect.inp')
      !
      write(1,"('&PRT zform=.true., /')")
      write(1,"('Expectation values, fort.11,properties=3')")
      !nvib = neval_
      !
      lpot = nalf_ ! 38
      npropin = nparams
      nprt = 0
      nv1 = 0
      !
      write(1,"(4i5)") lpot,npropin,nprt,nv1
      do i = 1,npropin 
        !
        if (ivar(i)>0) then
          !
          write(1,"(i5)") 1
          !
        else
          !
          write(1,"(i5)") 0
          !    
        endif
        !
      enddo
      !
      close(1,status='keep')
      !
    end subroutine create_the_xpect_job_file
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
      real(8),intent(in)    :: sigma(:),eps(:)
      real(8),intent(inout) :: wt(:)
      !
      integer            :: npts,i,nrow,nused
      real               :: da1,da2,wtsum,da
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
         if (wt(nrow)>0) wt(nrow) = 1.0d0/( sigma(nrow)**2 + eps(nrow)**2 )
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


