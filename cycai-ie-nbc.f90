! chang2013 in 2D with import export
! Integration using pseudospectral methods
! The integration lattice is n x m which is doubled for Neumann boundatry conditions.
! 2d FFTs are calculated calling the DFT routines of Intel Math Kernel Library.
! In real space the field is stored in a Real(8) matrix of dimension 2n x 2m.
! The field in Fourier space is in CCE format, which is a Double Complex matrix of dimension (n+1)*m.
!    Columns 1 ... n correspond to qy = 0 ... (n-1)*delta_ky and columns n+1...n to qy=-n*delta_ky ... -delta_ky
!    Rows 1 ... n+1 correspond to qx = 0 ... n*delta_kx
!    Elements in row 1 do not have imaginary part, but the program does not make use of this fact.
!    For n even (usual case) elements in row n+1 does not have imaginary part (the program does not make use of this).
!
! This version allows for complex propagators in Fourier space.
!
! Written by Daniel Ruiz-Reynés 
! Version: 19 June 2024
Use MKL_DFTI
implicit none
include 'dim.h'


!Linear and nonlinear propagators in Fourier space for delta_t and 2 delta_t
double complex, dimension (n+1,2*m) :: l_prop,nl_prop,l_prop2,nl_prop2
double complex, dimension (n+1,2*m) :: l_prop_S,nl_prop_S,l_prop2_S,nl_prop2_S

!2DFields in Fourier space and equivalent 1d vectors needed for MKL FFTs
double complex, dimension (n+1,2*m) :: e_fourier,enl_fourier,field_fourier
double complex, dimension ((n+1)*2*m) :: e_fourier_equiv,enl_fourier_equiv,field_fourier_equiv
double complex, dimension (n+1,2*m) :: S_fourier,Snl_fourier
double complex, dimension ((n+1)*2*m) :: S_fourier_equiv,Snl_fourier_equiv

equivalence (e_fourier,e_fourier_equiv)
equivalence (enl_fourier,enl_fourier_equiv)
equivalence (S_fourier,S_fourier_equiv)
equivalence (Snl_fourier,Snl_fourier_equiv)
equivalence (field_fourier,field_fourier_equiv)

!2D Fields in real space and equivelent 1d vectors needed for MKL FFTs
real(8), dimension (2*n,2*m)  :: e,S,agrad,field, derv, derv2
real(8), dimension (n,m)  :: half_e, half_S
real(8), dimension (2*n*2*m) :: e_equiv,S_equiv,field_equiv

equivalence (e,e_equiv)
equivalence (S,S_equiv)
equivalence (field,field_equiv)

real(8), dimension (2*n) :: x
!Wavevectors
real(8), dimension (2*m) :: qy
real(8), dimension (n+1) :: qx

!Nuclei parameters
integer, dimension(2*n_nuclei) :: i_i
real(8), dimension(2*n_nuclei) :: x_i,eps_i,s_i
!Model and integration parameters
character(len=32) :: parameter_file,ICn_file,ICS_file,FCn_file,FCS_file,n_file,S_file,empty,nuclei_file
real(8) :: ks,a_deg,b_deg,ec50_deg,n_deg,K_deg 
real(8) :: a_cdc25,b_cdc25,ec50_cdc25,n_cdc25,K_cdc25,K_neb
real(8) :: a_wee1,b_wee1,ec50_wee1,n_wee1,K_wee1
real(8) :: Da,Di
real(8) :: ec50_neb,n_neb
real(8) :: delta_x,delta_y,delta_t,t_write,t_final,noise
integer :: seed
namelist /model_params/ ks,a_cdc25,b_cdc25,ec50_cdc25,n_cdc25,a_wee1,b_wee1,ec50_wee1,n_wee1,a_deg,b_deg,ec50_deg,n_deg,Da,Di
namelist /nuclear_params/ ec50_neb,n_neb
namelist /integration_params/ delta_x,delta_y,delta_t,t_write,t_final
namelist /init_params/ noise,seed

!Auxiliar variables
double complex :: alfa,alfa_S
integer :: l,j,k,i,n_write,idx_write
real(8) :: delta_kx,delta_ky,dran_g,q2,pi,nonlinear_scale
real(8), dimension(n+1,2*m) :: nablas

!Auxiliar variables for fft using mkl.
type(DFTI_DESCRIPTOR), POINTER :: punt_f,punt_b
integer :: status
integer, dimension (2) :: matrix_size
integer, dimension (3) ::  strides_in,strides_out


!Get 5 arguments of execution.
empty = ''
call getarg(1,parameter_file)
write(*,*) parameter_file

call getarg(2,ICn_file)
write(*,*) ICn_file
call getarg(3,ICS_file)
write(*,*) ICS_file

call getarg(4,nuclei_file)
write(*,*) nuclei_file

call getarg(5,n_file)
if (n_file.eq.empty) n_file = 'a.dat'
write(*,*) n_file
call getarg(6,S_file)
if (S_file.eq.empty) S_file = 'i.dat'
write(*,*) S_file

call getarg(7,FCn_file)
if (FCn_file.eq.empty) FCn_file = 'fa.dat'
write(*,*) FCn_file
call getarg(8,FCS_file)
if (FCS_file.eq.empty) FCS_file = 'fi.dat'
write(*,*) FCS_file


!******** Reading input parameters *****************************************
open (unit=10,file=parameter_file,status='old')
read(10,nml=model_params)
read(10,nml=nuclear_params)
read(10,nml=integration_params)
read(10,nml=init_params)
close (unit=10)

print*, 'ks    =',ks
print*, 'a_deg,b_deg,ec50_deg,n_deg    =',a_deg,b_deg,ec50_deg,n_deg
print*, 'a_cdc25,b_cdc25,ec50_cdc25,n_cdc25    =',a_cdc25,b_cdc25,ec50_cdc25,n_cdc25
print*, 'a_wee1,b_wee1,ec50_wee1,n_wee1    =',a_wee1,b_wee1,ec50_wee1,n_wee1
print*, 'Da    =',Da
print*, 'Di    =',Di
print*, 'N_nuclei = ', N_nuclei
print*, 'ec50_neb,n_neb = ', ec50_neb,n_neb
print*, 'dx       =',delta_x
print*, 'dy       =',delta_y
print*, 'dt       =',delta_t
print*, 'tf       =',t_final



!******** Setting propagators for integration in Fourier space *************
!Setting q vectors.
pi=dacos(-1.d0)
delta_kx = pi/(n*delta_x)
delta_ky = pi/(m*delta_y)
do k=1,n
  qx(k) = (k-1)*delta_kx
enddo
qx(n+1) = -n*delta_kx

do k=1,m
  qy(k) = (k-1)*delta_ky
  qy(k+m) = (-m+k-1)*delta_ky
enddo

print*, 'dkx       =',delta_kx
print*, 'dky       =',delta_ky

!Setting propagators
!eq and eq2: propagators linear term over a time step delta_t and 2*delta_t
!feq and feq2: propagators nonlinear term over delta_t and 2*delta_t scaled with the strengh of nonlinear term and Fourier transform normalization.
nonlinear_scale = 1.d0/(2*n*2*m)
do j = 1, n+1
  do k = 1, 2*m
    q2 = qx(j)*qx(j) + qy(k)*qy(k)
    nablas(j,k) = -q2
    alfa = dcmplx(+ Da*q2 ,0.d0)
    alfa_S = dcmplx(+ Di*q2,0.d0)
    if ((cdabs(alfa).eq.0.d0)) then
      l_prop(j,k)  = dcmplx(1.d0,0.d0)
      nl_prop(j,k) = dcmplx(nonlinear_scale*delta_t,0.d0)
      l_prop2(j,k) = dcmplx(1.d0,0.d0)
      nl_prop2(j,k)= dcmplx(2.d0*nonlinear_scale*delta_t,0.d0)
    else
      l_prop(j,k)  = cdexp(-alfa*delta_t)
      nl_prop(j,k) = nonlinear_scale*(1.d0-l_prop(j,k))/alfa
      l_prop2(j,k) = cdexp(-2.d0*alfa*delta_t)
      nl_prop2(j,k)= nonlinear_scale*(1.d0-l_prop2(j,k))/alfa
    endif
    if ((cdabs(alfa_S).eq.0.d0)) then
      l_prop_S(j,k)  = dcmplx(1.d0,0.d0)
      nl_prop_S(j,k) = dcmplx(nonlinear_scale*delta_t,0.d0)
      l_prop2_S(j,k) = dcmplx(1.d0,0.d0)
      nl_prop2_S(j,k)= dcmplx(2.d0*nonlinear_scale*delta_t,0.d0)
    else
      l_prop_S(j,k)  = cdexp(-alfa_S*delta_t)
      nl_prop_S(j,k) = nonlinear_scale*(1.d0-l_prop_S(j,k))/alfa_S
      l_prop2_S(j,k) = cdexp(-2.d0*alfa_S*delta_t)
      nl_prop2_S(j,k)= nonlinear_scale*(1.d0-l_prop2_S(j,k))/alfa_S
    endif
  enddo
enddo
!******** Neumann Boundary Conditions ********************************
l_prop = dcmplx(real(l_prop),0d0)
nl_prop = dcmplx(real(nl_prop),0d0)
l_prop2 = dcmplx(real(l_prop2),0d0)
nl_prop2 = dcmplx(real(nl_prop2),0d0)
l_prop_S = dcmplx(real(l_prop_S),0d0)
nl_prop_S = dcmplx(real(nl_prop_S),0d0)
l_prop2_S = dcmplx(real(l_prop2_S),0d0)
nl_prop2_S = dcmplx(real(nl_prop2_S),0d0)


!******** Seting handles for FFT transforms ********************************
matrix_size(1)= 2*n  !Matrix size in real space in x direction.
matrix_size(2)= 2*m  !Matrix size in real space in y direction.
strides_in(1) = 0 !Strides for the transform from real space into fourier space.
strides_in(2) = 1
strides_in(3) = 2*n
strides_out(1) = 0     !Strides for the transform from Fourier space into real space.
strides_out(2) = 1
strides_out(3) = n+1
!Creating a pointer for the forward transform from nxn real matrix into a n+1x2n complex matrix
status = DftiCreateDescriptor(punt_f, DFTI_DOUBLE, DFTI_REAL, 2, matrix_size)
status = DftiSetValue(punt_f, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
status = DftiSetValue(punt_f, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
status = DftiSetValue(punt_f, DFTI_INPUT_STRIDES, strides_in)
status = DftiSetValue(punt_f, DFTI_OUTPUT_STRIDES, strides_out)
status = DftiCommitDescriptor(punt_f)
!Creating a pointer for the backward transform. It is the same as for the forward but the input and output strides have to be reversed.
!See mkl/examples/dftfreal_2d_cce_double_ex2.f90
status = DftiCopyDescriptor (punt_f,punt_b)
status = DftiSetValue(punt_b, DFTI_INPUT_STRIDES, strides_out)
status = DftiSetValue(punt_b, DFTI_OUTPUT_STRIDES, strides_in)
status = DftiCommitDescriptor(punt_b)

!******** Initialization of the field *************************************
call dran_ini(seed) !Initialization of random number generator.

open(unit=20,file=ICn_file,status='unknown',form='unformatted')
read(20) half_e
close(20)
open(unit=20,file=ICS_file,status='unknown',form='unformatted')
read(20) half_S
close(20)
open(unit=30,file=nuclei_file,status='unknown',form='unformatted')
read(30) x_i(1:n_nuclei)
read(30) eps_i(1:n_nuclei)
read(30) s_i(1:n_nuclei)
close(unit=30)

do k=1,m
  do j=1,n
   half_e(j,k) = half_e(j,k)+noise*dran_g()
   half_S(j,k) = half_S(j,k)+noise*dran_g()
  enddo
enddo

do k=1,m ! Mirror ICs
  do j=1,n
   e(j,k) = half_e(j,k)
   e(2*n -(j-1),k) = half_e(j,k)
   e(j,2*m - (k-1)) = half_e(j,k)
   e(2*n -(j-1),2*m - (k-1)) = half_e(j,k)
   S(j,k) = half_S(j,k)
   S(2*n -(j-1),k) = half_S(j,k)
   S(j,2*m - (k-1)) = half_S(j,k)
   S(2*n -(j-1),2*m - (k-1)) = half_S(j,k)
  enddo
enddo

do j=1,n ! Mirror x
  x(j) = (j-1)*delta_x
  x(2*n -(j-1)) = (j-1)*delta_x
enddo
i_i(1:n_nuclei) = nint(x_i(1:n_nuclei)/delta_x) + 1 
do j = 1, n_nuclei ! Mirror i_i, eps_i, s_i
  i_i(2*n_nuclei -(j-1)) = i_i(j)
  eps_i(2*n_nuclei -(j-1)) = eps_i(j)
  s_i(2*n_nuclei -(j-1)) = s_i(j)
enddo

K_deg = ec50_deg**n_deg
K_cdc25 = ec50_cdc25**n_cdc25
K_wee1 = ec50_wee1**n_wee1
K_neb = ec50_neb**n_neb

!Opening data file to store field evoltion and writing initial condition.
open(unit=40,file=n_file,status='unknown',form='unformatted')
write(40) e(1:n,1:m)
open(unit=50,file=s_file,status='unknown',form='unformatted')
write(50) S(1:n,1:m)
write(*,*) 0.d0,e(n/2,m/2),S(n/2,m/2)!writing field at central point on screen.

!Compute the initial field in Fourier space and scaling it to 1/n^2
status = DftiComputeForward(punt_f,e_equiv,e_fourier_equiv)
e_fourier = e_fourier/(2*n*2*m)
status = DftiComputeForward(punt_f,S_equiv,S_fourier_equiv)
S_fourier = S_fourier/(2*n*2*m)
 
!******** Integration ******************************************************
n_write = nint(t_write/(2.0*delta_t))

do idx_write = 1, nint(t_final/t_write)  !Loop over times in which writing is written
    
  do i = 1, n_write !Integration between times in which data is written  

    !Compute nonlinear part e
    forall (j=1:n+1) field_fourier(j,:) = dcmplx(0,1d0)*qx(j)*e_fourier(j,:)
    status = DftiComputeBackward(punt_b,field_fourier_equiv,field_equiv)

    call derivatives_of_v(e,x,K_neb,n_neb,i_i,eps_i,s_i, derv,derv2)

    field = ks + (a_cdc25 + b_cdc25*e**n_cdc25/(K_cdc25 + e**n_cdc25))*S &
          - (a_wee1 + b_wee1*K_wee1/(K_wee1 + e**n_wee1))*e  &
          - (a_deg + b_deg*e**n_deg/(K_deg + e**n_deg))*e &
          + field*derv + e*derv2

    status = DftiComputeForward(punt_f,field_equiv,field_fourier_equiv)    !Transform the nonlinear terms into Fourier space.
    !Update e dt in fourier space
    enl_fourier = l_prop*e_fourier + nl_prop*field_fourier  !Propagation of the field in Fourier space over a time step delta_t

    !Compute nonlinear part S
    forall (j=1:n+1) field_fourier(j,:) = dcmplx(0,1d0)*qx(j)*S_fourier(j,:)
    status = DftiComputeBackward(punt_b,field_fourier_equiv,field_equiv)

    call derivatives_of_v(S,x,K_neb,n_neb,i_i,eps_i,s_i, derv,derv2)

    field = - (a_cdc25 + b_cdc25*e**n_cdc25/(K_cdc25 + e**n_cdc25))*S &
          + (a_wee1 + b_wee1*K_wee1/(K_wee1 + e**n_wee1))*e  &
          - (a_deg + b_deg*e**n_deg/(K_deg + e**n_deg))*S &
          + field*derv + S*derv2
    status = DftiComputeForward(punt_f,field_equiv,field_fourier_equiv)    !Transform the nonlinear terms into Fourier space.
    !update S dt in fourier space
    Snl_fourier = l_prop_S*S_fourier + nl_prop_S*field_fourier  !Propagation of the field in Fourier space over a time step delta_t

    !Update e,S real space
    status = DftiComputeBackward(punt_b,enl_fourier_equiv,e_equiv)
    status = DftiComputeBackward(punt_b,Snl_fourier_equiv,S_equiv)

    ! -------- Compute field for 2dt ----------

    !Compute nonlinear part e
    forall (j=1:n+1) field_fourier(j,:) = dcmplx(0,1d0)*qx(j)*enl_fourier(j,:)
    status = DftiComputeBackward(punt_b,field_fourier_equiv,field_equiv)

    call derivatives_of_v(e,x,K_neb,n_neb,i_i,eps_i,s_i, derv,derv2)

    field = ks + (a_cdc25 + b_cdc25*e**n_cdc25/(K_cdc25 + e**n_cdc25))*S &
          - (a_wee1 + b_wee1*K_wee1/(K_wee1 + e**n_wee1))*e  &
          - (a_deg + b_deg*e**n_deg/(K_deg + e**n_deg))*e &
          + field*derv + e*derv2
    status = DftiComputeForward(punt_f,field_equiv,field_fourier_equiv)    !Transform the nonlinear terms into Fourier space.
    !Update e 2dt in fourier space
    e_fourier = l_prop2*e_fourier + nl_prop2*field_fourier !Propagation of the field in Fourier space over a time step 2*delta_t

    !Compute nonlinear part S
    forall (j=1:n+1) field_fourier(j,:) = dcmplx(0,1d0)*qx(j)*Snl_fourier(j,:)
    status = DftiComputeBackward(punt_b,field_fourier_equiv,field_equiv)

    call derivatives_of_v(S,x,K_neb,n_neb,i_i,eps_i,s_i, derv,derv2)

    field = - (a_cdc25 + b_cdc25*e**n_cdc25/(K_cdc25 + e**n_cdc25))*S &
          + (a_wee1 + b_wee1*K_wee1/(K_wee1 + e**n_wee1))*e  &
          - (a_deg + b_deg*e**n_deg/(K_deg + e**n_deg))*S &
          + field*derv + S*derv2
    status = DftiComputeForward(punt_f,field_equiv,field_fourier_equiv)    !Transform the nonlinear terms into Fourier space.
    !Update S 2dt in fourier space
    S_fourier = l_prop2_S*S_fourier + nl_prop2_S*field_fourier  !Propagation of the field in Fourier space over a time step delta_t
    !Update e,S real space
    status = DftiComputeBackward(punt_b,e_fourier_equiv,e_equiv)
    status = DftiComputeBackward(punt_b,S_fourier_equiv,S_equiv)
    !Update e_fourier to remove acumulated error (otherwise leads to instabilities)
    status = DftiComputeForward(punt_f,e_equiv,e_fourier_equiv)
    e_fourier = e_fourier/(2*n*2*m)
    status = DftiComputeForward(punt_f,S_equiv,S_fourier_equiv)
    S_fourier = S_fourier/(2*n*2*m)

  enddo
  
  write(*,*) idx_write*t_write,e(n/2,m/2),S(n/2,m/2) !writing field at central point on screen.
  write(40) e(1:n,1:m)
  write(50) S(1:n,1:m)
enddo

!******** Finishing and storing final data  ********************************
close (unit=40) !Closing file where time field evolution is stored.
close (unit=50)
!Storing the final status of the field in double precision.

open(unit=60,file=FCn_file,status='unknown',form='unformatted')
write(60) e(1:n,1:m)
close(unit=60)
open(unit=70,file=FCS_file,status='unknown',form='unformatted')
write(70) S(1:n,1:m)
close(unit=70)

end

!-----------------------------------------------------------------------------------------
subroutine derivatives_of_v(field,x,K_neb,n_neb,i_i,eps_i,s_i, derv,derv2)
  include 'dim.h'
  integer :: i,j
  real(8), dimension (2*n,2*m) :: field, derv, derv2
  integer, dimension(2*n_nuclei) :: i_i
  real(8), dimension(2*n_nuclei) :: eps_i,s_i
  real(8) :: K_neb,n_neb,f_i
  real(8), dimension (2*n) :: x

  derv = 0d0
  derv2 = 0d0
  do j = 1, 2*n_nuclei
    f_i =  K_neb/(K_neb + field(i_i(j),0)**n_neb)
    derv(:,1)  = derv(:,1)  + f_i*eps_i(j)*dexp(-((x-x(i_i(j)))/s_i(j))**2/2)  &
    *(x-x(i_i(j)))/(s_i(j)**2)
    derv2(:,1) = derv2(:,1) + f_i*eps_i(j)*dexp(-((x-x(i_i(j)))/s_i(j))**2/2)  &
    *(1/(s_i(j)**2) -(x-x(i_i(j)))**2/(s_i(j)**4))
  enddo
  do k = 2, m
    derv(:,k) =  derv(:,1)
    derv(:,k) =  derv(:,1)
  enddo
  return
end
