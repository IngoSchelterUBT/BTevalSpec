subroutine padeSeries(w, wn, n, dt, a0, b0, a, b, mu)
  implicit none
  !------------------------------------------------------------------------------!
  !IO variables
  double precision, dimension(:), intent(in) :: w
  integer, intent(in) :: wn, n
  double precision, intent(in) :: dt, a0, b0
  double precision, dimension(:), intent(in) :: a, b !pade coefficients
  double complex, dimension(1:wn), intent(out) :: mu

  !------------------------------------------------------------------------------!
  !Local variables
  integer :: wi, k, i
  double complex :: z0, z1
  double complex, parameter :: imag = dcmplx(0.0,1.d0)
  double complex, dimension(1:wn) :: p, q
  double complex, dimension(1:n) :: z


  !------------------------------------------------------------------------------!
  !Pade approximation for all values of w
  do wi = 1,wn
    !----------------------------------------------------------------------------!
    !Get powers of the complex exponential function
    z0 = 1.d0
    z1 = exp(imag*w(wi)*dt)
    !z1 = exp(-imag*w(wi)*dt)
    z(1) = z1
    do i=2,n
      z(i) = z(i-1)*z1
    end do
  
    !----------------------------------------------------------------------------!
    !Create Pade approximation eq~(29) in [1]
    p(wi) = a0
    q(wi) = b0
    do k=1,n
      p(wi) = p(wi)+a(k)*z(k)
      q(wi) = q(wi)+b(k)*z(k)
    end do
    mu(wi) = p(wi)/q(wi)
  end do

end subroutine
