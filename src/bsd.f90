module bsd
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, bsd!"
  end subroutine say_hello
end module bsd
