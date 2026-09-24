// The two gfortran runtime entry points the k-space kernels can reach (allocation failure, a runtime check). With these
// here the library needs no Fortran runtime (libgfortran, libquadmath) on any platform: the kernels do no I/O and use
// no other library routine. Both report and stop, as gfortran's own would.
#include <cstdarg>
#include <cstdio>
#include <cstdlib>

extern "C" {

[[noreturn]] void _gfortran_os_error_at(const char* where, const char* message, ...) {
  std::fprintf(stderr, "CAPS (Fortran kernel) %s: ", where ? where : "");
  va_list ap;
  va_start(ap, message);
  std::vfprintf(stderr, message, ap);
  va_end(ap);
  std::fprintf(stderr, "\n");
  std::abort();
}

[[noreturn]] void _gfortran_runtime_error(const char* message, ...) {
  std::fprintf(stderr, "CAPS (Fortran kernel): ");
  va_list ap;
  va_start(ap, message);
  std::vfprintf(stderr, message, ap);
  va_end(ap);
  std::fprintf(stderr, "\n");
  std::abort();
}

}
