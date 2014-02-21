// Some utils and wrappers used in the ninja-library.


#ifndef NINJA_TMP_UTILS_HH
#define NINJA_TMP_UTILS_HH

#include <ninja/types.hh>


#if HAVE_DECLTYPE

// If decltype() is available, use the C++11 std::common_type to
// choose return type in mixed types operations.

#include <type_traits>

namespace ninja{
  namespace details{

    template<typename T, typename U>
    struct common_type {
      typedef typename std::common_type<T,U>::type type;
    };

  }
}

#else  // HAVE_DECLTYPE

// If decltype() is not available, use some template structs in order
// to choose return type in mixed types operations.

#include <ninja/types.hh>

namespace ninja {
  namespace details {

    // Choose return type in mixed floating point operations
    template<typename T, typename U>
    struct common_type;
    template<> struct common_type<ninja::Complex, ninja::Real> {
      typedef ninja::Complex type;
    };
    template<>
    struct common_type<ninja::Real,ninja::Complex> {
      typedef ninja::Complex type;
    };
    template<>
    struct common_type<ninja::Complex,ninja::Complex> {
      typedef ninja::Complex type;
    };
    template<typename T>
    struct common_type<T,T> {
      typedef T type;
    };

  } // namespace details
} // namespace ninja

#endif // ! HAVE_DECLTYPE

#endif // NINJA_TMP_UTILS_HH
