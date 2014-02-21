// This -*- C++ -*- header file contains the definition of the classes
// ZeroFloat and ZeroFloatArray.
//
// ZeroFloat is an empty class which behaves like the floating point
// numer 0.0.  The compiler should be able to automatically set to
// zero some (sub)expressions containing the ZeroFloat type at compile
// time.  It is used for the Masslesss type instantation of the
// Amplitde methods.
//
// ZeroFloatArray is an empty class which behaves like an array of
// ZeroFloat objects.


#ifndef ZERO_FLOAT_HH
#define ZERO_FLOAT_HH

#include <iostream>

namespace ninja {

  class ZeroFloat {
  public:
    ZeroFloat() {}

    // +0 = 0
    ZeroFloat operator + () const
    {
      return ZeroFloat();
    }
    
    // -0 = 0
    ZeroFloat operator - () const
    {
      return ZeroFloat();
    }
    
  private:
    // no-data
  };

  // 0+0 = 0
  // x+0 = 0
  // 0+x = 0
  inline ZeroFloat operator + (ZeroFloat, ZeroFloat)
  {
    return ZeroFloat();
  }
  template<typename T>
  inline T operator + (const T & x, ZeroFloat)
  {
    return x;
  }
  template<typename T>
  inline T operator + (ZeroFloat, const T & x)
  {
    return x;
  }

  // 0-0 = 0
  // x-0 = x
  // 0-x = -x
  inline ZeroFloat operator - (ZeroFloat, ZeroFloat)
  {
    return ZeroFloat();
  }
  template<typename T>
  inline T operator - (const T & x, ZeroFloat)
  {
    return x;
  }
  template<typename T>
  inline T operator - (ZeroFloat, const T & x)
  {
    return -x;
  }

  // 0*0 = 0
  // 0*x = 0
  // x*0 = 0
  inline ZeroFloat operator * (ZeroFloat, ZeroFloat)
  {
    return ZeroFloat();
  }
  template<typename T>
  inline ZeroFloat operator * (const T &, ZeroFloat)
  {
    return ZeroFloat();
  }
  template<typename T>
  inline ZeroFloat operator * (ZeroFloat, const T &)
  {
    return ZeroFloat();
  }

  // 0/x = 0
  template<typename T>
  inline ZeroFloat operator / (ZeroFloat, const T &)
  {
    return ZeroFloat();
  }

  // Comparosions with 0
  template<typename T>
  inline bool operator == (ZeroFloat, const T & x)
  {
    return (x==T(0));
  }
  template<typename T>
  inline bool operator == (const T & x, ZeroFloat)
  {
    return (x==T(0));
  }
  template<typename T>
  inline bool operator != (ZeroFloat, const T & x)
  {
    return (x!=T(0));
  }
  template<typename T>
  inline bool operator != (const T & x, ZeroFloat)
  {
    return (x!=T(0));
  }
  template<typename T>
  inline bool operator < (ZeroFloat, const T & x)
  {
    return (x>T(0));
  }
  template<typename T>
  inline bool operator < (const T & x, ZeroFloat)
  {
    return (x<T(0));
  }

  // Standard output stream
  inline std::ostream & operator << (std::ostream & os, ZeroFloat)
  {
    return os << "(zero)";
  }


  // The following class can be used in order to mimic an array of zeroes
  class ZeroFloatArray {

  public:

    ZeroFloatArray() {}

    template<typename T> ZeroFloatArray(const T *) {}

    ZeroFloat operator [] (int)
    {
      return ZeroFloat();
    }

  private:
    // no-data
  };

  
  // Specialize ninja::const_pointer<ZeroFloat> (defined in ninja/types.hh)
  template<typename X> struct const_pointer;
  template<> struct const_pointer<ZeroFloat> {
    typedef ZeroFloatArray type;
  };


  namespace details {

    // Specialize ninja::details::common_type<...> (defined in
    // tmp_utils.hh internal header)
    template<typename T, typename U> struct common_type;
    template<typename T>
    struct common_type<ZeroFloat, T> {
      typedef T type;
    };
    template<typename T>
    struct common_type<T,ZeroFloat> {
      typedef T type;
    };
    template<>
    struct common_type<ZeroFloat,ZeroFloat> {
      typedef ZeroFloat type;
    };


    // compile-time check: if a type is massless
    template <typename MassType>
    struct MasslessTypeError {};
    template <>
    struct MasslessTypeError<ZeroFloat> {
      static void MassType_must_be_massless() {}
    };

  } // namespace details

} // namespace ninja

#endif // ZERO_FLOAT_HH
