// -*- C++ -*- header file included in the source files generated by
// the script ninjanumgen.py.  It is supposed to work around some
// issue with the evaluation of the output of Form in C++.

#include <ninja/num_defs.hh>

namespace {

  typedef ninja::Complex NinjaAbbrType;

  //template<unsigned N>
  //struct NinjaFormAbbr {
  //    ninja::Complex data_[N];
  //    ninja::Complex & operator[](unsigned i)
  //    {
  //      return data_[i-1];
  //    }
  //    ninja::Complex & operator[](unsigned i) const
  //    {
  //      return data_[i-1];
  //    }
  //};

  // Integer powers.  Inline implementation, avoids calls to pow
  // function.
  template <typename T>
  inline T powi1(const T & x)
  {
    return x;
  }
  template <typename T>
  inline T powi2(const T & x)
  {
    return x*x;
  }
  template <typename T>
  inline T powi3(const T & x)
  {
    return x*x*x;
  }
  template <typename T>
  inline T powi4(const T & x) 
  {
    T temp = x*x; // x^2
    return temp*temp;
  }
  template <typename T>
  inline T powi5(const T & x) 
  {
    T temp = x*x;
    return temp*temp*x;
  }
  template <typename T>
  inline T powi6(const T & x) 
  {
    T temp = x*x*x;
    return temp*temp;
  }
  template <typename T>
  inline T powi7(const T & x) 
  {
    T temp = x*x*x;
    return temp*temp*x;
  }
  template <typename T>
  inline T powi8(const T & x) 
  {
    T temp = x*x; // x^2
    temp *= temp; // x^4
    return temp*temp;
  }
  template <typename T>
  inline T powi9(const T & x) 
  {
    T temp = x*x*x;
    return temp*temp*temp;
  }
  template <typename T>
  inline T powi(T x, unsigned y) 
  {
    T res(1);
    while (y) {
      if (y & 1)
        res *= x;
      x *= x;
      y >>= 1;
    }
    return res;
  }

  // Minkosvki product
  template <typename T, typename U>
  inline ninja::Complex ninjaMP(const T & p, const U & u)
  {
    return ninja::mp(p,u);
  }

} // namespace
