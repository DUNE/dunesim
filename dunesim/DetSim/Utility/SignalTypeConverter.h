// SignalTypeConverter.h
//
// David Adams
// November 2015
//
// Utility class to convert between int and floating point types.
// Float to int conversion is rounding to nearest int:
//   (-1.5 , 0.5] --> -1
//   (-0.5 , 0.5) --> 0
//   [ 0.5 , 1.5) --> 1
//
// Future development may add parameters which specify different
// conversion strategies.
//
// Examples of usage.
//   float val = 1.2;
//   SignalTypeConverter sigcon.
//   int ival1 = sigcon.floatToInt(val);
//   short ival2 = sigcon.convert<short>(val);

#include <type_traits>

class SignalTypeConverter {

public:

  // Ctor.
  SignalTypeConverter();

  // Convert type T1 to type T2.
  // Uses floatToInt if If T1 is floating and T2 is integral.
  // Otherwise, standard type conversion.
  // Usage: int ival = convert<int>(myfloat);
  template<typename T2, typename T1>
  T2 convert(T1 val) const {
    if ( std::is_floating_point<T1>() && std::is_integral<T2>() ) {
      return floatToInt<T2, T1>(val);
    } else {
      return T2(val);
    }
  }

  // Conversion from floating point to int with arbitrary types.
  template<typename INT =int, typename FLOAT>
  INT floatToInt(FLOAT val) const {
    if ( val < 0.0 ) return -floatToInt<INT, FLOAT>(-val);
    return INT(val + 0.5);
  }

};
