#ifndef MAXHS_TYPENAMES_H
#define MAXHS_TYPENAMES_H

/*********************************************************************
Copyright (c) 2021, Fahiem Bacchus

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*********************************************************************/
#include <typeinfo>
#include <string>
using namespace std::string_literals;

template <typename T> struct TypeName {
  // Generally only some compilers will give a sensible name for
  // an arbitary type. So this struct should be specialised with the types
  // you need names for in your program
  static const std::string get() { return typeid(T).name(); }
};

// a specialization for each type of those you want to support
// and don't like the string returned by typeid
template <> struct TypeName<int32_t> {
  static const std::string get() { return "int32_t"s; }
};

template <> struct TypeName<uint32_t> {
  static const std::string get() { return "uint32_t"s; }
};

template <> struct TypeName<int64_t> {
  static const std::string get() { return "int64_t"s; }
};

template <> struct TypeName<uint64_t> {
  static const std::string get() { return "uint64_t"s; }
};

template <> struct TypeName<long double> {
  static const std::string get() { return "long double"s; }
};

template <> struct TypeName<double> {
  static const std::string get() { return "double"s; }
};

template <> struct TypeName<float> {
  static const std::string get() { return "float"s; }
};

template <> struct TypeName<bool> {
  static const std::string get() { return "bool"s; }
};

template <> struct TypeName<char> {
  static const std::string get() { return "char"s; }
};

template <> struct TypeName<std::string> {
  static const std::string get() { return "string"s; }
};

#endif
