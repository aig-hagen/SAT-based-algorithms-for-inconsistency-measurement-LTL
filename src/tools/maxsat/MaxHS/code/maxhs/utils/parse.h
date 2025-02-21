/***********[parse.h]
Copyright (c) 2022 Fahiem Bacchus

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

***********/
#ifndef PARSE_H
#define PARSE_H

#include <string>

template <typename T>
void parse_num(const std::string& ln, const char*& pos, T& num) {
  auto eol = ln.data() + ln.size();
  while (pos < eol && std::isspace(*pos)) ++pos;
  auto [end_of_num, ec] = std::from_chars(pos, eol, num);
  pos = end_of_num;
  if (ec != std::errc())
    throw std::runtime_error("Invalid line:\"" + ln + "\"\n");
}

#if (__cpp_lib_to_chars < 201611L)
inline void parse_num(const std::string& ln, const char*& pos, double& num) {
  char* end_of_num;
  num = std::strtod(pos, &end_of_num);
  if (pos == end_of_num) {
    throw std::runtime_error("Invalid line:\"" + ln + "\"\n");
  }
  if (errno == ERANGE)
    throw std::runtime_error("Number out of range on line:\n\"" + ln + "\"\n");
  pos = end_of_num;
}
inline void parse_num(const std::string& ln, const char*& pos, long double& num) {
  char* end_of_num;
  num = std::strtold(pos, &end_of_num);
  if (pos == end_of_num) {
    throw std::runtime_error("Invalid line:\"" + ln + "\"\n");
  }
  if (errno == ERANGE)
    throw std::runtime_error("Number out of range on line:\n\"" + ln + "\"\n");
  pos = end_of_num;
}


#endif

#endif
