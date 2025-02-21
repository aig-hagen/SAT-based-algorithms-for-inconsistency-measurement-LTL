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
#ifndef MAXHS_OPTIONS_H
#define MAXHS_OPTIONS_H

#include <charconv>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include "maxhs/utils/io.h"
#include "maxhs/utils/typenames.h"

namespace MaxHS {

//================================================================================
// Options (pural) is an abstract class that gives the interface for
// all types options. Specific type options are Option<T> singular.

// USAGE:
/*
  - call Options::set_usage_help with an already formatted string. This string
  will be printed if the command line contains -h, -help, or --help-verb

  - set up individual Option<T> objects. If more than one file needs
    access to the option, extern the definition and declare it in the
    other files that need it...or declare shared parameters in
    maxhs/utils/params.h and include that in all files.

  - call Options::print_option_settings if you want to print out the state of
    all options.

  - call Options::parse_options(argc, argv) to process the command line (CL).
      1) if CL contains -h, -help, -help-verb usage help will be
         printed to cout.

      2) every argv item will be processed. If it is a valid option it
         will be removed from argv and the corresponding Option<T>
         will be set. argc and argv on return will contain all
         unparsable options. Note options start with '-' so arguments
         without initial '-' will be preserved in the order given.
         E.g., if the program expects [Option<T> ... Option<T>]
         <filename> and CL contains only valid options argc on return
         will be 1 and argv will contain only <filename>

         You can check the remaining elements of argv to determine if
         an error occured.

         Additionally, if CL contains the argv item "--" the
         subsequent items of argv are preserved and the "--" is
         removed. "--" is used to indicate the end of options on the
         CL. Note that if invalid options before "--" are encountered
         they also will be in argv on return and will be located
         before "--"

    3) false will be returned if CL contains any invalid options
       (i.e., argv items starting with '-') or -h, -help, -help-verb
 */

constexpr bool options_DEBUG{false};

class Options {
 public:
  // set string for usage help
  static void set_usage_help(const std::string& str) { Options::usage = str; }
  static bool parse_options(int& argc, char** argv);
  static void print_option_settings(std::ostream& out,
                                    const std::string_view& prefix = "c ");
  static constexpr int No_Limit{-1};
  static constexpr int Unsigned_No_Limit{0};

 private:
  static std::string usage;
  static void print_usage(bool verbose = false);
  static std::vector<Options*>& get_options() {
    static std::vector<Options*> options;
    return options;
  }
  static std::string& get_usage() { return usage; }

 protected:
  // specific Option<T> (notice Option singular not Options pural)
  virtual void print_setting(std::ostream& out) = 0;
  static const std::string format_description(const std::string& desc);
  virtual bool parse(std::string_view cl_str) = 0;
  virtual void help(bool verbose = false) {
    std::cout << std::left << std::setw(20) << "-" + name << "=<" << type_name
              << ">\n";
    if (verbose) std::cout << Options::format_description(description) << '\n';
  }

  const std::string name;
  const std::string description;
  const std::string category;
  const std::string type_name;

  bool parse_name_value(std::string_view cl_arg) {
    // check if arg matches cl arg of the form -<name>=<value>
    if (cl_arg[0] != '-' || cl_arg.substr(1, name.size()) != name) return false;
    // name matches
    if (cl_arg[name.size() + 1] != '=') {
      std::cout << "Error: option -" << name << " requires value (in format "
                << name << "=<value>)\n";
      return false;
    }
    return true;
  }

  Options(const std::string& c, const std::string& n, const std::string& d,
          const std::string& t)
      : name(n), description(d), category(c), type_name(t) {
    DEBUG_PR<options_DEBUG>("adding option ", name, " to get_options()\n");
    get_options().push_back(this);
    DEBUG_PR<options_DEBUG>("# options ", get_options().size(), '\n');
  }
  bool operator<(const Options& opt) const {
    auto cmp = category.compare(opt.category);
    return (cmp < 0 || (cmp == 0 && name < opt.name));
  }
  friend bool opt_lt(const Options*, const Options*);
  static void print_usage_and_exit(bool verbose);
};  // namespace MaxHS

inline bool opt_lt(const Options* x, const Options* y) { return *x < *y; }

//==================================================================================
// Specific Option types
template <typename T> class Option : public Options {
  static_assert(std::is_arithmetic_v<T>, "Cannot define Option of this type\n");

 protected:
  T value;
  T lb;
  T ub;

 public:
  Option(const std::string& cat, const std::string& nm, const std::string& desc,
         T val, T lower = std::numeric_limits<T>::min(),
         T upper = std::numeric_limits<T>::max())
      : Options(cat, nm, desc, TypeName<T>::get()),
        value(val),
        lb(lower),
        ub(upper) {}

  operator T(void) const { return value; }
  Option& operator=(T x) {
    value = x;
    return *this;
  }

  bool parse(std::string_view cl_arg) {
    if (!parse_name_value(cl_arg)) return false;
    T parsed_value;
    /*
    // Works with gcc-11 and clang but from_chars for floats
    // not yet implemented in gcc < 11.
    auto [ptr, ec] =
      std::from_chars((const char*) cl_arg.data() + name.size() + 2,
                      (const char*) cl_arg.data() + cl_arg.size(),
    parsed_value); if (ec == std::errc::invalid_argument) { std::cout << "Error:
    (" << cl_arg << ") option " << name
                << " supplied with invalid value (should be of type "
                << type_name << ")\n";
      return false;
    }
    if (ec == std::errc::result_out_of_range) {
      std::cout << "Error: (" << cl_arg << ") option " << name
                << " supplied with value that cannot fit into type "
                << type_name << '\n';
      return false;
    }
    if (ec != std::errc()) {
      std::cout << "Error: (" << cl_arg << ") while parsing value for " << name
                << '\n';
      return false;
    }
    */
    std::istringstream ss(std::string(cl_arg.substr(name.size() + 2)));
    ss >> parsed_value;
    if (!ss) {
      std::cout << "Error: (" << cl_arg << ") option " << name
                << " supplied with invalied type (should be of type"
                << type_name << ")\n";
    }
    if (parsed_value > ub) {
      std::cout << "Error: (" << cl_arg << ") option " << name
                << " given too large a value (max = " << ub << ")\n";
      return false;
    }
    if (parsed_value < lb) {
      std::cout << "Error: (" << cl_arg << ") option " << name
                << " given too small a value (min = " << lb << ")\n";
      return false;
    }
    DEBUG_PR<options_DEBUG>("parsed numeric option ", name, " got ", name, "=",
                            parsed_value, '\n');
    value = parsed_value;
    return true;
  }

  void help(bool verbose = false) {
    auto prec = cout.precision();
    auto flags = cout.flags();
    
    std::cout << "  " << std::left << std::setw(25) << "-" + name << " = "
              << std::left << std::setw(8) << std::setprecision(1) << std::fixed << value;
    cout.precision(prec);
    cout.setf(flags);
    
    std::cout << "[";
    if (lb == std::numeric_limits<T>::min() && lb != 0)
      std::cout << "min";
    else
      std::cout << lb;
    std::cout << "--";
    if (ub == std::numeric_limits<T>::max())
      cout << "max";
    else
      std::cout << ub;
    std::cout << "] " << '<' << type_name << ">\n";
    if (verbose) std::cout << Options::format_description(description) << '\n';
  }
  void print_setting(std::ostream& out) {
    out << "-" << name << " = " << value;
  }
};

//==================================================================================================
// Std::String option:

template <> class Option<std::string> : public Options {
 protected:
  std::string value;

 public:
  Option(const std::string& cat, const std::string& nm, const std::string& desc,
         const std::string& val)
      : Options(cat, nm, desc, TypeName<std::string>::get()), value(val) {}

  operator const std::string &(void) const { return value; }
  Option& operator=(const std::string& x) {
    value = x;
    return *this;
  }
  Option& operator=(std::string&& x) {
    value = x;
    return *this;
  }
  bool parse(std::string_view cl_arg) {
    if (!parse_name_value(cl_arg)) return false;
    value = cl_arg.substr(name.size() + 2, cl_arg.size() - name.size() - 2);
    DEBUG_PR<options_DEBUG>("parsed string option ", name, " got ", name, "=",
                            value, '\n');
    return true;
  }
  void help(bool verbose = false) {
    std::cout << "  " << std::left << std::setw(25) << "-" + name << " = "
              << std::left << std::setw(8) << value;
    std::cout << " <" << type_name << ">\n";
    if (verbose) std::cout << Options::format_description(description) << '\n';
  }
  void print_setting(std::ostream& out) {
    out << "-" << name << " = " << value;
  }
};

//==================================================================================================
// Bool option:

template <> class Option<bool> : public Options {
  bool value;

 public:
  Option(const std::string& cat, const std::string& nm, const std::string& desc,
         bool val)
      : Options(cat, nm, desc, TypeName<bool>::get()), value(val) {}

  operator bool(void) const { return value; }
  Option& operator=(bool b) {
    value = b;
    return *this;
  }

  bool parse(std::string_view cl_arg) {
    if (cl_arg.size() < name.size() + 1) return false;
    if (cl_arg[0] != '-') return false;
    bool negated = false;
    if (cl_arg.substr(1, 3) == "no-") {
      cl_arg.remove_prefix(4);
      negated = true;
    } else
      cl_arg.remove_prefix(1);
    if (cl_arg != name) return false;
    value = !negated;
    DEBUG_PR<options_DEBUG>("parsed bool option ", name, " got ", name, "=",
                            value, '\n');
    return true;
  }

  void help(bool verbose = false) {
    std::cout << "  " << std::left << std::setw(25) << "-" + name << " = "
              << std::left << std::setw(8) << (value ? "on" : "off") << "[-"
              << name << ", -no-" << name;
    std::cout << "] "
              << "<" << type_name << ">\n";
    if (verbose) std::cout << Options::format_description(description) << '\n';
  }

  virtual void print_setting(std::ostream& out) {
    out << "-" << name << " = " << (value ? "true" : "false");
  }
};

template <> class Option<char> : public Options {
  char value;

 public:
  Option(const std::string& cat, const std::string& nm, const std::string& desc,
         char val)
      : Options(cat, nm, desc, TypeName<char>::get()), value(val) {}

  operator char(void) const { return value; }
  Option& operator=(char b) {
    value = b;
    return *this;
  }

  bool parse(std::string_view cl_arg) {
    if (!parse_name_value(cl_arg)) return false;
    if (cl_arg.size() != name.size() + 3) {
      std::cout
          << "Error: (" << cl_arg << ") option " << name
          << " supplied with invalid value (should be single character)\n";
      return false;
    }
    value = cl_arg[name.size() + 2];
    DEBUG_PR<options_DEBUG>("parsed char option ", name, " got ", name, "=",
                            value, '\n');
    return true;
  }

  void help(bool verbose = false) {
    std::cout << "  " << std::left << std::setw(25) << "-" + name << " = "
              << std::left << std::setw(8) << value;
    std::cout << " <" << type_name << ">\n";
    if (verbose) std::cout << Options::format_description(description) << '\n';
  }
  void print_setting(std::ostream& out) {
    out << "-" << name << " = " << value;
  }
};

//=================================================================================================
}  // namespace MaxHS

#endif
