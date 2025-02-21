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

#include "maxhs/utils/options.h"
#include <algorithm>
#include <cctype>
#include <iostream>
#include <limits>
#include <string>

using namespace MaxHS;
using std::cout;
using std::string;
using std::string_view;
using std::vector;

std::string Options::usage{};

bool Options::parse_options(int& argc, char** argv) {
  /* -help, -h, or -help-verb print help. Process all command line
     options that are known. These are removed from argv. Other
     command line argument (positional or unknown options)
     remain. False is returned if -h, -help, -help-verb or an invalid
     option (an arg starting with '-' is encountered). Stop parsing if
     there is a '--' arg. The '--' arg is removed and the remain args
     kept.
   */
  bool err{false};
  int j{1};
  for (int i{1}; i < argc; ++i) {
    const string_view arg{argv[i]};
    DEBUG_PR<options_DEBUG>("parsing i = ", i, " j=", j, " arg = ", arg, '\n');
    if (arg == "-help-verb" || arg == "--help-verb") {
      print_usage(true);
      err = true;
    } else if (arg == "-help" || arg == "-h" || arg == "--help") {
      print_usage(false);
      err = true;
    } else if (arg == "--") {
      for (++i; i < argc; ++i) argv[j++] = argv[i];
      break;
    } else if (arg.front() == '-') {
      DEBUG_PR<options_DEBUG>("checking option\n");
      bool parsed = std::any_of(
          Options::get_options().begin(), Options::get_options().end(),
          [&arg](Options* opt) {
            DEBUG_PR<options_DEBUG>("checking ", arg, " against ", opt->name, '\n');
            return opt->parse(arg);
          });
      if (!parsed) {  // unknown option (keep)
        cout << "Error: can't parse option " << arg << " (ignoring).\n";
        argv[j++] = argv[i];
        err = true;
      }
    } else {
      // non-option arg keep
      argv[j++] = argv[i];
    }
  }
  argc = j;
  return !err;
}

void Options::print_usage(bool verbose) {
  std::cout << Options::get_usage();
  if (!Options::get_usage().empty()) std::cout << '\n';
  sort(Options::get_options().begin(), Options::get_options().end(), opt_lt);

  string_view prev_cat{};
  for (const auto opt : Options::get_options()) {
    string_view cat{opt->category};
    if (cat != prev_cat) cout << '\n' << cat << " OPTIONS:\n" << '\n';
    opt->help(verbose);
    prev_cat = cat;
  }

  cout << "\nOTHER OPTIONS:\n\n";
  cout << "  --             pass subsequent arguments to program without "
          "parsing\n";
  cout << "  -h | -help    Print help message.\n";
  cout << "  -help-verb    Print verbose help message.\n\n";
}

const string Options::format_description(const string& d) {
  constexpr int line_size{61};
  const string prefix{"         "s};
  if (d.size() <= line_size)
    return ("         " + d);
  else {
    string fmt{};
    int cur_line_size{0};
    for (auto c : d) {
      if (cur_line_size == 0) {
        if (isspace(c))  // skip space at start of new line
          continue;
        else
          fmt += prefix;  // reached non-space at start: add prefix
      }
      if (cur_line_size < line_size || !isspace(c)) {
        // don't break line at non-space
        fmt += c;
        ++cur_line_size;
      } else {
        // line long enough and reached space in description
        // note this space char is skipped
        fmt += '\n';
        cur_line_size = 0;
      }
    }
    return fmt;
  }
}

void Options::print_option_settings(std::ostream& out,
                                    const string_view& prefix) {
  constexpr int Opts_per_line{1};
  if (Options::get_options().empty()) return;
  sort(Options::get_options().begin(), Options::get_options().end(), opt_lt);
  out << prefix << "PARAMETER SETTINGS\n";
  string_view cur_cat = Options::get_options().front()->category;
  out << prefix << '\n' << prefix << cur_cat << " options:\n";
  int on_line{0};
  for (auto opt : Options::get_options()) {
    if (cur_cat != opt->category) {
      cur_cat = opt->category;
      if (on_line != 0) out << '\n' << prefix << '\n';
      out << prefix << '\n' << prefix << cur_cat << " options:\n";
      on_line = 0;
    }
    if (on_line == 0) out << prefix;
    opt->print_setting(out);
    ++on_line;
    if (on_line != Opts_per_line)
      out << ",  ";
    else {
      out << '\n';
      on_line = 0;
    }
  }
}
