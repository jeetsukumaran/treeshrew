#include <iomanip>
#include <iostream>
#include "utility.hpp"

namespace treeshrew {

void treeshrew_assertion_failed(char const * expr, char const * function, char const * file, long line) {
    treeshrew_print(std::cerr,
            "\n(treeshrew) Assertion failed:",
            "\n  expr: " , expr,
            "\n  func: " , function,
            "\n  file: " , file,
            "\n  line: " , line);
    std::exit(1);
}

void treeshrew_assert_approx_eq_failed(char const * x,
        double val_x,
        const char * y,
        double val_y,
        char const * function,
        char const * file,
        long line) {
    treeshrew_print(std::cerr,
            "\n(treeshrew) Approximately equal assertion failed:",
            "\n  " , x , " (" , std::fixed , std::setprecision(20) , val_x , ") approx equal to " , y , " (" , std::fixed , std::setprecision(20) , val_y , ")",
            "\n  func: " , function,
            "\n  file: " , file,
            "\n  line: " , line);
    std::exit(1);
}

} // namespace treeshrew

