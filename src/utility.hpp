#ifndef TREESHREW_UTILITY_HPP
#define TREESHREW_UTILITY_HPP

#include <iostream>
#include <cmath>

namespace treeshrew {

//////////////////////////////////////////////////////////////////////////////
// Asserts

#define TREESHREW_ASSERT_RAISES_EXCEPTION 1

#if defined(IGNORE_TREESHREW_ASSERT) || defined(NDEBUG)
#   define TREESHREW_ASSERT(expr)
#   define TREESHREW_ASSERT_APPROX_EQUAL(x, y)
#else
    void treeshrew_assertion_failed(char const * expr, char const * function, char const * file, long line);
#   define TREESHREW_ASSERT(expr)  if (!(expr)) treeshrew::treeshrew_assertion_failed((const char *)#expr, (const char *)__FUNCTION__, __FILE__, __LINE__)
#   define TREESHREW_ASSERT_APPROX_EQUAL(x, y)  if (fabs(((x)-(y))/(x)) > 1.0e-6) std::cerr << std::fixed << std::setprecision(20) << (x) << ' ' << std::fixed << std::setprecision(20) << (y) << '\n'; if (fabs(((x)-(y))/(x)) > 1.0e-6) treeshrew::treeshrew_assert_approx_eq_failed((const char *)#x, x, (const char *)#y, y, (const char *)__FUNCTION__, __FILE__, __LINE__)
#endif

#define TREESHREW_NDEBUG_ASSERT(expr)  if (!(expr)) treeshrew::treeshrew_assertion_failed((const char *)#expr, (const char *)__FUNCTION__, __FILE__, __LINE__)
#define TREESHREW_NDEBUG_ASSERT_APPROX_EQUAL(x, y)  if (fabs(((x)-(y))/(x)) > 1.0e-6) std::cerr << std::fixed << std::setprecision(20) << (x) << ' ' << std::fixed << std::setprecision(20) << (y) << '\n'; if (fabs(((x)-(y))/(x)) > 1.0e-6) treeshrew::treeshrew_assert_approx_eq_failed((const char *)#x, x, (const char *)#y, y, (const char *)__FUNCTION__, __FILE__, __LINE__)

#define TREESHREW_DEBUG_GENE_OUTPUT(v)
#define TREESHREW_DEBUG_GENE_DUMP_VAR(v)
#define TREESHREW_DEBUG_JUMP_OUTPUT(v)
#define TREESHREW_DEBUG_JUMP_DUMP_VAR(v)
#if defined(TREESHREW_DEBUG_PRINTING) && TREESHREW_DEBUG_PRINTING
#   include <iostream>
#   define TREESHREW_DEBUG_OUTPUT(v) (std::cerr << v)
#   define TREESHREW_DEBUG_DUMP_VAR(x) (std::cerr << #x" = " << x << std::endl);
#   if defined(TREESHREW_DEBUG_GENE) && TREESHREW_DEBUG_GENE
#       undef TREESHREW_DEBUG_GENE_OUTPUT(v)
#       define TREESHREW_DEBUG_GENE_OUTPUT(v) (std::cerr << v)
#       undef TREESHREW_DEBUG_GENE_DUMP_VAR(v)
#       define TREESHREW_DEBUG_GENE_DUMP_VAR(x) (std::cerr << #x" = " << x << " (debug_gene_output)" << std::endl);
#   endif
#   if defined(TREESHREW_DEBUG_JUMP) && TREESHREW_DEBUG_JUMP
#       undef TREESHREW_DEBUG_JUMP_OUTPUT(v)
#       define TREESHREW_DEBUG_JUMP_OUTPUT(v) (std::cerr << v)
#       undef TREESHREW_DEBUG_JUMP_DUMP_VAR(v)
#       define TREESHREW_DEBUG_JUMP_DUMP_VAR(x) (std::cerr << #x" = " << x << " (debug_jump_output)" << std::endl);
#   endif
#else
#   define TREESHREW_DEBUG_OUTPUT(v)
#   define TREESHREW_DEBUG_DUMP_VAR(x)
#endif

void treeshrew_assertion_failed(char const * expr, char const * function, char const * file, long line);
void treeshrew_assert_approx_eq_failed(char const * x, double val_x, const char * y, double val_y, char const * function, char const * file, long line);

//////////////////////////////////////////////////////////////////////////////
// Printing

template <typename S>
void treeshrew_print(S&) {}

template <typename S, typename T, typename... Types>
void treeshrew_print(S& stream, const T& arg1, const Types&... args) {
        stream << arg1;          // print first argument
        treeshrew_print(stream, args...);   // call treeshrew_stderr() for remaining arguments }
}

template <typename... Types>
void treeshrew_abort(const Types&... args) {
    treeshrew_print(std::cerr, "-treeshrew: ", args...);
    std::cerr << std::endl;
    exit(1);
}

} // treeshrew

#endif
