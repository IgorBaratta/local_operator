#include <cstddef>           // for std::size_t
#include <type_traits>       // for std::conditional_t

template <std::size_t Precision, std::size_t BatchSize>
struct VectorExtensions
{
    using S = std::conditional_t<Precision == 8, double, float>;
#if defined(__clang__)
    typedef S scalarB __attribute__((ext_vector_type(BatchSize)));
#elif defined(__GNUC__) || defined(__GNUG__)
    typedef S scalarB __attribute__((vector_size(BatchSize * sizeof(S))));
#else
#error "Compiler not supported"
#endif
    using T = std::conditional_t<BatchSize == 1, S, scalarB>;
};
