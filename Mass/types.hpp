#include <cstddef>     // for std::size_t
#include <type_traits> // for std::conditional_t

template <std::size_t Precision, std::size_t BatchSize>
struct VectorExtensions
{
    #if PRECISION == 8
        typedef double T __attribute__((vector_size(BatchSize * sizeof(double))));
    #else
        typedef float T __attribute__((vector_size(BatchSize * sizeof(float))));
    #endif
    using S = std::conditional_t<Precision == 8, double, float>;
};

template <std::size_t Precision>
struct VectorExtensions<Precision, 1>
{
    using T = std::conditional_t<Precision == 8, double, float>;
    using S = std::conditional_t<Precision == 8, double, float>;
};
