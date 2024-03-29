#include <vector>
#include <array>
#include <iostream>

template <typename T>
std::vector<T> create_geometry(int num_batches, int batch_size, int geom_size)
{
    std::array<double, 24> coords = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0,
                                     1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0};
    std::vector<T> geometry(num_batches * geom_size);
#ifdef USE_VECTOR_EXTENSIONS
    for (int c = 0; c < num_batches; c++)
        for (int i = 0; i < geom_size; i++)
            for (int j = 0; j < batch_size; j++)
                geometry[c * geom_size + i][j] = coords[i];
#else
    for (int c = 0; c < num_batches; c++)
        for (int i = 0; i < geom_size; i++)
            geometry[c * geom_size + i] = coords[i];
#endif
    return geometry;
}