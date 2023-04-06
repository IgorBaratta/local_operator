#include <vector>
#include <array>
#include <iostream>

template <typename T, typename S>
std::vector<T> create_geometry(int num_batches, int batch_size, int geom_size)
{
    std::vector<S>
        coords = {0., 0., 0., 0., 1., 1., 1., 1.,  // X coordinates
                  0., 0., 1., 1., 0., 0., 1., 1.,  // Y coordinates
                  0., 5., 0., 5., 0., 5., 0., 5.}; // Z coordinates

    std::vector<T> geometry(num_batches * geom_size);
    if constexpr (std::is_same<T, S>::value)
    {
        for (int c = 0; c < num_batches; c++)
            for (int i = 0; i < geom_size; i++)
                geometry[c * geom_size + i] = coords[i];
    }
    else
    {
        for (int c = 0; c < num_batches; c++)
            for (int i = 0; i < geom_size; i++)
                for (int j = 0; j < batch_size; j++)
                    geometry[c * geom_size + i][j] = static_cast<S>(coords[i]);
    }
    return geometry;
}
