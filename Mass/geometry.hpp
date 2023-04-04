#include <vector>
#include <array>
#include <iostream>

template <typename T, typename S>
std::vector<T> create_geometry(int num_batches, int batch_size, int geom_size)
{

    // Geometry of the cell
    //  [[0. 0. 0.]
    //  [1. 0. 0.]
    //  [1. 1. 0.]
    //  [1. 1. 1.]]
    // Coords stores in column major format XXXXYYYYZZZZ
    std::array<double, 12> coords = {0., 1., 1., 1.,
                                     0., 0., 1., 1.,
                                     0., 0., 0., 1.};

    std::vector<T> geometry(num_batches * geom_size);
    if constexpr (std::is_same<T, double>::value or std::is_same<T, float>::value)
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