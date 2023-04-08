#include <vector>
#include <array>
#include <iostream>

template <typename T, typename S>
std::vector<T> create_geometry(int num_batches, int batch_size, int cub_nq, bool precompute)
{
    if (precompute)
    {
        std::vector<T> geometry(num_batches * cub_nq * 6);
        std::array<S, 6> G = {1, 0, 0, 1, 0, 1};
        for (int b = 0; b < num_batches; b++)
            for (int i = 0; i < cub_nq; i++)
                for (std::size_t ig = 0; ig < G.size(); ig++)
                    geometry[b * cub_nq * 6 + ig * cub_nq + i] = static_cast<T>(G[ig]);

        return geometry;
    }
    else
    {
        std::vector<S>
            coords = {0., 0., 0., 0., 1., 1., 1., 1.,  // X coordinates
                      0., 0., 1., 1., 0., 0., 1., 1.,  // Y coordinates
                      0., 5., 0., 5., 0., 5., 0., 5.}; // Z coordinates

        std::vector<T> geometry(num_batches * coords.size());
        if constexpr (std::is_same<T, S>::value)
        {
            for (int c = 0; c < num_batches; c++)
                for (std::size_t i = 0; i < coords.size(); i++)
                    geometry[c * coords.size() + i] = coords[i];
        }
        else
        {
            for (int c = 0; c < num_batches; c++)
                for (std::size_t i = 0; i < coords.size(); i++)
                    for (int j = 0; j < batch_size; j++)
                        geometry[c * coords.size() + i][j] = static_cast<S>(coords[i]);
        }
        return geometry;
    }
}
