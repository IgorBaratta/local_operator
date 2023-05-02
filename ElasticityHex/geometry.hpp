#include <vector>
#include <array>

template <typename T, typename S>
std::vector<T> create_geometry(int num_batches, int cub_nq, bool precompute)
{
    constexpr int size_quad = 9;
    std::vector<T> geometry(num_batches * cub_nq * size_quad);
    std::array<S, size_quad> G = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    for (int b = 0; b < num_batches; b++)
        for (int i = 0; i < cub_nq; i++)
            for (std::size_t ig = 0; ig < G.size(); ig++)
                geometry[b * cub_nq * size_quad + ig * cub_nq + i] = static_cast<T>(G[ig]);
    return geometry;
}
