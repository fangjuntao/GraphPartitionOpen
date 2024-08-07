#ifndef KMEANS_SIMPLE_MATRIX_HPP
#define KMEANS_SIMPLE_MATRIX_HPP
 
namespace kmeans {
 
template<typename Data_, typename Index_, typename Dim_ = int>
class SimpleMatrix {
public:
    SimpleMatrix(Dim_ num_dimensions, Index_ num_observations, const Data_* data) : 
        my_num_dim(num_dimensions), my_num_obs(num_observations), my_data(data), my_long_num_dim(num_dimensions) {}
 
private:
    Dim_ my_num_dim;
    Index_ my_num_obs;
    const Data_* my_data;
    size_t my_long_num_dim;
 
public:
    typedef Data_ data_type;
 
    typedef Index_ index_type;
 
    typedef Dim_ dimension_type;
 
    struct RandomAccessWorkspace{};
 
    struct ConsecutiveAccessWorkspace {
        ConsecutiveAccessWorkspace(index_type start) : at(start) {}
        size_t at;
    };
 
    struct IndexedAccessWorkspace {
        IndexedAccessWorkspace(const index_type* sequence) : sequence(sequence) {}
        const index_type* sequence;
        size_t at = 0;
    };
 
public:
    Index_ num_observations() const {
        return my_num_obs;
    }
 
    dimension_type num_dimensions() const {
        return my_num_dim;
    }
 
    RandomAccessWorkspace create_workspace() const {
        return RandomAccessWorkspace();
    }
 
    ConsecutiveAccessWorkspace create_workspace(index_type start, index_type) const {
        return ConsecutiveAccessWorkspace(start);
    }
 
    IndexedAccessWorkspace create_workspace(const index_type* sequence, index_type) const {
        return IndexedAccessWorkspace(sequence);
    }
 
    const data_type* get_observation(Index_ i, [[maybe_unused]] RandomAccessWorkspace& workspace) const {
        return my_data + static_cast<size_t>(i) * my_long_num_dim; // avoid overflow during multiplication.
    } 
 
    const data_type* get_observation(ConsecutiveAccessWorkspace& workspace) const {
        return my_data + (workspace.at++) * my_long_num_dim; // everything is already a size_t.
    } 
 
    const data_type* get_observation(IndexedAccessWorkspace& workspace) const {
        return my_data + static_cast<size_t>(workspace.sequence[workspace.at++]) * my_long_num_dim; // avoid overflow during multiplication.
    } 
};
 
}
 
#endif
