#ifndef KMEANS_MOCK_MATRIX_HPP
#define KMEANS_MOCK_MATRIX_HPP
 
namespace kmeans {
 
class MockMatrix {
public:
    MockMatrix(int num_dim, int num_obs, const double* data) : my_num_dim(num_dim), my_num_obs(num_obs), my_data(data), my_long_num_dim(num_dim) {}
public:
    typedef double data_type;
 
    typedef int index_type;
 
    typedef int dimension_type;
 
private:
    dimension_type my_num_dim;
    index_type my_num_obs;
    const data_type* my_data;
    size_t my_long_num_dim;
 
public:
    index_type num_observations() const {
        return my_num_obs;
    }
 
    dimension_type num_dimensions() const {
        return my_num_dim;
    }
 
public:
    struct RandomAccessWorkspace {};
 
    RandomAccessWorkspace create_workspace() const {
        return RandomAccessWorkspace();
    }
 
    struct ConsecutiveAccessWorkspace {
        ConsecutiveAccessWorkspace(index_type start) : at(start) {}
        size_t at;
    };
 
    ConsecutiveAccessWorkspace create_workspace(index_type start, [[maybe_unused]] index_type length) const {
        return ConsecutiveAccessWorkspace(start);
    }
 
    struct IndexedAccessWorkspace {
        IndexedAccessWorkspace(const index_type* sequence) : sequence(sequence) {}
        const index_type* sequence;
        size_t at = 0;
    };
 
    IndexedAccessWorkspace create_workspace(const index_type* sequence, [[maybe_unused]] index_type length) const {
        return IndexedAccessWorkspace(sequence);
    }
 
public:
    const data_type* get_observation(int i, [[maybe_unused]] RandomAccessWorkspace& workspace) const {
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