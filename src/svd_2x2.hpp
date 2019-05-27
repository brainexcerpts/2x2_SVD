#ifndef TOOL_BOX_SVD_2x2_HPP__
#define TOOL_BOX_SVD_2x2_HPP__

#include <cassert>
#include "mat2.hpp"
#include "vec2.hpp"

// =============================================================================
namespace Tbx {
// =============================================================================

/**
 * @class SVD_2x2
 * @brief compute a singular value decomposition (SVD) of a 2x2 matrix
 *
 * The SVD decompose the matrix M = U * S * V^t.
 * With S the diagonal of the singular values of M
 *
 * @tparam comput_U : wether we compute or not the matrix U of the the SVD
 * @tparam comput_V : wether we compute or not the matrix V of the the SVD
 *
 * @note Code adapted from Eigen 3 (http://eigen.tuxfamily.org) from the general
 * Jacobi SVD of a NxM matrix.
 */
template< bool compute_U, bool compute_V >
class SVD_2x2 {
public:
    typedef float Real;

    SVD_2x2(const Mat2& matrix) { compute( matrix ); }

    void compute(const Mat2& matrix);

    const Mat2& matrix_u() const
    {
        assert( compute_u() );
        return _matrix_u;
    }

    const Mat2& matrix_v() const
    {
        assert( compute_v() );
        return _matrix_v;
    }

    ///  Singular values are always sorted in decreasing order.
    Vec2 singular_values() const { return _singular_values; }

    /// @return if user asked for computation of U
    inline bool compute_u() const { return compute_U; }
    /// @return if user asked for computation of V
    inline bool compute_v() const { return compute_U; }

    /// @return the number of singular values that are not exactly 0
    int non_zero_singular_values() const { return _non_zero_singular_val; }

private:
    Mat2 _matrix_u;
    Mat2 _matrix_v;
    Vec2 _singular_values;

    Mat2 _work_mat;
    int _non_zero_singular_val;
};

}// END Tbx NAMESPACE ==========================================================

#include "svd_2x2.inl"

#endif // TOOL_BOX_SVD_2x2_HPP__
