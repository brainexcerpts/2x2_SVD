
#include "svd_2x2.hpp"

// =============================================================================
namespace Tbx {
// =============================================================================

// =============================================================================
namespace Details {
// =============================================================================

/**
  * @class JacobiRotation
  * @brief Rotation given by a cosine-sine pair.
  *
  * This class represents a Jacobi or Givens rotation.
  * This is a 2D rotation in the plane \c J of angle \f$ \theta \f$ defined by
  * its cosine \c c and sine \c s as follows:
  * \f$ J = \left ( \begin{array}{cc} c & \overline s \\ -s  & \overline c \end{array} \right ) \f$
  *
  * You can apply the respective counter-clockwise rotation to a column vector
  * \c v by* applying its adjoint on the left: \f$ v = J^* v \f$
  * \endcode
  */
template<typename Real>
class Jacobi_rotation_2x2 {
public:
    /// Default constructor without any initialization.
    Jacobi_rotation_2x2() { }

    /// Construct a planar rotation from a cosine-sine pair (\a c, \c s).
    Jacobi_rotation_2x2(const Real& c, const Real& s) : _cos(c), _sin(s) {}

    Real& cos() { return _cos; }
    Real& sin() { return _sin; }
    Real cos() const { return _cos; }
    Real sin() const { return _sin; }

    /// Concatenates two planar rotation
    Jacobi_rotation_2x2 operator*(const Jacobi_rotation_2x2& jr)
    {
        return Jacobi_rotation_2x2(_cos * jr._cos - _sin * jr._sin,
                                   _cos * jr._sin + _sin * jr._cos);
    }

    /// Returns the transposed transformation
    Jacobi_rotation_2x2 transpose() const { return Jacobi_rotation_2x2(_cos, -_sin); }

    /// Makes \c *this as a Jacobi rotation \c J such that applying \a J on both
    /// the right and left sides of the 2x2 selfadjoint matrix
    /// \f$ B = \left ( \begin{array}{cc} \text{this}_{pp} & \text{this}_{pq} \\ (\text{this}_{pq})^* & \text{this}_{qq} \end{array} \right )\f$
    /// yields a diagonal matrix \f$ A = J^* B J \f$
    inline bool make_jacobi(Mat2& m, int p, int q)
    {
        return make_jacobi( m(p, p), m(p, q), m(q, q) );
    }

    /// Makes \c *this as a Jacobi rotation \a J such that applying \a J on both
    /// the right and left sides of the selfadjoint 2x2 matrix
    /// \f$ B = \left ( \begin{array}{cc} x & y \\ \overline y & z \end{array} \right )\f$
    /// yields a diagonal matrix \f$ A = J^* B J \f$
    bool make_jacobi(Real x, Real y, Real z);

    /// Applies the rotation in the plane \a j to the rows \a p and \a q
    /// of \c m, i.e., it computes B = J * B,
    /// with \f$ B = \left ( \begin{array}{cc} \text{m.row}(p) \\ \text{m.row}(q) \end{array} \right ) \f$.
    inline void apply_on_the_left(Mat2& m, int p, int q)
    {
        for(int i = 0; i < 2; ++i)
        {
            Real xi = m(p, i);
            Real yi = m(q, i);
            m(p, i) =  cos() * xi + sin() * yi;
            m(q, i) = -sin() * xi + cos() * yi;
        }

    }

    /// Applies the rotation in the plane \a j to the columns \a p and \a q of
    /// \c m, i.e., it computes B = B * J with
    /// \f$ B = \left ( \begin{array}{cc} \text{m.col}(p) & \text{m.col}(q) \end{array} \right ) \f$.
    inline void apply_on_the_right(Mat2& m, int p, int q)
    {
        for(int i = 0; i < 2; ++i)
        {
            Real xi = m(i, p);
            Real yi = m(i, q);
            m(i, p) =  cos() * xi - sin() * yi;
            m(i, q) =  sin() * xi + cos() * yi;
        }
    }

private:
    Real _cos, _sin;
};

// -----------------------------------------------------------------------------

template<typename Real>
bool Jacobi_rotation_2x2<Real>::make_jacobi(Real x, Real y, Real z)
{
    if(y == Real(0))
    {
        _cos = Real(1);
        _sin = Real(0);
        return false;
    }
    else
    {
        Real tau = (x-z) / ( Real(2) * std::abs(y) );
        Real w = std::sqrt(tau*tau + Real(1));
        Real t;
        if(tau > Real(0))
        {
            t = Real(1) / (tau + w);
        }
        else
        {
            t = Real(1) / (tau - w);
        }
        Real sign_t = t > Real(0) ? Real(1) : Real(-1);
        Real n = Real(1) / std::sqrt( t*t + Real(1) );
        _sin = -sign_t * ( y / std::abs( y ) ) * std::abs( t ) * n;
        _cos = n;
        return true;
    }
}

// -----------------------------------------------------------------------------

template< typename Real >
void real_2x2_jacobi_svd(
        const Mat2& matrix,
        Jacobi_rotation_2x2<Real>* j_left,
        Jacobi_rotation_2x2<Real>* j_right)
{
    Mat2 m;

    m(0, 0) = matrix(1, 1); m(0, 1) = matrix(1, 0);
    m(1, 0) = matrix(0, 1); m(1, 1) = matrix(0, 0);

    Jacobi_rotation_2x2<Real> rot;

    Real t = m(0, 0) + m(1, 1);
    Real d = m(1, 0) - m(0, 1);

    if(t == Real(0))
    {
        rot.cos() = Real(0);
        rot.sin() = d > Real(0) ? Real(1) : Real(-1);
    }
    else
    {
        Real u = d / t;
        rot.cos() = Real(1) / std::sqrt(Real(1) + (u*u));
        rot.sin() = rot.cos() * u;
    }

    rot.apply_on_the_left(m, 0, 1);

    j_right->make_jacobi(m, 0, 1);

    *j_left = rot * j_right->transpose();
}

} // END Details namespace -----------------------------------------------------

template<bool U, bool V>
void SVD_2x2<U, V>::compute(const Mat2& matrix)
{
    // currently we stop when we reach precision 2*epsilon as the last bit of
    // precision can require an unreasonable number of iterations,
    // only worsening the precision of U and V as we accumulate more rotations
    const Real precision = Real(2) * std::numeric_limits<Real>::epsilon();

    // limit for very small denormal numbers to be considered zero in order to
    // avoid infinite loops
    const Real consider_null = Real(2) * std::numeric_limits<Real>::denorm_min();

    _work_mat = matrix;

    if( compute_u() ) _matrix_u = Mat2::identity();
    if( compute_v() ) _matrix_v = Mat2::identity();

    // The main Jacobi SVD iteration
    // if this 2x2 matrix is not diagonal already...
    // notice that this comparison will evaluate to false if any NaN is involved,
    // similarly, small denormal numbers are considered zero.
    Real max_diag = std::max( std::abs( _work_mat(1, 1) ),
                              std::abs( _work_mat(0, 0) ) );

    Real threshold = std::max(consider_null, precision * max_diag);

    Real max_anti_diag = std::max( std::abs( _work_mat(1, 0) ),
                                   std::abs( _work_mat(0, 1) ) );

    if( max_anti_diag > threshold)
    {
        Details::Jacobi_rotation_2x2<Real> j_left, j_right;
        Details::real_2x2_jacobi_svd(_work_mat, &j_left, &j_right);

        // accumulate resulting Jacobi rotations
        j_left.apply_on_the_left(_work_mat, 1, 0);
        if( compute_u() ) (j_left.transpose()).apply_on_the_right(_matrix_u, 1, 0);

        j_right.apply_on_the_right(_work_mat, 1, 0);
        if( compute_v() ) j_right.apply_on_the_right(_matrix_v, 1, 0);
    }

    // The work matrix is now diagonal,
    // so ensure it's positive so its diagonal entries are the singular values
    for(int i = 0; i < 2; ++i)
    {
        Real a = std::abs( _work_mat(i, i) );
        _singular_values[i] = a;
        if( compute_u() && ( a != Real(0) ) )
            _matrix_u.set_col(i, _matrix_u.col(i) *  _work_mat(i, i) / a);
    }

    // Sort singular values in descending order
    // and compute the number of nonzero singular values
    _non_zero_singular_val = 2;

    if( _singular_values.x < _singular_values.y )
    {
        if(_singular_values.x == Real(0))
        {
            // y strictly greater so not null
            _non_zero_singular_val = 1;
            return;
        }
        std::swap(_singular_values.x, _singular_values.y);
        if( compute_u() ) _matrix_u.swap_cols( 0 );
        if( compute_v() ) _matrix_v.swap_cols( 0 );
    }
    else
    {
        if(_singular_values.x == Real(0))
            _non_zero_singular_val = 0; // x greater or equal so y is also null

    }

    return;
}

}// END Tbx NAMESPACE ==========================================================
