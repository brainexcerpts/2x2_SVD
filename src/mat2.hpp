#ifndef TOOL_BOX_MAT2_HPP__
#define TOOL_BOX_MAT2_HPP__

#include "vec2.hpp"
#include <cmath>
#include <stdlib.h>
#include <stdio.h>

// =============================================================================
namespace Tbx {
// =============================================================================

/**
 * @name Mat2
 * @brief Handling 2x2 matrix
 *
 * @see Transfo Vec3
 */
struct Mat2 {

    /// Linear matrix storage with <b> rows first (i.e row major) </b>
    /// Using Mat2 with OpenGL can be done by transposing first:
    /// @code
    ///     Transfo tr;
    ///     // OpenGL is column major !
    ///     glMultMatrixf( (GLfloat)(tr.transpose().m) );
    /// @endcode
    float m[4];

    // -------------------------------------------------------------------------
    /// @name constructors
    // -------------------------------------------------------------------------

     inline Mat2() {   }

    /// Row first filling of the matrix
    inline
    Mat2(float a_, float b_,
         float c_, float d_)
    {
        m[0] = a_; m[1] = b_;
        m[2] = c_; m[3] = d_;
    }

    inline
    Mat2(const Vec2& x, const Vec2& y)
    {
        m[0] = x.x; m[1] = y.x;
        m[2] = x.y; m[3] = y.y;
    }

    // -------------------------------------------------------------------------
    /// @name operators overload
    // -------------------------------------------------------------------------

    //----------------
    // Multiplications
    //----------------

    inline Vec2 operator*(const Vec2& v) const
    {
        const float x = v.x * m[0] + v.y * m[1];
        const float y = v.x * m[2] + v.y * m[3];
        return Vec2(x, y);
    }

    inline Mat2 operator*(const Mat2& m) const
    {
        return Mat2(this->m[0] * m[0] + this->m[1] * m[2],
                    this->m[0] * m[1] + this->m[1] * m[3],

                    this->m[2] * m[0] + this->m[3] * m[2],
                    this->m[2] * m[1] + this->m[3] * m[3] );
    }

    inline Mat2 operator*(float x) const
    {
        return Mat2(m[0] * x, m[1] * x,
                    m[2] * x, m[3] * x);
    }

    inline Mat2& operator*=(float x)
    {
        m[0] *= x; m[1] *= x;
        m[2] *= x; m[3] *= x;
        return *this;
    }

    inline friend Mat2 operator*(const float x_, const Mat2& mat)
    {
        return Mat2(x_ * mat[0], x_ * mat[1],
                    x_ * mat[2], x_ * mat[3]);
    }

    //----------
    // Divisions
    //----------

    inline Mat2 operator/(float x) const
    {
        return Mat2(m[0] / x, m[1] / x,
                    m[2] / x, m[3] / x);
    }

    inline Mat2& operator/=(float x)
    {
        m[0] /= x; m[1] /= x;
        m[2] /= x; m[3] /= x;
        return *this;
    }

    inline friend Mat2 operator/(const float x_, const Mat2& mat)
    {
        return Mat2(x_ / mat[0], x_ / mat[1],
                    x_ / mat[2], x_ / mat[3]);
    }

    //----------
    // Additions
    //----------


    inline Mat2 operator+(const Mat2& m) const
    {
        return Mat2(this->m[0] + m[0], this->m[1] + m[1],
                    this->m[2] + m[2], this->m[3] + m[3]);
    }

    inline Mat2 operator+(float x) const
    {
        return Mat2(m[0] + x, m[1] + x,
                    m[2] + x, m[3] + x);
    }

    inline friend Mat2 operator+(const float x_, const Mat2& mat)
    {
        return Mat2(x_ + mat[0], x_ + mat[1],
                    x_ + mat[2], x_ + mat[3]);
    }

    inline Mat2& operator+=(float x)
    {
        m[0] += x; m[1] += x;
        m[2] += x; m[3] += x;
        return *this;
    }

    //--------------
    // Substractions
    //--------------


    inline Mat2 operator-(const Mat2& m) const
    {
        return Mat2(this->m[0] - m[0], this->m[1] - m[1],
                    this->m[2] - m[2], this->m[3] - m[3]);
    }

    inline Mat2 operator-() const
    {
        return Mat2(-m[0], -m[1],
                    -m[2], -m[3]);
    }

    inline Mat2 operator-(float x) const
    {
        return Mat2(m[0] - x, m[1] - x,
                    m[2] - x, m[3] - x);
    }

    inline friend Mat2 operator-(const float x_, const Mat2& mat)
    {
        return Mat2(x_ - mat[0], x_ - mat[1],
                    x_ - mat[2], x_ - mat[3]);
    }

    inline Mat2& operator-=(float x)
    {
        m[0] -= x; m[1] -= x;
        m[2] -= x; m[3] -= x;
        return *this;
    }

    // -------------------------------------------------------------------------
    /// @name Accessors
    // -------------------------------------------------------------------------

    //----------------
    // Access elements
    //----------------

    inline const float& operator()(int row, int column) const
    {
        assert(row >= 0 && row < 2);
        assert(column >= 0 && column < 2);
        return m[column + row*2];
    }

    inline float& operator()(int row, int column)
    {
        assert(row >= 0 && row < 2);
        assert(column >= 0 && column < 2);
        return m[column + row*2];
    }

    inline float& operator[](int idx) {
        assert(idx >= 0 && idx < 4);
        return m[idx];
    }

    inline const float& operator[](int idx) const {
        assert(idx >= 0 && idx < 4);
        return m[idx];
    }

    //----------------
    // Column access
    //----------------

    inline void set_col(int i, const Vec2& c) {
        assert(i >= 0 && i < 2);
        m[i    ] = c.x;
        m[i + 2] = c.y;
    }

    inline Vec2 col(int i) const {
        assert(i >= 0 && i < 2);
        return Vec2( m[i], m[i + 2]);
    }

    /// Swap the ith column with the other
    inline void swap_cols(int i)
    {
        assert(i >= 0 && i < 2);
        // Next column:
        const int j = (~i) & 0x01; // equivalent to (i+1) % 2

        // Swap first row
        float tmp = m[i];
        m[i] = m[j];
        m[j] = tmp;

        // Swap second row
        tmp = m[i + 2];
        m[i + 2] = m[j + 2];
        m[j + 2] = tmp;
    }

     inline Vec2 x() const { return Vec2(m[0], m[2]); }
     inline Vec2 y() const { return Vec2(m[1], m[3]); }

    //------------
    // Row access
    //------------

    inline void set_row(int i, const Vec2& c) {
        assert(i >= 0 && i < 2);
        m[i*2    ] = c.x;
        m[i*2 + 1] = c.y;
    }

    inline Vec2 row(int i) const {
        assert(i >= 0 && i < 2);
        return Vec2( m[i*2], m[i*2 + 1]);
    }

    /// Swap the ith row with the other
    inline void swap_rows(int i)
    {
        assert(i >= 0 && i < 2);
        // Next column:
        const int j = ((~i) & 0x01) * 2; // equivalent to ((i+1) % 2) * 2
        i *= 2;

        // Swap first row
        float tmp = m[i];
        m[i] = m[j];
        m[j] = tmp;

        // Swap second row
        tmp = m[i + 1];
        m[i + 1] = m[j + 1];
        m[j + 1] = tmp;
    }

    // -------------------------------------------------------------------------
    /// @name operations
    // -------------------------------------------------------------------------

    inline float det() const { return m[0] * m[3] - m[1] * m[2]; }

    /// @return the matrix with normalized x, y column vectors

    inline Mat2 normalized() const {
        return Mat2(x().normalized(), y().normalized());
    }

    inline Mat2 inverse() const
    {
        float idet = 1.f / det();
        return Mat2(  m[3], -m[1],
                     -m[2],  m[0]) * idet;
    }


    inline Mat2 transpose() const
    {
        return Mat2(m[0], m[2],
                    m[1], m[3]);
    }


    inline void set_abs()
    {
        m[0] = fabs(m[0]);
        m[1] = fabs(m[1]);
        m[2] = fabs(m[2]);
        m[3] = fabs(m[3]);
    }

     inline float max_elt() const {
        return fmaxf( fmaxf(m[0],m[1]), fmaxf(m[2],m[3]) );
    }

     inline float min_elt() const {
        return fminf( fminf(m[0],m[1]), fminf(m[2],m[3]) );
    }

    //--------------------------------------------------------------------------
    /// @name Static constructors
    //--------------------------------------------------------------------------


    static inline Mat2 identity()
    {
        return Mat2(1.f, 0.f,
                    0.f, 1.f);
    }

};

}// END Tbx NAMESPACE ==========================================================

#endif // TOOL_BOX_MAT2_HPP__
