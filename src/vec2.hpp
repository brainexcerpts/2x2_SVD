#ifndef TOOL_BOX_VEC2_HPP__
#define TOOL_BOX_VEC2_HPP__

#include <cmath>
#include <iostream>
#include <stdio.h>

// =============================================================================
namespace Tbx {
// =============================================================================


/** @brief 2D float vector type
*/
struct Vec2 {

    float x, y;

    // -------------------------------------------------------------------------
    /// @name Constructors
    // -------------------------------------------------------------------------

    Vec2() { x = 0.f; y = 0.f; }

    Vec2(float x_, float y_) { x = x_; y = y_; }

    Vec2(float v) { x = v; y = v; }

    static inline Vec2 unit_x(){ return Vec2(1.f, 0.f); }
    static inline Vec2 unit_y(){ return Vec2(0.f, 1.f); }
    static inline Vec2 zero()  { return Vec2(0.f, 0.f); }


    static inline Vec2 unit_scale(){ return Vec2(1.f, 1.f); }

    inline void set(float x_, float y_) { x = x_; y = y_; }

    static Vec2 random(float r){
        float r2 = 2.f * r;
        float x_ = rand() * 1.f /RAND_MAX;
        float y_ = rand() * 1.f /RAND_MAX;
        return Vec2(x_ * r2 - r, y_ * r2 - r);
    }

    // -------------------------------------------------------------------------
    /// @name Overload operators
    // -------------------------------------------------------------------------

    // ----------
    // Additions
    // ----------


    Vec2 operator+(const Vec2 &v_) const { return Vec2(x+v_.x, y+v_.y); }

    Vec2& operator+= (const Vec2 &v_) {
        x += v_.x;
        y += v_.y;
        return *this;
    }

    Vec2 operator+(float f_) const { return Vec2(x+f_, y+f_); }

     inline friend Vec2 operator+(const float d_, const Vec2& vec) { return Vec2(d_+vec.x, d_+vec.y); }

    Vec2& operator+= (float f_) {
        x += f_;
        y += f_;
        return *this;
    }

    // -------------
    // Substractions
    // -------------

    /// substraction
    Vec2 operator-(const Vec2 &v_) const { return Vec2(x-v_.x, y-v_.y); }

    Vec2& operator-= (const Vec2& v_) {
        x -= v_.x;
        y -= v_.y;
        return *this;
    }

    /// opposite vector
    Vec2 operator-() const { return Vec2(-x, -y); }

    Vec2 operator-(float f_) const { return Vec2(x-f_, y-f_); }

    inline friend Vec2 operator-(const float d_, const Vec2& vec) { return Vec2(d_-vec.x, d_-vec.y); }

    Vec2& operator-= (float f_) {
        x -= f_;
        y -= f_;
        return *this;
    }

    // -------------
    // Comparisons
    // -------------

    bool operator!= (const Vec2 &v_) const {
        return (x != v_.x) | (y != v_.y);
    }

    bool operator==(const Vec2& d_)  const {
        return (x == d_.x) && (y == d_.y);
    }

    // TODO operator "<" for maps

    // -------------
    // Divisions
    // -------------

    Vec2 operator/(const float d_) const {
        return Vec2(x/d_, y/d_);
    }

    Vec2& operator/=(const float d_) {
        x /= d_;
        y /= d_;
        return *this;
    }

    Vec2 operator/(const Vec2 &v_) const { return Vec2(x/v_.x, y/v_.y); }

    Vec2& operator/=(const Vec2& d_) {
        x /= d_.x;
        y /= d_.y;
        return *this;
    }

    // ----------------
    // Multiplication
    // ----------------

    /// rhs scalar multiplication
    Vec2 operator*(const float d_) const { return Vec2(x*d_, y*d_); }

    /// lhs scalar multiplication
     inline friend
    Vec2 operator*(const float d_, const Vec2& vec) { return Vec2(d_*vec.x, d_*vec.y); }

    Vec2& operator*=(const float d_) {
        x *= d_;
        y *= d_;
        return *this;
    }

    Vec2 operator*(const Vec2 &v_) const { return Vec2(x*v_.x, y*v_.y); }

    Vec2& operator*=(const Vec2& d_) {
        x *= d_.x;
        y *= d_.y;
        return *this;
    }

    // -------------------------------------------------------------------------
    /// @name Operators on vector
    // -------------------------------------------------------------------------

    /// @return vector perpendicular
    Vec2 perp() const { return Vec2(-y, x); }

    /// product of all components
    float product() const { return x*y; }

    /// sum of all components
    float sum() const { return x+y; }

    /// Average all components
    float average() const { return (x+y) * 0.5f; }

    /// semi dot product
    Vec2 mult(const Vec2& v) const {
        return Vec2(x*v.x, y*v.y);
    }

    /// dot product
    float dot(const Vec2 &v_) const {
        return x * v_.x + y * v_.y;
    }

    /// @return signed angle between [-PI; PI] starting from 'this' to 'v_'
    float signed_angle(const Vec2 &v_) const {
        return atan2( x * v_.y - y * v_.x, x * v_.x + y * v_.y );
    }

    /// absolute value of the dot product
    float abs_dot(const Vec2 &v_) const {
        return fabsf(x * v_.x + y * v_.y);
    }

    /// norm squared
    float norm_squared() const {
        return dot(*this);
    }

    /// normalization
    Vec2 normalized() const {
        return (*this) * (1.f/sqrtf(norm_squared()));
    }

    /// normalization
    float normalize() {
        float l = sqrtf(norm_squared());
        float f = 1.f / l;
        x *= f;
        y *= f;
        return l;
    }

    /// normalization
    float safe_normalize(const float eps = 1e-10f) {
        float l = sqrtf(norm_squared());
        if(l > eps){
            float f = 1.f / l;
            x *= f;
            y *= f;
            return l;
        } else {
            x = 1.f;
            y = 0.f;
            return 0.f;
        }
    }

    /// norm
    float norm() const {
        return sqrtf(norm_squared());
    }

    /// value of the min coordinate
    float get_min() const {
        return fminf(x,y);
    }

    /// value of the max coordinate
    float get_max() const {
        return fmaxf(x,y);
    }

    /// clamp each vector values
    Vec2 clamp(float min_v, float max_v) const {
        return Vec2( fminf( fmaxf(x, min_v), max_v),
                     fminf( fmaxf(y, min_v), max_v));
    }

    /// clamp each vector values
    Vec2 clamp(const Vec2& min_v, const Vec2& max_v) const {
        return Vec2( min( max(x, min_v.x), max_v.x),
                     min( max(y, min_v.y), max_v.y));
    }

    /// floorf every components
    Vec2 floor() const {
        return Vec2( floorf(x), floorf(y) );
    }

    /// rotate of 0 step to the left (present for symmetry)
    Vec2 perm_x() const {
        return Vec2(x, y);
    }

    /// rotate of 1 step to the left (so that y is the first coordinate)
    Vec2 perm_y() const {
        return Vec2(y, x);
    }

    /// @return the vector to_project projected on the line defined by the
    /// direction '*this'
    /// @warning don't forget to normalize the vector before calling this !
    Vec2 proj_on_line(const Vec2& to_project) const {
        return (*this) * (*this).dot( to_project );
    }

    // -------------------------------------------------------------------------
    /// @name Accessors
    // -------------------------------------------------------------------------


    inline const float& operator[](int i) const{
        assert( i < 2);
        return ((float*)this)[i];
    }


    inline float& operator[](int i){
        assert( i < 2);
        return ((float*)this)[i];
    }

    /// Conversion returns the memory address of the vector.
    /// Very convenient to pass a Vec pointer as a parameter to OpenGL:
    /// @code
    /// Vec2 pos, normal;
    /// glNormal2fv(normal);
    /// glVertex2fv(pos);
    /// @endcode
     operator const float*() const { return &x; }

    /// Conversion returns the memory address of the vector. (Non const version)
     operator float*() { return &x; }

    // -------------------------------------------------------------------------
    /// @name Print vector
    // -------------------------------------------------------------------------

    inline void print() const {
        printf("%f, %f\n", x, y);
    }

    inline friend
    std::ostream& operator<< ( std::ostream& ofs, const Vec2& v2 )
    {
        ofs << v2.x << ", " << v2.y << "; ";
        return ofs;
    }
};

}// END Tbx NAMESPACE ==========================================================


#endif // TOOL_BOX_VEC2_HPP__
