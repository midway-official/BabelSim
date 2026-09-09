#pragma once

#include <cmath>
#include <cstddef>
#include <ostream>

namespace babelsim {

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    constexpr double& operator[](std::size_t i) { return (&x)[i]; }
    constexpr double operator[](std::size_t i) const { return (&x)[i]; }

    constexpr Vec3& operator+=(const Vec3& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }
    constexpr Vec3& operator-=(const Vec3& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }
    constexpr Vec3& operator*=(double scale) {
        x *= scale;
        y *= scale;
        z *= scale;
        return *this;
    }
    constexpr Vec3& operator/=(double scale) {
        x /= scale;
        y /= scale;
        z /= scale;
        return *this;
    }
};

// 向量场梯度按分量存储：rows[0]=grad(Ux) 等。三个 Vec3 行使每个 cell 的 tensor 连续，
// 不需要 Eigen 对齐或专用分配器。
struct Tensor3 {
    Vec3 rows[3]{};

    constexpr Vec3& operator[](std::size_t i) { return rows[i]; }
    constexpr const Vec3& operator[](std::size_t i) const { return rows[i]; }
};

constexpr Vec3 operator+(Vec3 lhs, const Vec3& rhs) { return lhs += rhs; }
constexpr Vec3 operator-(Vec3 lhs, const Vec3& rhs) { return lhs -= rhs; }
constexpr Vec3 operator-(const Vec3& value) {
    return {-value.x, -value.y, -value.z};
}
constexpr Vec3 operator*(Vec3 value, double scale) { return value *= scale; }
constexpr Vec3 operator*(double scale, Vec3 value) { return value *= scale; }
constexpr Vec3 operator/(Vec3 value, double scale) { return value /= scale; }

constexpr double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

constexpr Vec3 cross(const Vec3& a, const Vec3& b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    };
}

constexpr double squaredNorm(const Vec3& value) { return dot(value, value); }
inline double norm(const Vec3& value) { return std::sqrt(squaredNorm(value)); }
inline bool isFinite(const Vec3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y) &&
        std::isfinite(value.z);
}

inline Tensor3 operator+(Tensor3 a, const Tensor3& b) {
    for (int i = 0; i < 3; ++i) a[i] += b[i];
    return a;
}
inline Tensor3 operator-(Tensor3 a, const Tensor3& b) {
    for (int i = 0; i < 3; ++i) a[i] -= b[i];
    return a;
}
inline Tensor3 operator*(double a, Tensor3 b) {
    for (int i = 0; i < 3; ++i) b[i] *= a;
    return b;
}
inline Tensor3& operator+=(Tensor3& a, const Tensor3& b) { return a = a + b; }
inline Vec3 dot(const Tensor3& tensor, const Vec3& vector) {
    return {dot(tensor[0], vector), dot(tensor[1], vector), dot(tensor[2], vector)};
}
inline Tensor3 transpose(const Tensor3& a) {
    Tensor3 result;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) result[i][j] = a[j][i];
    return result;
}
inline double trace(const Tensor3& a) { return a[0][0] + a[1][1] + a[2][2]; }
inline double symmetricBoundaryValue(double value, Vec3) { return value; }
inline Vec3 symmetricBoundaryValue(Vec3 value, Vec3 normal) {
    return value - dot(value, normal) * normal;
}
inline Tensor3 symmetricBoundaryValue(const Tensor3& value, Vec3 n) {
    Tensor3 reflected;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    reflected[i][j] += ((i == k ? 1.0 : 0.0) - 2*n[i]*n[k]) *
                        value[k][l] * ((l == j ? 1.0 : 0.0) - 2*n[l]*n[j]);
    return 0.5 * (value + reflected);
}

inline std::ostream& operator<<(std::ostream& output, const Vec3& value) {
    return output << '(' << value.x << ", " << value.y << ", " << value.z << ')';
}

}  // babelsim 命名空间
