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

inline std::ostream& operator<<(std::ostream& output, const Vec3& value) {
    return output << '(' << value.x << ", " << value.y << ", " << value.z << ')';
}

}  // babelsim 命名空间
