#pragma once

#include "babelsim/field.h"

namespace babelsim::detail {
// 存储维护接口，不安装到 Solver SDK。它只提供热路径所需的借用视图，不改变容量。
struct FieldAccess {
    template <typename T> static BoundaryCondition<T> condition(const Field<T>& field, Index face,
        double flux = std::numeric_limits<double>::quiet_NaN()) { return field.faceCondition(face, flux); }
    template <typename T> static T trace(const Field<T>& field, Index face,
        double flux = std::numeric_limits<double>::quiet_NaN()) { return field.boundaryTrace(face, flux); }
    template <typename T> static void setTrace(Field<T>& field, Index face, T value) {
        field.m_boundary_trace.at(face) = value;
    }
    // 无输入边界微分约束的单元结果使用相邻单元常值外推；调用前同步单元值。
    template <typename T> static void extrapolateTrace(Field<T>& field) {
        if (!field.m_calculated_boundary) return;
        for (Index face = 0; face < field.m_mesh->faceCount(); ++face)
            if (field.m_mesh->boundaryFace(face))
                field.m_boundary_trace[face] = field.m_values[field.m_mesh->owner(face)];
    }
    template <typename T> static T* data(Field<T>& field) {
        field.m_halo_valid = false;
        return field.m_values.data();
    }
    template <typename T> static const T* data(const Field<T>& field) { return field.m_values.data(); }
    template <typename T> static const std::vector<T>& values(const Field<T>& field) { return field.m_values; }
    template <typename T> static bool haloValid(const Field<T>& field) {
        return field.m_halo_valid;
    }
    template <typename T> static void markHaloValid(Field<T>& field) {
        field.m_halo_valid = true;
    }
};
template <typename T> T* fieldData(Field<T>& field) { return FieldAccess::data(field); }
template <typename T> const T* fieldData(const Field<T>& field) { return FieldAccess::data(field); }
template <typename T> const std::vector<T>& fieldValues(const Field<T>& field) { return FieldAccess::values(field); }
template <typename T> bool haloValid(const Field<T>& field) { return FieldAccess::haloValid(field); }
template <typename T> void markHaloValid(Field<T>& field) { FieldAccess::markHaloValid(field); }
}  // babelsim::detail 命名空间
