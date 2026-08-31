#pragma once

#include "babelsim/field.h"

namespace babelsim::detail {
// 存储维护接口，不安装到 Solver SDK。它只提供热路径所需的借用视图，不改变容量。
struct FieldAccess {
    template <typename T> static T* data(Field<T>& field) { return field.m_values.data(); }
    template <typename T> static const T* data(const Field<T>& field) { return field.m_values.data(); }
    template <typename T> static const std::vector<T>& values(const Field<T>& field) { return field.m_values; }
};
template <typename T> T* fieldData(Field<T>& field) { return FieldAccess::data(field); }
template <typename T> const T* fieldData(const Field<T>& field) { return FieldAccess::data(field); }
template <typename T> const std::vector<T>& fieldValues(const Field<T>& field) { return FieldAccess::values(field); }
}  // babelsim::detail 命名空间
