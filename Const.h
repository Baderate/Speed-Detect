#pragma once

// 圆周率，统一用于角度/弧度换算。
constexpr double PI = 3.141592653589793238462643383279502884;

// WGS84 椭球长半轴，单位 m。
constexpr double EARTH_SEMI_MAJOR = 6378137.0;

// WGS84 第一偏心率平方。
constexpr double EARTH_FIRST_ECCENTRICITY_SQUARED = 0.00669437999013;

// 算法内部历史变量名，保留为常量别名，避免公式里混入多套数值。
constexpr double earth_R = EARTH_SEMI_MAJOR;
constexpr double earth_e2 = EARTH_FIRST_ECCENTRICITY_SQUARED;

