#include "Head.h"

namespace
{
bool IsInExcludeWindows(
    long long t_ms,
    const std::vector<std::pair<long long, long long>>& excludeWindows)
{
    for (const auto& win : excludeWindows) {
        if (t_ms >= win.first && t_ms <= win.second) {
            return true;
        }
    }
    return false;
}

bool IsInKeepWindow(long long t_ms, long long startMs, long long endMs)
{
    return t_ms >= startMs && t_ms <= endMs;
}

LocalRef BuildLocalRefFromCoords(const std::vector<CoorData>& coords)
{
    double sumLat = 0.0;
    double sumLon = 0.0;
    int count = 0;

    for (const auto& c : coords) {
        if (std::isfinite(c.BLH[0]) && std::isfinite(c.BLH[1])) {
            sumLat += c.BLH[0];
            sumLon += c.BLH[1];
            ++count;
        }
    }

    LocalRef ref;
    if (count > 0) {
        ref.lat0_rad = (sumLat / count) * PI / 180.0;
        ref.lon0_rad = (sumLon / count) * PI / 180.0;
        ref.cosLat0 = std::cos(ref.lat0_rad);
    }
    return ref;
}

void BlhToXY(double latDeg, double lonDeg, const LocalRef& ref, double& x, double& y)
{
    const double latRad = latDeg * PI / 180.0;
    const double lonRad = lonDeg * PI / 180.0;
    y = earth_R * (latRad - ref.lat0_rad);
    x = earth_R * (lonRad - ref.lon0_rad) * ref.cosLat0;
}
}

/**
 * 函数功能：
 *     计算两个经纬度点之间的球面距离。
 *
 * 输入：
 *     lat1、lon1：第一个点的纬度和经度，单位 deg。
 *     lat2、lon2：第二个点的纬度和经度，单位 deg。
 *
 * 输出：
 *     double：两点间距离，单位 m。
 *
 * 调用位置：
 *     SD.cpp 的速度检测评估、真值误差统计会调用。
 *
 * 数据流作用：
 *     BLH 坐标 -> 空间距离误差。
 */
double haversine(double lat1, double lon1, double lat2, double lon2)
{
    const double phi1 = lat1 * PI / 180.0;
    const double phi2 = lat2 * PI / 180.0;
    const double deltaPhi = (lat2 - lat1) * PI / 180.0;
    const double deltaLambda = (lon2 - lon1) * PI / 180.0;

    const double a = std::sin(deltaPhi / 2.0) * std::sin(deltaPhi / 2.0)
        + std::cos(phi1) * std::cos(phi2)
        * std::sin(deltaLambda / 2.0) * std::sin(deltaLambda / 2.0);
    const double c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
    return earth_R * c;
}

/**
 * 函数功能：
 *     将 GPS 周和周内秒转换为 UTC 毫秒时间戳。
 *
 * 输入：
 *     gpsWeek：GPS 周。
 *     gpsWeekSecond：GPS 周内秒。
 *
 * 输出：
 *     long long：UTC 时间戳，单位 ms。
 *
 * 调用位置：
 *     ConvertCoorToCarPoint、FilterCoordsByExcludeWindows、ComputeCoorTimestampSec 调用。
 *
 * 数据流作用：
 *     coor 文件中的 GPS 时间 -> 全项目统一 UTC 毫秒时间。
 */
long long GpsTimeToUtcTimestampMs(int gpsWeek, double gpsWeekSecond)
{
    const double unixSeconds = gpsWeek * 604800.0 + gpsWeekSecond + 315964800.0;
    return static_cast<long long>(std::llround(unixSeconds * 1000.0));
}

/**
 * 函数功能：
 *     将 CoorData 中的 GPS 时间转换为 UTC 秒级时间戳。
 *
 * 输入：
 *     coor：GNSS 原始观测点。
 *
 * 输出：
 *     long long：UTC 时间戳，单位 s。
 *
 * 调用位置：
 *     MM.cpp、HMM.cpp、SD.cpp 中用于计算相邻历元时间间隔。
 *
 * 数据流作用：
 *     CoorData 时间字段 -> 匹配算法中的 dt / 分链判断 / 输出时间戳。
 */
long long ComputeCoorTimestampSec(const CoorData& coor)
{
    return static_cast<long long>(GpsTimeToUtcTimestampMs(coor.week, coor.ws) / 1000LL);
}

/**
 * 函数功能：
 *     删除已知异常时间窗内的 GNSS 点。
 *
 * 输入：
 *     coords：原始 GNSS 观测数组。
 *
 * 输出：
 *     std::vector<CoorData>：删除异常时间窗后的 GNSS 观测数组。
 *
 * 调用位置：
 *     main.cpp 的状态机 MM 模式和 HMM 模式在匹配前调用。
 *
 * 数据流作用：
 *     原始 CoorData -> 过滤后的 CoorData，避免固定异常片段影响地图匹配。
 */
std::vector<CoorData> FilterCoordsByExcludeWindows(const std::vector<CoorData>& coords)
{
    const long long keepStartMs = 1678329750000LL;
    const long long keepEndMs = 1678341540000LL;

    const std::vector<std::pair<long long, long long>> excludeWindows =
    {
        {1678330500000LL, 1678330620000LL},
        {1678331790000LL, 1678331850000LL},
        {1678333560000LL, 1678333680000LL},
        {1678340400000LL, 1678341000000LL}
    };

    std::vector<CoorData> filtered;
    filtered.reserve(coords.size());

    int invalidTimeCount = 0;
    int outsideKeepWindowCount = 0;
    int excludedWindowCount = 0;
    for (const auto& pt : coords) {
        const long long utcMs = GpsTimeToUtcTimestampMs(pt.week, pt.ws);
        if (utcMs < 0) {
            ++invalidTimeCount;
            continue;
        }
        if (!IsInKeepWindow(utcMs, keepStartMs, keepEndMs)) {
            ++outsideKeepWindowCount;
            continue;
        }
        if (IsInExcludeWindows(utcMs, excludeWindows)) {
            ++excludedWindowCount;
            continue;
        }
        filtered.push_back(pt);
    }

    const int removedCount =
        invalidTimeCount + outsideKeepWindowCount + excludedWindowCount;
    std::cout << "时间过滤完成：大时间段外删除 = " << outsideKeepWindowCount
        << "，小异常窗口删除 = " << excludedWindowCount
        << "，无效时间删除 = " << invalidTimeCount
        << "，保留点数量 = " << filtered.size() << std::endl;

    std::cout << "硬排除完成：删除点数量 = " << removedCount
        << "，保留点数量 = " << filtered.size() << std::endl;

    return filtered;
}

/**
 * 函数功能：
 *     将单个 GNSS 原始点转换为统一车辆状态点 CarPoint。
 *
 * 输入：
 *     coor：GNSS 原始观测点。
 *     ref：局部平面坐标参考原点；为空时不计算 x/y。
 *
 * 输出：
 *     CarPoint：包含 UTC 毫秒时间、BLH、ENU 速度、质量评分和可选 x/y。
 *
 * 调用位置：
 *     ConvertCoorToCarPoints 内部逐点调用。
 *
 * 数据流作用：
 *     CoorData -> CarPoint，是 SD、输出和最终结果展示的统一格式转换。
 */
CarPoint ConvertCoorToCarPoint(const CoorData& coor, const LocalRef* ref)
{
    CarPoint carPoint;
    carPoint.timestamp = GpsTimeToUtcTimestampMs(coor.week, coor.ws);
    carPoint.BLH[0] = coor.BLH[0];
    carPoint.BLH[1] = coor.BLH[1];
    carPoint.BLH[2] = coor.BLH[2];
    carPoint.V_ENU[0] = coor.velE;
    carPoint.V_ENU[1] = coor.velN;
    carPoint.V_ENU[2] = coor.velU;
    carPoint.quality_score = coor.quality_score;
    carPoint.trusted = coor.trusted;
    carPoint.trust_reason = coor.trust_reason;

    if (ref != nullptr) {
        BlhToXY(carPoint.BLH[0], carPoint.BLH[1], *ref, carPoint.x, carPoint.y);
    }

    return carPoint;
}

/**
 * 函数功能：
 *     将 GNSS 原始点序列转换为车辆状态点序列。
 *
 * 输入：
 *     coords：GNSS 原始观测数组。
 *     ref：可选局部平面坐标参考原点；为空时按轨迹均值自动建立。
 *
 * 输出：
 *     std::vector<CarPoint>：统一车辆状态点数组。
 *
 * 调用位置：
 *     1. main.cpp 的 SD 模式生成检测输入。
 *     2. HMM.cpp 启用 SD 时先转换为 CarPoint 做速度检测。
 *
 * 数据流作用：
 *     CoorData 数组 -> CarPoint 数组，供 SD 和结果输出使用。
 */
std::vector<CarPoint> ConvertCoorToCarPoints(const std::vector<CoorData>& coords, const LocalRef* ref)
{
    LocalRef localRef;
    const LocalRef* refToUse = ref;
    if (refToUse == nullptr) {
        localRef = BuildLocalRefFromCoords(coords);
        refToUse = &localRef;
    }

    std::vector<CarPoint> points;
    points.reserve(coords.size());
    for (const auto& c : coords) {
        points.push_back(ConvertCoorToCarPoint(c, refToUse));
    }
    return points;
}



