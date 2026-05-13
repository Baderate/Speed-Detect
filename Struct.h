#pragma once

#include <limits>
#include <string>
#include <vector>

#include "Const.h"

inline double NaNValue()
{
    return std::numeric_limits<double>::quiet_NaN();
}

// ======================== 全项目共享结构 ========================

// SD 主流程模式。后续增加新的 SD 算法时，只扩展这里和 SD.cpp 内部分发逻辑。
enum class SDMode
{
    Off = 0,                         // 关闭 SD。
    ThreeSigma = 1,                  // 基于 3sigma 的统计异常检测。
    SpeedPositionPrediction = 2      // 用前一时刻速度和位置预测当前点，再与 GNSS 比较。
};

// 统一车辆状态点：后续 SD、MM、HMM 和结果输出都以它作为最终表达。
struct CarPoint
{
    long long timestamp = 0;             // UTC 时间戳，单位 ms。
    double BLH[3] = { NaNValue(), NaNValue(), NaNValue() };
    double x = NaNValue();               // 局部平面坐标 x，单位 m。
    double y = NaNValue();               // 局部平面坐标 y，单位 m。
    double V_ENU[3] = { NaNValue(), NaNValue(), NaNValue() };
    double quality_score = 1.0;
    bool trusted = true;
    int trust_reason = 0;
};

// coor 文件中的原始 GNSS 观测。字段顺序尽量贴合文件列，避免读取时丢信息。
struct CoorData
{
    int year = 0;
    int month = 0;
    int day = 0;
    int hour = 0;
    int minute = 0;
    double second = 0.0;

    int week = 0;
    double ws = 0.0;
    int nSV = 0;
    double sigma0 = 0.0;

    double BLH[3] = { 0.0, 0.0, 0.0 };
    double V_ENU[3] = { 0.0, 0.0, 0.0 };

    double velU = 0.0;
    double velN = 0.0;
    double velE = 0.0;

    double phaU = 0.0;
    double phaN = 0.0;
    double phaE = 0.0;

    double sigma[2] = { 0.0, 0.0 };
    double sigma1 = 0.0;
    double sigma2 = 0.0;
    int LObs = 0;
    double clk1 = 0.0;
    double trop = 0.0;

    double CBP2[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    double ISB[4] = { 0.0, 0.0, 0.0, 0.0 };

    // coor 后续新增列先保留在这里，避免解析器因格式扩展失效。
    std::vector<double> extraCols;

    bool trusted = true;
    int trust_reason = 0;
    double quality_score = 1.0;
};

// txt 真值或匹配结果点。redius 保留原数据拼写，避免和已有文件字段脱节。
struct TxtData
{
    std::string dataType;
    long long timestamp = 0;             // UTC 时间戳，单位 ms。
    double BLH[3] = { NaNValue(), NaNValue(), NaNValue() };
    double redius = NaNValue();
};

// ======================== 路网共享结构 ========================
// 这些结构由 File / Road / MM / HMM / Util 共用，因此保留在 Struct.h 中。

// 路网折线点，BLH 顺序为纬度、经度、高程。
struct RoadNode
{
    double BLH[3] = { NaNValue(), NaNValue(), NaNValue() };
};

// 一条道路中心线或虚拟车道线。
struct RoadPolyline
{
    int id = 0;
    int laneId = 0;
    double laneOffset = 0.0;
    std::string name;
    std::vector<RoadNode> points;
};

// 路网集合。sourcePath 记录可直接交给算法核心解析的 KML 路径。
struct RoadNetwork
{
    std::string sourcePath;
    bool hasVirtualLanes = false;
    std::vector<RoadPolyline> polylines;
};

// 小范围经纬度和平面坐标互转的局部参考原点。
struct LocalRef
{
    double lat0_rad = 0.0;
    double lon0_rad = 0.0;
    double cosLat0 = 1.0;
};

// ======================== EKF 配置结构 ========================
// EKF 状态为 [x, y, velE, velN]，GNSS 作为主观测，MM/HMM 作为道路伪观测。
struct EKFOptions
{
    double initialPositionStdM = 8.0;        // 初始位置标准差，单位 m。
    double initialVelocityStdMps = 3.0;      // 初始速度标准差，单位 m/s。

    double processPositionNoiseM = 0.8;      // 位置过程噪声，单位 m。
    double processVelocityNoiseMps = 0.6;    // 速度过程噪声，单位 m/s。

    double gnssPositionNoiseM = 6.0;         // GNSS 位置观测标准差，单位 m。
    double gnssVelocityNoiseMps = 1.5;       // GNSS 速度观测标准差，单位 m/s。
    double roadPositionNoiseM = 2.0;         // MM/HMM 道路伪观测标准差，单位 m。

    double minQuality = 0.05;                // 质量评分下限，避免噪声无限放大。
    double untrustedNoiseScale = 25.0;       // SD 标记不可信时的噪声放大倍数。
    long long maxPseudoTimeDiffMs = 1000;    // 伪观测和 GNSS 的最大时间差，单位 ms。

    bool useGnssVelocity = true;             // 是否使用 GNSS 速度作为速度观测。
    bool useRoadPseudo = true;               // 是否使用 MM/HMM 结果作为道路伪观测。
};

// 算法内部使用的道路段平面表示。
struct RoadSegmentXY
{
    double lat1 = 0.0;
    double lon1 = 0.0;
    double lat2 = 0.0;
    double lon2 = 0.0;
    double x1 = 0.0;
    double y1 = 0.0;
    double x2 = 0.0;
    double y2 = 0.0;
    int id = 0;
    int clusterId = 0;
    double heading = 0.0;
};

// ======================== SD-MM 共享结构 ========================
// 这些结构只服务于 SD.cpp 和 MM.cpp 之间的接口传递，因此单独分组说明。

enum SDMMStage
{
    SD_MM_STAGE_EPOCH = 1,
    SD_MM_STAGE_CANDIDATE = 2,
    SD_MM_STAGE_TRANSITION = 3,
    SD_MM_STAGE_UPDATE = 4
};

struct SDMMRequest
{
    int stage = SD_MM_STAGE_EPOCH;
    bool use_sd = false;
    int mode = 0;

    bool has_prev = false;
    bool has_speed_meas = false;
    bool has_vel = false;
    bool is_static_epoch = false;
    bool use_vel_heading = false;
    bool has_pred_step = false;

    double x_obs = 0.0;
    double y_obs = 0.0;
    double x_dr = 0.0;
    double y_dr = 0.0;
    double prev_x = 0.0;
    double prev_y = 0.0;
    double vx_last = 0.0;
    double vy_last = 0.0;
    double vE_meas = 0.0;
    double vN_meas = 0.0;
    double speed_meas = 0.0;
    double dt = 1.0;

    double W_OBS = 0.0;
    double W_DR = 0.0;
    double W_CONT = 0.0;
    double W_VEL_DIR = 0.0;
    double road_heading = 0.0;

    double route_dist = -1.0;
    double pred_step = 0.0;
    double fallback_step = 0.0;
    double trans_pred_beta = 12.0;

    double best_x = 0.0;
    double best_y = 0.0;
    double old_prev_x = 0.0;
    double old_prev_y = 0.0;
};

struct SDMMResult
{
    bool enabled = false;
    bool has_speed = false;
    bool is_static_epoch = false;
    bool use_vel_heading = false;

    double vE = 0.0;
    double vN = 0.0;
    double speed = 0.0;

    double x_dr = 0.0;
    double y_dr = 0.0;
    double x_center = 0.0;
    double y_center = 0.0;

    bool trusted = true;
    int trust_reason = 0;
    double quality_score = 1.0;

    double W_OBS = 0.0;
    double W_DR = 0.0;
    double W_CONT = 0.0;
    double W_VEL_DIR = 0.0;
    double vel_dir_cost = 0.0;

    double pred_step = 0.0;
    double trans_pred_cost = 0.0;

    double vx_last = 0.0;
    double vy_last = 0.0;
    bool has_vel = false;
};

// ======================== 单模块内部结构说明 ========================
// 以下结构只在单个模块中使用，因此不放在 Struct.h，避免公共头不断膨胀：
// 1. Road.cpp：
//    TrailPointLL、TrailPointXY。
// 2. MM.cpp：
//    GraphEdge、RoadSegmentExt、SeqCandidate、ProvisionalState、EpochMeta。
// 3. HMM.cpp：
//    GraphEdge、SegExt、HMMCandidate、ChainRange。
// 4. SD.cpp：
//    SDMMResidualStats、SDMMResidualJudge。


