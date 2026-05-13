#pragma once

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Const.h"
#include "Struct.h"

// ======================== File.cpp：文件读写 ========================
/**
 * @brief 读取 coor 格式 GNSS 原始轨迹文件。
 * @param filename coor 文件路径。
 * @return 解析后的 GNSS 原始观测序列；文件打开失败或格式无效时返回空数组。
 */
std::vector<CoorData> ReadCoorFile(const std::string& filename);

/**
 * @brief 读取 txt 格式参考真值文件。
 * @param filename txt 文件路径。
 * @return 参考真值点序列；缺失字段保留为 NaN。
 */
std::vector<TxtData> ReadTxtFile(const std::string& filename);

/**
 * @brief 读取轨迹文本文件，并直接转换为 CarPoint 序列。
 * @param filename 轨迹文件路径，支持 MM/HMM 匹配结果文件和 CarPoint 输出文件两类格式。
 * @return 车辆状态点序列；若输入文件没有 x/y，则函数内部会自动补齐局部平面坐标。
 */
std::vector<CarPoint> ReadCarTrajectory(const std::string& filename);

/**
 * @brief 读取地图匹配结果文件，并转换为 TxtData 结构。
 * @param filename 地图匹配结果文件路径。
 * @param defaultDataType 文件行未显式给出类型时使用的数据类型。
 * @param defaultRedius 文件行未显式给出半径时使用的默认半径。
 * @return 可用于误差评估或 EKF 对齐的匹配点序列。
 */
std::vector<TxtData> ReadRoadMatchedFile(
    const std::string& filename,
    const std::string& defaultDataType = "GPS_matched",
    double defaultRedius = 0.5);

/**
 * @brief 读取 KML 路网中心线。
 * @param filename KML 文件路径。
 * @return RoadNetwork 路网结构，sourcePath 会记录输入文件路径。
 */
RoadNetwork ReadKmlRoadNetwork(const std::string& filename);

/**
 * @brief 输出统一车辆轨迹文件。
 * @param filename 输出文件路径。
 * @param points 待输出的车辆状态点序列。
 * @return 写入成功返回 true，否则返回 false。
 */
bool WriteCarTrajectory(const std::string& filename, const std::vector<CarPoint>& points);

/**
 * @brief 将 RoadNetwork 写为 KML 折线文件。
 * @param filename 输出 KML 文件路径。
 * @param road 待输出的路网或虚拟车道路网。
 * @return 写入成功返回 true，否则返回 false。
 */
bool WriteRoadNetworkToKml(const std::string& filename, const RoadNetwork& road);

/**
 * @brief 将 RoadNetwork 写为调试用 TXT 文件。
 * @param filename 输出 TXT 文件路径。
 * @param road 待输出的路网或虚拟车道路网。
 * @return 写入成功返回 true，否则返回 false。
 */
bool WriteRoadNetworkToTxt(const std::string& filename, const RoadNetwork& road);

// ======================== Util.cpp：通用工具 ========================
/**
 * @brief 使用 Haversine 公式计算两个经纬度点的球面距离。
 * @param lat1 第一个点纬度，单位 deg。
 * @param lon1 第一个点经度，单位 deg。
 * @param lat2 第二个点纬度，单位 deg。
 * @param lon2 第二个点经度，单位 deg。
 * @return 两点间距离，单位 m。
 */
double haversine(double lat1, double lon1, double lat2, double lon2);

/**
 * @brief 将 GPS 周和周内秒转换为 UTC 毫秒时间戳。
 * @param gpsWeek GPS 周。
 * @param gpsWeekSecond GPS 周内秒。
 * @return UTC 时间戳，单位 ms。
 */
long long GpsTimeToUtcTimestampMs(int gpsWeek, double gpsWeekSecond);

/**
 * @brief 计算 CoorData 对应的 UTC 秒级时间戳。
 * @param coor GNSS 原始观测点。
 * @return UTC 时间戳，单位 s。
 */
long long ComputeCoorTimestampSec(const CoorData& coor);

/**
 * @brief 删除固定异常时间窗内的 GNSS 观测点。
 * @param coords 原始 GNSS 观测序列。
 * @return 删除异常时间窗后的 GNSS 观测序列。
 */
std::vector<CoorData> FilterCoordsByExcludeWindows(const std::vector<CoorData>& coords);

/**
 * @brief 将单个 CoorData 转换为统一 CarPoint。
 * @param coor GNSS 原始观测点。
 * @param ref 可选局部坐标参考原点；为空时只填 BLH、速度和质量字段。
 * @return 转换后的车辆状态点。
 */
CarPoint ConvertCoorToCarPoint(const CoorData& coor, const LocalRef* ref = nullptr);

/**
 * @brief 将 GNSS 原始观测序列转换为统一 CarPoint 序列。
 * @param coords GNSS 原始观测序列。
 * @param ref 可选局部坐标参考原点；为空时根据输入轨迹自动估计。
 * @return 转换后的车辆状态点序列。
 */
std::vector<CarPoint> ConvertCoorToCarPoints(const std::vector<CoorData>& coords, const LocalRef* ref = nullptr);

// ======================== Road.cpp：路网和虚拟车道 ========================
/**
 * @brief 基于道路中心线生成虚拟车道路网。
 * @param inputRoad 原始路网中心线。
 * @param useVirtualLane 是否生成虚拟车道；为 false 时返回原路网副本。
 * @return 与输入同结构的 RoadNetwork，hasVirtualLanes 标记是否包含虚拟车道。
 */
RoadNetwork GenerateVirtualLaneNetwork(const RoadNetwork& inputRoad, bool useVirtualLane = true);

/**
 * @brief 从输入 KML 直接生成虚拟车道 KML 和 TXT 文件。
 * @param inputKmlFile 原始路网 KML 文件路径。
 * @param outputKmlFile 输出虚拟车道 KML 文件路径。
 * @param outputTxtFile 输出虚拟车道 TXT 文件路径。
 * @return 两个输出文件都写入成功时返回 true。
 */
bool GenerateVirtualLaneFiles(
    const std::string& inputKmlFile,
    const std::string& outputKmlFile,
    const std::string& outputTxtFile);

/**
 * @brief 获取当前虚拟车道配置下的横向偏移数量。
 * @param useVirtualLane 是否启用虚拟车道。
 * @return 启用时返回虚拟车道数量；关闭时返回 1。
 */
int GetVirtualLaneOffsetCount(bool useVirtualLane);

/**
 * @brief 获取指定虚拟车道的横向偏移。
 * @param useVirtualLane 是否启用虚拟车道。
 * @param laneIndex 车道索引，越界时会被截断到有效范围。
 * @return 横向偏移距离，单位 m。
 */
double GetVirtualLaneOffset(bool useVirtualLane, int laneIndex);

// ======================== SD.cpp：速度检测 ========================
/**
 * @brief SD 主流程入口，并直接修改 CarPoint 的质量字段。
 * @param points 车辆状态点序列，函数会修改 trusted、trust_reason、quality_score。
 * @param mode SD 算法模式，目前支持 3sigma 和前一时刻速度位置预测。
 */
void RunSpeedDetection(std::vector<CarPoint>& points, SDMode mode);

/**
 * @brief SD 验证入口，内部调用 SD 主流程，再与参考真值对齐输出评估指标。
 * @param coords GNSS 原始观测序列。
 * @param truthdata 参考真值点序列。
 * @param mode SD 算法模式，目前支持 3sigma 和前一时刻速度位置预测。
 * @param truthAbnormalThresholdM 真值异常阈值，单位 m。
 * @param debugOutputFile 调试明细输出路径；为空时不输出调试文件。
 */
void EvaluateSpeedDetection(
    const std::vector<CoorData>& coords,
    const std::vector<TxtData>& truthdata,
    SDMode mode = SDMode::SpeedPositionPrediction,
    double truthAbnormalThresholdM = 5.0,
    const std::string& debugOutputFile = "");

/**
 * @brief 兼容旧调用的独立 SD 召回率验证流程，内部转调 EvaluateSpeedDetection。
 * @param coords GNSS 原始观测序列。
 * @param truthdata 参考真值点序列。
 * @param sdThresholdM 旧接口保留参数，当前主流程阈值由 SD 模式内部控制。
 * @param truthAbnormalThresholdM 真值误差异常阈值，单位 m。
 * @param debugOutputFile 调试明细输出路径。
 */
void RunSDModule_RecallOnly(
    const std::vector<CoorData>& coords,
    const std::vector<TxtData>& truthdata,
    double sdThresholdM,
    double truthAbnormalThresholdM,
    const std::string& debugOutputFile);

/**
 * @brief MM/HMM 内部使用的速度辅助评价接口。
 * @param c 当前 GNSS 点；非历元阶段可传 nullptr。
 * @param req 本次评价请求，包含阶段、速度、候选、转移等上下文。
 * @param out 输出评价结果，包括搜索中心、质量评分、代价修正等。
 * @return SD 启用且完成对应阶段评价时返回 true；关闭时返回 false。
 */
bool SD_MMEvaluate(const CoorData* c, const SDMMRequest& req, SDMMResult& out);

// ======================== MM.cpp：状态机地图匹配 ========================
/**
 * @brief 状态机地图匹配主接口。
 * @param road 路网数据；sourcePath 指向实际供算法解析的 KML 文件。
 * @param gnssData GNSS 原始观测序列。
 * @param enableSD 是否启用 SD 辅助。
 * @param sdMode SD 模式位掩码。
 * @param outputFile 匹配结果 TXT 输出路径，同时会生成同名 KML。
 * @return 匹配后的车辆状态点序列。
 */
std::vector<CarPoint> MapMarching(
    const RoadNetwork& road,
    const std::vector<CoorData>& gnssData,
    bool enableSD = false,
    int sdMode = 15,
    const std::string& outputFile = "out/mm_matched_points.txt");

/**
 * @brief 状态机地图匹配核心流程。
 * @param coords GNSS 原始观测序列，函数会写入质量字段。
 * @param kmlFile 路网 KML 文件路径。
 * @param outputFile 匹配结果 TXT 输出路径。
 */
void RunStateMachineMapMatching(
    std::vector<CoorData>& coords,
    const std::string& kmlFile,
    const std::string& outputFile);

// ======================== HMM.cpp：隐马尔可夫地图匹配 ========================
/**
 * @brief HMM 地图匹配主接口。
 * @param road 路网数据；sourcePath 指向实际供算法解析的 KML 文件。
 * @param gnssData GNSS 原始观测序列。
 * @param enableSD 是否启用 SD 质量辅助。
 * @param outputFile 匹配结果 TXT 输出路径。
 * @return 匹配后的车辆状态点序列。
 */
std::vector<CarPoint> HiddenMarkovMapMatching(
    const RoadNetwork& road,
    const std::vector<CoorData>& gnssData,
    bool enableSD = false,
    const std::string& outputFile = "out/hmm_matched_points.txt");

/**
 * @brief 经典 HMM/Viterbi 地图匹配核心流程。
 * @param coords GNSS 原始观测序列；启用 SD 时会读取质量字段。
 * @param kmlFile 路网 KML 文件路径。
 * @param outputFile 匹配结果 TXT 输出路径。
 * @param enableSD 是否使用 SD 写入的质量信息调整发射/转移代价。
 */
void RunClassicHmmMapMatching(
    std::vector<CoorData>& coords,
    const std::string& kmlFile,
    const std::string& outputFile,
    bool enableSD = false);

// ======================== EKF.cpp：融合模块 ========================
/**
 * @brief EKF 融合主入口。
 * @param gnssData GNSS 原始观测序列，作为主观测来源。
 * @param matchedPoints MM 或 HMM 输出的车辆状态点序列，作为道路伪观测来源。
 * @param options EKF 参数配置。
 * @return 融合后的车辆状态点序列，状态量为 [x, y, velE, velN]。
 */
std::vector<CarPoint> RunEKFFusion(
    const std::vector<CoorData>& gnssData,
    const std::vector<CarPoint>& matchedPoints,
    const EKFOptions& options = EKFOptions());
