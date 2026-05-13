#include "Head.h"

#include <clocale>
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>

namespace
{
const std::string DATA_DIR = "data\\";
const std::string OUT_DIR = "out\\";

const std::string DEFAULT_COOR_FILE = DATA_DIR + "p500068.23coor";
const std::string DEFAULT_SD_COOR_FILE = DATA_DIR + "vivo068.23coor";
const std::string DEFAULT_TRUTH_FILE = DATA_DIR + "20230309.txt";
const std::string DEFAULT_ROAD_FILE = DATA_DIR + "road.kml";

void ConfigureConsoleEncoding()
{
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    std::setlocale(LC_ALL, ".UTF-8");
}

int ReadSDSwitch()
{
    int enable = 0;
    std::cout << "是否启用 SD 辅助？(0: 关闭, 1: 开启): ";
    std::cin >> enable;
    return enable == 0 ? 0 : 1;
}

int ReadEKFSDNoiseSwitch()
{
    int enable = 1;
    std::cout << "是否使用 SD 质量信息调整 EKF 观测噪声？(0: 不使用, 1: 使用): ";
    std::cin >> enable;
    return enable == 0 ? 0 : 1;
}

int ReadEKFSource()
{
    int source = 1;
    std::cout << "选择 EKF 道路伪观测来源（1: MM结果, 2: HMM结果）: ";
    std::cin >> source;
    return source == 2 ? 2 : 1;
}

int ReadEKFSDMatchedSwitch()
{
    int enable = 0;
    std::cout << "EKF 道路伪观测是否使用 SD 后地图匹配结果？(0: 普通MM/HMM, 1: SD后MM/HMM): ";
    std::cin >> enable;
    return enable == 0 ? 0 : 1;
}

SDMode ReadSDMode()
{
    int mode = 2;
    std::cout << "选择 SD 模式（1: 3sigma, 2: 前一时刻速度位置预测）: ";
    std::cin >> mode;

    if (mode == 1) {
        return SDMode::ThreeSigma;
    }
    if (mode == 2) {
        return SDMode::SpeedPositionPrediction;
    }

    std::cout << "SD 模式输入无效，默认使用前一时刻速度位置预测。" << std::endl;
    return SDMode::SpeedPositionPrediction;
}

void ApplySDQualityToCoords(std::vector<CoorData>& coords)
{
    std::vector<CarPoint> sdPoints = ConvertCoorToCarPoints(coords);
    RunSpeedDetection(sdPoints, SDMode::SpeedPositionPrediction);

    const size_t n = std::min(coords.size(), sdPoints.size());
    for (size_t i = 0; i < n; ++i) {
        coords[i].trusted = sdPoints[i].trusted;
        coords[i].trust_reason = sdPoints[i].trust_reason;
        coords[i].quality_score = sdPoints[i].quality_score;
    }
}

std::vector<CarPoint> ReadEKFMatchedPoints(int source, bool useSDMatched)
{
    const bool useHMM = (source == 2);
    const std::string prefix = useHMM
        ? (useSDMatched ? "hmm_sd" : "hmm")
        : (useSDMatched ? "mm_sd" : "mm");
    const std::string carFile = OUT_DIR + prefix + "_car_points.txt";
    const std::string matchedFile = OUT_DIR + prefix + "_matched_points.txt";

    std::vector<CarPoint> points = ReadCarTrajectory(carFile);
    if (points.empty()) {
        std::cout << "未读取到 " << carFile
            << "，尝试读取 " << matchedFile << std::endl;
        points = ReadCarTrajectory(matchedFile);
    }
    return points;
}

std::string BuildEKFOutputFile(int source, bool useSDMatched, bool useSDNoise)
{
    const std::string sourceName = (source == 2) ? "hmm" : "mm";
    const std::string matchTag = useSDMatched ? "sdmatch" : "normalmatch";
    const std::string noiseTag = useSDNoise ? "sdnoise" : "fixednoise";

    return OUT_DIR + "ekf_" + sourceName + "_" + matchTag + "_" + noiseTag + "_points.txt";
}

/**
 * 函数功能：
 *     主流程公共准备步骤：读取原始路网，并同步生成虚拟车道文件。
 *
 * 输入：
 *     roadFile：原始 KML 路网文件路径。
 *     road：输出参数，保存原始路网。
 *     virtualRoad：输出参数，保存虚拟车道路网。
 *
 * 输出：
 *     bool：路网读取和虚拟车道文件写出都成功时返回 true。
 *
 * 调用位置：
 *     main.cpp 的模式 2 和模式 3。
 *
 * 数据流作用：
 *     data\road.kml -> RoadNetwork / virtual_lane_lines.kml / virtual_lane_lines.txt。
 */
bool PrepareRoadWithVirtualLane(
    const std::string& roadFile,
    RoadNetwork& road,
    RoadNetwork& virtualRoad)
{
    road = ReadKmlRoadNetwork(roadFile);
    if (road.polylines.empty()) {
        std::cerr << "错误：路网为空，无法继续。" << std::endl;
        return false;
    }

    virtualRoad = GenerateVirtualLaneNetwork(road, true);
    virtualRoad.sourcePath = OUT_DIR + "virtual_lane_lines.kml";

    bool okKml = WriteRoadNetworkToKml(virtualRoad.sourcePath, virtualRoad);
    bool okTxt = WriteRoadNetworkToTxt(OUT_DIR + "virtual_lane_lines.txt", virtualRoad);
    if (!okKml || !okTxt) {
        std::cerr << "错误：虚拟车道文件生成失败。" << std::endl;
        return false;
    }

    return true;
}
}

int main()
{
    ConfigureConsoleEncoding();

    std::cout << "选择模式（1: SD效果验证, 2: 状态机MM, 3: HMM, 4: EKF融合）: ";

    int mode = 0;
    std::cin >> mode;

    if (mode == 1)
    {
        std::vector<CoorData> coords = ReadCoorFile(DEFAULT_SD_COOR_FILE);
        std::vector<TxtData> truthdata = ReadTxtFile(DEFAULT_TRUTH_FILE);
        SDMode sdMode = ReadSDMode();

        EvaluateSpeedDetection(
            coords,
            truthdata,
            sdMode,
            5.0,
            OUT_DIR + "sd_predict_recall_debug.txt");

        return 0;
    }

    if (mode == 2)
    {
        int enableSD = ReadSDSwitch();
        std::vector<CoorData> coordsRaw = ReadCoorFile(DEFAULT_COOR_FILE);
        std::vector<CoorData> coords = FilterCoordsByExcludeWindows(coordsRaw);

        if (coords.empty()) {
            std::cerr << "错误：轨迹为空，无法继续。" << std::endl;
            return -1;
        }

        RoadNetwork road;
        RoadNetwork virtualRoad;
        if (!PrepareRoadWithVirtualLane(DEFAULT_ROAD_FILE, road, virtualRoad)) {
            return -1;
        }

        const std::string matchedOutputFile = enableSD != 0
            ? OUT_DIR + "mm_sd_matched_points.txt"
            : OUT_DIR + "mm_matched_points.txt";
        const std::string carOutputFile = enableSD != 0
            ? OUT_DIR + "mm_sd_car_points.txt"
            : OUT_DIR + "mm_car_points.txt";

        std::vector<CarPoint> result = MapMarching(
            road,
            coords,
            enableSD != 0,
            15,
            matchedOutputFile);

        WriteCarTrajectory(carOutputFile, result);
        std::cout << "状态机 MM 完成，结果已输出到 "
            << matchedOutputFile << " 和 " << carOutputFile << std::endl;
        return 0;
    }

    if (mode == 3)
    {
        int enableSD = ReadSDSwitch();
        std::vector<CoorData> coordsRaw = ReadCoorFile(DEFAULT_COOR_FILE);
        std::vector<CoorData> coords = FilterCoordsByExcludeWindows(coordsRaw);

        if (coords.empty()) {
            std::cerr << "错误：轨迹为空，无法继续。" << std::endl;
            return -1;
        }

        RoadNetwork road;
        RoadNetwork virtualRoad;
        if (!PrepareRoadWithVirtualLane(DEFAULT_ROAD_FILE, road, virtualRoad)) {
            return -1;
        }

        const std::string matchedOutputFile = enableSD != 0
            ? OUT_DIR + "hmm_sd_matched_points.txt"
            : OUT_DIR + "hmm_matched_points.txt";
        const std::string carOutputFile = enableSD != 0
            ? OUT_DIR + "hmm_sd_car_points.txt"
            : OUT_DIR + "hmm_car_points.txt";

        std::vector<CarPoint> result = HiddenMarkovMapMatching(
            virtualRoad,
            coords,
            enableSD != 0,
            matchedOutputFile);

        WriteCarTrajectory(carOutputFile, result);
        std::cout << "HMM 完成，结果已输出到 "
            << matchedOutputFile << " 和 " << carOutputFile << std::endl;
        return 0;
    }

    if (mode == 4)
    {
        int source = ReadEKFSource();
        int useSDMatched = ReadEKFSDMatchedSwitch();
        int enableSD = ReadEKFSDNoiseSwitch();

        std::vector<CoorData> coordsRaw = ReadCoorFile(DEFAULT_COOR_FILE);
        std::vector<CoorData> coords = FilterCoordsByExcludeWindows(coordsRaw);
        if (coords.empty()) {
            std::cerr << "错误：轨迹为空，无法继续。" << std::endl;
            return -1;
        }

        if (enableSD != 0) {
            ApplySDQualityToCoords(coords);
        }

        std::vector<CarPoint> matchedPoints = ReadEKFMatchedPoints(
            source,
            useSDMatched != 0);
        if (matchedPoints.empty()) {
            std::cerr << "错误：EKF 道路伪观测为空，请先运行对应的 MM 或 HMM 模式。" << std::endl;
            return -1;
        }

        EKFOptions options;
        std::vector<CarPoint> result = RunEKFFusion(coords, matchedPoints, options);
        const std::string ekfOutputFile = BuildEKFOutputFile(
            source,
            useSDMatched != 0,
            enableSD != 0);
        if (!WriteCarTrajectory(ekfOutputFile, result)) {
            return -1;
        }

        std::cout << "EKF 融合完成，结果已输出到 "
            << ekfOutputFile << std::endl;
        return 0;
    }

    std::cerr << "输入有误。" << std::endl;
    return -1;
}

