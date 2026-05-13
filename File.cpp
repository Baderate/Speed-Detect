#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <regex>

#include "Head.h"

namespace
{
LocalRef BuildLocalRefFromCarPoints(const std::vector<CarPoint>& points)
{
    double sumLat = 0.0;
    double sumLon = 0.0;
    int count = 0;

    for (const auto& point : points) {
        if (std::isfinite(point.BLH[0]) && std::isfinite(point.BLH[1])) {
            sumLat += point.BLH[0];
            sumLon += point.BLH[1];
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

void FillCarPointXY(std::vector<CarPoint>& points)
{
    LocalRef ref = BuildLocalRefFromCarPoints(points);

    for (auto& point : points) {
        if (std::isfinite(point.x) && std::isfinite(point.y)) {
            continue;
        }

        if (!std::isfinite(point.BLH[0]) || !std::isfinite(point.BLH[1])) {
            continue;
        }

        const double latRad = point.BLH[0] * PI / 180.0;
        const double lonRad = point.BLH[1] * PI / 180.0;
        point.y = earth_R * (latRad - ref.lat0_rad);
        point.x = earth_R * (lonRad - ref.lon0_rad) * ref.cosLat0;
    }
}
}

/**
 * 函数功能：
 *     读取 txt 格式参考真值文件，并转换为项目统一的 TxtData 数组。
 *
 * 输入：
 *     filename：txt 真值文件路径，例如 data\20230309.txt。
 *
 * 输出：
 *     std::vector<TxtData>：参考真值点数组。
 *     每个点包含 dataType、UTC 毫秒时间戳、BLH 坐标和 redius。
 *
 * 调用位置：
 *     1. main.cpp 的 SD 效果验证模式读取参考真值。
 *     2. SD.cpp 的 EvaluateSpeedDetection / RunSDModule_RecallOnly 用它和 GNSS 结果做时间对齐评估。
 *
 * 数据流作用：
 *     data\*.txt -> TxtData -> SD 评估/误差统计。
 */
std::vector<TxtData> ReadTxtFile(const std::string& filename)
{
    std::vector<TxtData> txtDataVec;  // 用于存放读取结果
    std::ifstream inFile(filename);

    if (!inFile) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return txtDataVec;
    }

    std::string line;
    while (std::getline(inFile, line)) {
        TxtData item;

        // 若每行外层带有大括号，则先去掉
        if (!line.empty() && line.front() == '{') {
            line = line.substr(1, line.length() - 2);
        }

        // 每行按“键:值”形式解析
        std::stringstream ss(line);
        std::string field;
        while (std::getline(ss, field, ',')) {
            // 去掉字段两端空白字符
            field = field.erase(0, field.find_first_not_of(" \t"));
            field = field.erase(field.find_last_not_of(" \t") + 1);

            size_t pos = field.find(":");
            if (pos != std::string::npos) {
                std::string key = field.substr(0, pos);
                std::string value = field.substr(pos + 1);

                // 去掉字符串值两端的单引号
                if (!value.empty() && value.front() == '\'' && value.back() == '\'') {
                    value = value.substr(1, value.length() - 2);
                }

                // 按字段名写入结构体
                if (key == "'dataType'") {
                    item.dataType = value;
                }
                else if (key == "'timestamp'") {
                    item.timestamp = std::stoll(value);
                }
                else if (key == "'y'") {
                    item.BLH[0] = std::stod(value);  // 纬度
                }
                else if (key == "'x'") {
                    item.BLH[1] = std::stod(value);  // 经度
                }
                else if (key == "'z'") {
                    item.BLH[2] = std::stod(value);  // 高度
                }
                else if (key == "'redius'") {
                    item.redius = std::stod(value);
                }
            }
        }

        txtDataVec.push_back(item);
    }

    inFile.close();
    printf("真值文件（txt）读取完毕\n");
    std::cout << "真值路径点数量:" << txtDataVec.size() << std::endl;
    return txtDataVec;
}

/**
 * 函数功能：
 *     读取 coor 格式 GNSS 原始轨迹文件，并转换为项目统一的 CoorData 数组。
 *
 * 输入：
 *     filename：coor 文件路径，例如 data\p500068.23coor。
 *
 * 输出：
 *     std::vector<CoorData>：GNSS 原始观测数组。
 *     主要字段包括 GPS 时间、BLH 坐标、ENU 速度、质量信息和保留扩展列。
 *
 * 调用位置：
 *     1. main.cpp 的 SD 效果验证模式读取 GNSS 轨迹。
 *     2. main.cpp 的状态机 MM 模式读取待匹配轨迹。
 *     3. main.cpp 的 HMM 模式读取待匹配轨迹。
 *
 * 数据流作用：
 *     data\*.coor -> CoorData -> SD / MM / HMM。
 */
std::vector<CoorData> ReadCoorFile(const std::string& filename)
{
    std::vector<CoorData> data;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return data;
    }

    std::string line;

    // 跳过表头
    if (!std::getline(file, line)) {
        std::cerr << "文件为空: " << filename << std::endl;
        return data;
    }

    int lineNo = 1;

    while (std::getline(file, line)) {
        ++lineNo;
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::vector<double> cols;
        double val = 0.0;

        while (iss >> val) {
            cols.push_back(val);
        }

        // 至少要保证前 33 列在
        if (cols.size() < 33) {
            std::cerr << "第 " << lineNo << " 行字段数不足，跳过: " << line << std::endl;
            continue;
        }

        CoorData entry;

        entry.year = static_cast<int>(cols[0]);
        entry.month = static_cast<int>(cols[1]);
        entry.day = static_cast<int>(cols[2]);
        entry.hour = static_cast<int>(cols[3]);
        entry.minute = static_cast<int>(cols[4]);
        entry.second = cols[5];

        entry.week = static_cast<int>(cols[6]);
        entry.ws = cols[7];
        entry.nSV = static_cast<int>(cols[8]);
        entry.sigma0 = cols[9];

        entry.BLH[0] = cols[10];   // B
        entry.BLH[1] = cols[11];   // L
        entry.BLH[2] = cols[12];   // H

        // ===== 速度：VU1 VN1 VE1 =====
        entry.velU = cols[13];
        entry.velN = cols[14];
        entry.velE = cols[15];

        // ===== 相位相关量：VU2 VN2 VE2 =====
        entry.phaU = cols[16];
        entry.phaN = cols[17];
        entry.phaE = cols[18];

        // 内部数组沿用 [N, E, U] 顺序，保留原始速度分量的同时方便算法读取。
        entry.V_ENU[0] = entry.velN;
        entry.V_ENU[1] = entry.velE;
        entry.V_ENU[2] = entry.velU;

        entry.sigma1 = cols[19];
        entry.sigma2 = cols[20];
        entry.sigma[0] = entry.sigma1;
        entry.sigma[1] = entry.sigma2;

        entry.LObs = static_cast<int>(cols[21]);
        entry.clk1 = cols[22];
        entry.trop = cols[23];

        for (int i = 0; i < 5; ++i) {
            entry.CBP2[i] = cols[24 + i];
        }

        for (int i = 0; i < 4; ++i) {
            entry.ISB[i] = cols[29 + i];
        }

        // 后面如果还有附加列，先保存
        if (cols.size() > 33) {
            entry.extraCols.assign(cols.begin() + 33, cols.end());
        }

        data.push_back(entry);
    }

    file.close();

    printf("轨迹文件读取完毕\n");
    std::cout << "轨迹路径点数量:" << data.size() << std::endl;
    return data;
}

/**
 * 函数功能：
 *     读取已经输出的轨迹文本文件，并直接转换为 CarPoint 数组。
 *
 * 输入：
 *     filePath：轨迹文件路径。
 *     支持两类格式：
 *     1. GPS_matched timestamp lat lon h radius
 *     2. timestamp lat lon h
 *     3. timestamp lat lon h x y vE vN vU quality trusted reason
 *
 * 输出：
 *     std::vector<CarPoint>：统一车辆状态点数组。
 *     若文件中没有 x/y，则函数内部自动根据 BLH 补齐局部平面坐标。
 *
 * 调用位置：
 *     1. MM.cpp 的 MapMarching 在核心匹配写出 txt 后，直接读回 CarPoint 结果。
 *     2. HMM.cpp 的 HiddenMarkovMapMatching 在 HMM 写出 txt 后，直接读回 CarPoint 结果。
 *
 * 数据流作用：
 *     out\*_matched_points.txt 或 out\*_car_points.txt -> CarPoint。
 */
std::vector<CarPoint> ReadCarTrajectory(const std::string& filePath)
{
    std::vector<CarPoint> points;

    std::ifstream fin(filePath.c_str());
    if (!fin.is_open())
    {
        std::cerr << "无法打开轨迹文件: " << filePath << std::endl;
        return points;
    }

    std::string line;
    int lineNo = 0;

    while (std::getline(fin, line))
    {
        lineNo++;

        if (line.empty())
            continue;

        std::stringstream ss(line);

        std::string firstToken;
        ss >> firstToken;

        if (!ss)
            continue;

        if (firstToken == "timestamp_ms") {
            continue;
        }

        CarPoint point;
        point.timestamp = 0;
        point.BLH[0] = NaNValue();
        point.BLH[1] = NaNValue();
        point.BLH[2] = NaNValue();
        point.x = NaNValue();
        point.y = NaNValue();
        point.V_ENU[0] = NaNValue();
        point.V_ENU[1] = NaNValue();
        point.V_ENU[2] = NaNValue();

        // 情况 1：
        // GPS_matched timestamp lat lon h radius
        // GPS_inv     timestamp lat lon h radius
        if (firstToken == "GPS_matched" ||
            firstToken == "GPS_inv" ||
            firstToken == "GPS" ||
            firstToken == "matched")
        {
            if (!(ss >> point.timestamp
                >> point.BLH[0]
                >> point.BLH[1]
                >> point.BLH[2]))
            {
                std::cerr << "[警告] 第 " << lineNo
                    << " 行轨迹格式无效，跳过：" << line << std::endl;
                continue;
            }
        }
        else
        {
            // 情况 2：
            // timestamp lat lon h
            try
            {
                point.timestamp = std::stoll(firstToken);
            }
            catch (...)
            {
                std::cerr << "[警告] 第 " << lineNo
                    << " 行轨迹格式无效，跳过：" << line << std::endl;
                continue;
            }

            if (!(ss >> point.BLH[0]
                >> point.BLH[1]
                >> point.BLH[2]))
            {
                std::cerr << "[警告] 第 " << lineNo
                    << " 行轨迹格式无效，跳过：" << line << std::endl;
                continue;
            }

            // 若文件中还包含 x/y、速度、质量等字段，则继续读取。
            ss >> point.x
                >> point.y
                >> point.V_ENU[0]
                >> point.V_ENU[1]
                >> point.V_ENU[2]
                >> point.quality_score;

            int trustedInt = 1;
            if (ss >> trustedInt) {
                point.trusted = (trustedInt != 0);
            }
            ss >> point.trust_reason;
        }

        points.push_back(point);
    }

    fin.close();
    FillCarPointXY(points);

    std::cout << "成功读取轨迹文件 " << filePath
        << " ：" << points.size()
        << " 个点" << std::endl;

    return points;
}

/**
 * 函数功能：
 *     读取地图匹配结果文件，并转换为 TxtData 数组。
 *
 * 输入：
 *     filename：地图匹配结果文件路径。
 *     defaultDataType：当文件行没有数据类型标签时使用的默认类型。
 *     defaultRedius：当文件行没有半径字段时使用的默认半径。
 *
 * 输出：
 *     std::vector<TxtData>：可作为“匹配后轨迹/伪真值”的点数组。
 *
 * 调用位置：
 *     目前主流程不直接调用，保留给后续误差评估、EKF 或匹配结果二次读取使用。
 *
 * 数据流作用：
 *     out\*_matched_points.txt -> TxtData -> 评估/融合模块。
 */
std::vector<TxtData> ReadRoadMatchedFile(
    const std::string& filename,
    const std::string& defaultDataType,
    double defaultRedius)
{
    std::vector<TxtData> data;

    std::ifstream fin(filename.c_str());
    if (!fin.is_open())
    {
        std::cerr << "无法打开地图匹配结果文件: " << filename << std::endl;
        return data;
    }

    std::string line;
    int lineNo = 0;

    while (std::getline(fin, line))
    {
        lineNo++;

        if (line.empty())
            continue;

        std::stringstream ss(line);

        TxtData d;
        d.dataType = defaultDataType;
        d.timestamp = 0;
        d.BLH[0] = 0.0;
        d.BLH[1] = 0.0;
        d.BLH[2] = 0.0;
        d.redius = defaultRedius;

        std::string firstToken;
        ss >> firstToken;

        if (!ss)
            continue;

        // 情况 1：
        // GPS_matched timestamp lat lon h radius
        if (firstToken == "GPS_matched" ||
            firstToken == "GPS_inv" ||
            firstToken == "GPS" ||
            firstToken == "matched")
        {
            d.dataType = firstToken;

            if (!(ss >> d.timestamp
                >> d.BLH[0]
                >> d.BLH[1]
                >> d.BLH[2]
                >> d.redius))
            {
                std::cerr << "[警告] 第 " << lineNo
                    << " 行数据格式无效，跳过：" << line << std::endl;
                continue;
            }
        }
        else
        {
            // 情况 2：
            // timestamp lat lon h radius
            // 或 timestamp lat lon h
            try
            {
                d.timestamp = std::stoll(firstToken);
            }
            catch (...)
            {
                std::cerr << "[警告] 第 " << lineNo
                    << " 行数据格式无效，跳过：" << line << std::endl;
                continue;
            }

            if (!(ss >> d.BLH[0] >> d.BLH[1] >> d.BLH[2]))
            {
                std::cerr << "[警告] 第 " << lineNo
                    << " 行数据格式无效，跳过：" << line << std::endl;
                continue;
            }

            if (!(ss >> d.redius))
                d.redius = defaultRedius;
        }

        data.push_back(d);
    }

    fin.close();

    std::cout << "[读取完成] " << filename
        << " 共读取 " << data.size()
        << " 个有效数据点" << std::endl;

    return data;
}

/**
 * 函数功能：
 *     读取 KML 路网文件中的 <coordinates> 折线，并转换为 RoadNetwork。
 *
 * 输入：
 *     filename：KML 路网文件路径，例如 data\road.kml。
 *
 * 输出：
 *     RoadNetwork：路网结构。
 *     sourcePath 保存原 KML 路径，polylines 保存每条道路折线。
 *
 * 调用位置：
 *     main.cpp 的 PrepareRoadWithVirtualLane 先调用本函数读取原始路网。
 *
 * 数据流作用：
 *     data\road.kml -> RoadNetwork -> Road.cpp 生成虚拟车道 / MM-HMM 使用路网路径。
 */
RoadNetwork ReadKmlRoadNetwork(const std::string& filename)
{
    RoadNetwork road;
    road.sourcePath = filename;

    std::ifstream fin(filename.c_str());
    if (!fin.is_open()) {
        std::cerr << "无法打开 KML 路网文件: " << filename << std::endl;
        return road;
    }

    std::string content(
        (std::istreambuf_iterator<char>(fin)),
        std::istreambuf_iterator<char>());
    fin.close();

    std::regex coordRegex("<coordinates>([^<]+)</coordinates>");
    std::smatch match;
    std::string::const_iterator it = content.cbegin();

    int roadId = 0;
    while (std::regex_search(it, content.cend(), match, coordRegex)) {
        std::stringstream ss(match[1].str());
        std::string token;

        RoadPolyline line;
        line.id = roadId;
        line.name = "road_" + std::to_string(roadId);

        while (ss >> token) {
            size_t p1 = token.find(',');
            if (p1 == std::string::npos) {
                continue;
            }

            size_t p2 = token.find(',', p1 + 1);
            RoadNode node;

            try {
                node.BLH[1] = std::stod(token.substr(0, p1));
                if (p2 == std::string::npos) {
                    node.BLH[0] = std::stod(token.substr(p1 + 1));
                    node.BLH[2] = 0.0;
                }
                else {
                    node.BLH[0] = std::stod(token.substr(p1 + 1, p2 - p1 - 1));
                    node.BLH[2] = std::stod(token.substr(p2 + 1));
                }
            }
            catch (...) {
                continue;
            }

            line.points.push_back(node);
        }

        if (line.points.size() >= 2) {
            road.polylines.push_back(line);
            ++roadId;
        }

        it = match.suffix().first;
    }

    std::cout << "KML 路网读取完毕，中心线数量: "
        << road.polylines.size() << std::endl;
    return road;
}

/**
 * 函数功能：
 *     将统一车辆状态点 CarPoint 写出为文本轨迹文件。
 *
 * 输入：
 *     filename：输出文件路径。
 *     points：车辆状态点数组。
 *
 * 输出：
 *     bool：写入成功返回 true，失败返回 false。
 *
 * 调用位置：
 *     1. main.cpp 的 SD 模式输出 sd_checked_points.txt。
 *     2. main.cpp 的 MM 模式输出 mm_car_points.txt。
 *     3. main.cpp 的 HMM 模式输出 hmm_car_points.txt。
 *
 * 数据流作用：
 *     CarPoint -> out\*_car_points.txt，方便检查每个点的坐标、速度和质量评分。
 */
bool WriteCarTrajectory(const std::string& filename, const std::vector<CarPoint>& points)
{
    std::ofstream fout(filename.c_str());
    if (!fout.is_open()) {
        std::cerr << "无法写入车辆轨迹文件: " << filename << std::endl;
        return false;
    }

    fout << "timestamp_ms lat lon h x y vE vN vU quality trusted reason\n";
    fout << std::fixed << std::setprecision(10);

    for (const auto& p : points) {
        fout << p.timestamp << " "
            << p.BLH[0] << " "
            << p.BLH[1] << " "
            << p.BLH[2] << " "
            << p.x << " "
            << p.y << " "
            << p.V_ENU[0] << " "
            << p.V_ENU[1] << " "
            << p.V_ENU[2] << " "
            << p.quality_score << " "
            << (p.trusted ? 1 : 0) << " "
            << p.trust_reason << "\n";
    }

    return true;
}

/**
 * 函数功能：
 *     将 RoadNetwork 写成 KML 折线文件。
 *
 * 输入：
 *     filename：输出 KML 文件路径。
 *     road：待输出的路网，可为原始路网或虚拟车道路网。
 *
 * 输出：
 *     bool：写入成功返回 true，失败返回 false。
 *
 * 调用位置：
 *     main.cpp 的 PrepareRoadWithVirtualLane 写出 out\virtual_lane_lines.kml。
 *
 * 数据流作用：
 *     RoadNetwork -> out\virtual_lane_lines.kml，用于地图显示和 HMM 虚拟车道输入。
 */
bool WriteRoadNetworkToKml(const std::string& filename, const RoadNetwork& road)
{
    std::ofstream fout(filename.c_str());
    if (!fout.is_open()) {
        std::cerr << "无法写入路网 KML 文件: " << filename << std::endl;
        return false;
    }

    fout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    fout << "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n";
    fout << "<Document>\n";
    fout << "<name>" << (road.hasVirtualLanes ? "Virtual Lane Lines" : "Road Network") << "</name>\n";
    fout << "<Style id=\"roadLine\"><LineStyle><color>ff00ffff</color><width>2</width></LineStyle></Style>\n";
    fout << std::fixed << std::setprecision(10);

    for (const auto& line : road.polylines) {
        if (line.points.size() < 2) {
            continue;
        }

        fout << "<Placemark>\n";
        fout << "<name>" << line.name << "</name>\n";
        fout << "<styleUrl>#roadLine</styleUrl>\n";
        fout << "<LineString><tessellate>1</tessellate><coordinates>\n";

        for (const auto& p : line.points) {
            fout << p.BLH[1] << "," << p.BLH[0] << "," << p.BLH[2] << "\n";
        }

        fout << "</coordinates></LineString>\n";
        fout << "</Placemark>\n";
    }

    fout << "</Document>\n";
    fout << "</kml>\n";
    return true;
}

/**
 * 函数功能：
 *     将 RoadNetwork 写成 TXT 明细文件。
 *
 * 输入：
 *     filename：输出 TXT 文件路径。
 *     road：待输出的路网，可为原始路网或虚拟车道路网。
 *
 * 输出：
 *     bool：写入成功返回 true，失败返回 false。
 *
 * 调用位置：
 *     main.cpp 的 PrepareRoadWithVirtualLane 写出 out\virtual_lane_lines.txt。
 *
 * 数据流作用：
 *     RoadNetwork -> out\virtual_lane_lines.txt，用于检查 road_id、lane_id、point_id 和 offset_m。
 */
bool WriteRoadNetworkToTxt(const std::string& filename, const RoadNetwork& road)
{
    std::ofstream fout(filename.c_str());
    if (!fout.is_open()) {
        std::cerr << "无法写入路网 TXT 文件: " << filename << std::endl;
        return false;
    }

    fout << "road_id lane_id point_id lat lon h offset_m\n";
    fout << std::fixed << std::setprecision(10);

    for (const auto& line : road.polylines) {
        for (size_t i = 0; i < line.points.size(); ++i) {
            const RoadNode& p = line.points[i];
            fout << line.id << " "
                << line.laneId << " "
                << i << " "
                << p.BLH[0] << " "
                << p.BLH[1] << " "
                << p.BLH[2] << " "
                << line.laneOffset << "\n";
        }
    }

    return true;
}

