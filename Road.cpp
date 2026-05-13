#include "Head.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <regex>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;

struct TrailPointLL
{
    double lat = 0.0;
    double lon = 0.0;
    double h = 0.0;
};

struct TrailPointXY
{
    double x = 0.0;
    double y = 0.0;
};

static const double TRAIL_R_EARTH = 6378137.0;

// =====================================================
// 虚拟车道参数
// =====================================================
// 5 条虚拟线：左2.4、左1.2、中心、右1.2、右2.4
static const int TRAIL_LANE_OFFSET_COUNT = 5;

static const double TRAIL_LANE_OFFSETS[TRAIL_LANE_OFFSET_COUNT] = {
    -2.4,
    -1.2,
     0.0,
     1.2,
     2.4
};

/**
 * 函数功能：
 *     根据开关返回候选车道偏移数量。
 *
 * 输入：
 *     useVirtualLane：是否启用虚拟车道。
 *
 * 输出：
 *     int：启用时返回 5 条偏移线数量；关闭时返回 1。
 *
 * 调用位置：
 *     MM.cpp 的候选生成阶段根据该数量遍历道路中心线/虚拟车道偏移。
 *
 * 数据流作用：
 *     虚拟车道开关 -> MM 候选数量控制。
 */
int GetVirtualLaneOffsetCount(bool useVirtualLane)
{
    return useVirtualLane ? TRAIL_LANE_OFFSET_COUNT : 1;
}

/**
 * 函数功能：
 *     返回指定虚拟车道的横向偏移量。
 *
 * 输入：
 *     useVirtualLane：是否启用虚拟车道。
 *     laneIndex：车道偏移索引。
 *
 * 输出：
 *     double：相对道路中心线的横向偏移，单位 m。
 *
 * 调用位置：
 *     1. MM.cpp 为同一条道路生成不同横向偏移候选。
 *     2. Road.cpp 生成虚拟车道 KML/TXT 时计算每条线的偏移。
 *
 * 数据流作用：
 *     laneIndex -> offset_m，用于虚拟车道几何生成和 MM 候选点构建。
 */
double GetVirtualLaneOffset(bool useVirtualLane, int laneIndex)
{
    if (!useVirtualLane)
        return 0.0;

    if (laneIndex < 0)
        laneIndex = 0;

    if (laneIndex >= TRAIL_LANE_OFFSET_COUNT)
        laneIndex = TRAIL_LANE_OFFSET_COUNT - 1;

    return TRAIL_LANE_OFFSETS[laneIndex];
}

static string trailMakeSiblingPath(
    const string& baseFile,
    const string& newName)
{
    size_t p = baseFile.find_last_of("/\\");

    if (p == string::npos)
        return newName;

    return baseFile.substr(0, p + 1) + newName;
}


// =====================================================
// 局部坐标转换
// =====================================================

static void trailMakeLocalRef(
    const vector<TrailPointLL>& pts,
    LocalRef& ref)
{
    double sumLat = 0.0;
    double sumLon = 0.0;
    int count = 0;

    for (size_t i = 0; i < pts.size(); ++i)
    {
        sumLat += pts[i].lat;
        sumLon += pts[i].lon;
        count++;
    }

    if (count == 0)
    {
        ref.lat0_rad = 0.0;
        ref.lon0_rad = 0.0;
        ref.cosLat0 = 1.0;
        return;
    }

    double lat0 = sumLat / (double)count;
    double lon0 = sumLon / (double)count;

    ref.lat0_rad = lat0 * PI / 180.0;
    ref.lon0_rad = lon0 * PI / 180.0;
    ref.cosLat0 = cos(ref.lat0_rad);
}

static void trailLL2XY(
    double lat_deg,
    double lon_deg,
    const LocalRef& ref,
    double& x,
    double& y)
{
    double lat_rad = lat_deg * PI / 180.0;
    double lon_rad = lon_deg * PI / 180.0;

    y = TRAIL_R_EARTH * (lat_rad - ref.lat0_rad);
    x = TRAIL_R_EARTH * (lon_rad - ref.lon0_rad) * ref.cosLat0;
}

static void trailXY2LL(
    double x,
    double y,
    const LocalRef& ref,
    double& lat_deg,
    double& lon_deg)
{
    double lat_rad = y / TRAIL_R_EARTH + ref.lat0_rad;
    double lon_rad = x / (TRAIL_R_EARTH * ref.cosLat0) + ref.lon0_rad;

    lat_deg = lat_rad * 180.0 / PI;
    lon_deg = lon_rad * 180.0 / PI;
}


// =====================================================
// KML 解析
// =====================================================

static vector<vector<TrailPointLL> > trailParseKMLLines(
    const string& inputKmlFile)
{
    vector<vector<TrailPointLL> > lines;

    ifstream fin(inputKmlFile.c_str());
    if (!fin.is_open())
    {
        cerr << "无法打开输入 KML: " << inputKmlFile << endl;
        return lines;
    }

    string content(
        (istreambuf_iterator<char>(fin)),
        istreambuf_iterator<char>());

    fin.close();

    regex re("<coordinates>([^<]+)</coordinates>");
    smatch m;
    string::const_iterator it = content.cbegin();

    while (regex_search(it, content.cend(), m, re))
    {
        string coordText = m[1].str();
        stringstream ss(coordText);
        string token;

        vector<TrailPointLL> line;

        while (ss >> token)
        {
            size_t p1 = token.find(',');
            if (p1 == string::npos)
                continue;

            size_t p2 = token.find(',', p1 + 1);

            double lon = 0.0;
            double lat = 0.0;
            double h = 0.0;

            try
            {
                if (p2 == string::npos)
                {
                    lon = stod(token.substr(0, p1));
                    lat = stod(token.substr(p1 + 1));
                    h = 0.0;
                }
                else
                {
                    lon = stod(token.substr(0, p1));
                    lat = stod(token.substr(p1 + 1, p2 - p1 - 1));
                    h = stod(token.substr(p2 + 1));
                }
            }
            catch (...)
            {
                continue;
            }

            TrailPointLL pt;
            pt.lat = lat;
            pt.lon = lon;
            pt.h = h;

            line.push_back(pt);
        }

        if (line.size() >= 2)
            lines.push_back(line);

        it = m.suffix().first;
    }

    return lines;
}


// =====================================================
// 虚拟车道生成
// =====================================================

static vector<TrailPointLL> trailBuildOffsetLine(
    const vector<TrailPointLL>& src,
    double offset)
{
    vector<TrailPointLL> out;

    if (src.size() < 2)
        return out;

    LocalRef ref;
    trailMakeLocalRef(src, ref);

    vector<TrailPointXY> xy(src.size());

    for (size_t i = 0; i < src.size(); ++i)
    {
        trailLL2XY(src[i].lat, src[i].lon, ref, xy[i].x, xy[i].y);
    }

    for (size_t i = 0; i < xy.size(); ++i)
    {
        double tx = 0.0;
        double ty = 0.0;

        if (i == 0)
        {
            tx = xy[1].x - xy[0].x;
            ty = xy[1].y - xy[0].y;
        }
        else if (i + 1 == xy.size())
        {
            tx = xy[i].x - xy[i - 1].x;
            ty = xy[i].y - xy[i - 1].y;
        }
        else
        {
            tx = xy[i + 1].x - xy[i - 1].x;
            ty = xy[i + 1].y - xy[i - 1].y;
        }

        double len = sqrt(tx * tx + ty * ty);
        if (len < 1e-8)
        {
            TrailPointLL p = src[i];
            out.push_back(p);
            continue;
        }

        double nx = -ty / len;
        double ny = tx / len;

        double ox = xy[i].x + offset * nx;
        double oy = xy[i].y + offset * ny;

        TrailPointLL p;
        trailXY2LL(ox, oy, ref, p.lat, p.lon);
        p.h = src[i].h;

        out.push_back(p);
    }

    return out;
}

static bool trailWriteVirtualLaneKML(
    const string& outputKmlFile,
    const vector<vector<vector<TrailPointLL> > >& allVirtualLines)
{
    ofstream fout(outputKmlFile.c_str());
    if (!fout.is_open())
    {
        cerr << "无法写入虚拟车道 KML: " << outputKmlFile << endl;
        return false;
    }

    fout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    fout << "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n";
    fout << "<Document>\n";
    fout << "<name>Virtual Lane Lines</name>\n";

    fout << "<Style id=\"laneLine\"><LineStyle><color>ff00ffff</color><width>2</width></LineStyle></Style>\n";

    fout << fixed << setprecision(10);

    for (size_t roadIdx = 0; roadIdx < allVirtualLines.size(); ++roadIdx)
    {
        for (size_t laneIdx = 0; laneIdx < allVirtualLines[roadIdx].size(); ++laneIdx)
        {
            const vector<TrailPointLL>& line = allVirtualLines[roadIdx][laneIdx];

            fout << "<Placemark>\n";
            fout << "<name>road_" << roadIdx << "_lane_" << laneIdx << "</name>\n";
            fout << "<styleUrl>#laneLine</styleUrl>\n";
            fout << "<LineString><tessellate>1</tessellate><coordinates>\n";

            for (size_t i = 0; i < line.size(); ++i)
            {
                fout << line[i].lon << "," << line[i].lat << "," << line[i].h << "\n";
            }

            fout << "</coordinates></LineString>\n";
            fout << "</Placemark>\n";
        }
    }

    fout << "</Document>\n";
    fout << "</kml>\n";

    fout.close();
    return true;
}

static bool trailWriteVirtualLaneTXT(
    const string& outputTxtFile,
    const vector<vector<vector<TrailPointLL> > >& allVirtualLines)
{
    ofstream fout(outputTxtFile.c_str());
    if (!fout.is_open())
    {
        cerr << "无法写入虚拟车道 TXT: " << outputTxtFile << endl;
        return false;
    }

    fout << fixed << setprecision(10);

    fout << "road_id lane_id point_id lat lon h offset_m\n";

    for (size_t roadIdx = 0; roadIdx < allVirtualLines.size(); ++roadIdx)
    {
        for (size_t laneIdx = 0; laneIdx < allVirtualLines[roadIdx].size(); ++laneIdx)
        {
            const vector<TrailPointLL>& line = allVirtualLines[roadIdx][laneIdx];

            double offset = GetVirtualLaneOffset(true, (int)laneIdx);

            for (size_t i = 0; i < line.size(); ++i)
            {
                fout << roadIdx << " "
                    << laneIdx << " "
                    << i << " "
                    << line[i].lat << " "
                    << line[i].lon << " "
                    << line[i].h << " "
                    << offset << "\n";
            }
        }
    }

    fout.close();
    return true;
}

/**
 * 函数功能：
 *     从 KML 文件读取道路中心线，并生成虚拟车道 KML/TXT 文件。
 *
 * 输入：
 *     inputKmlFile：原始路网 KML 文件路径。
 *     outputKmlFile：虚拟车道 KML 输出路径。
 *     outputTxtFile：虚拟车道 TXT 输出路径。
 *
 * 输出：
 *     bool：两个输出文件都成功写出时返回 true。
 *
 * 调用位置：
 *     可作为文件级接口独立生成虚拟车道；当前 main.cpp 更常用 GenerateVirtualLaneNetwork 后再写文件。
 *
 * 数据流作用：
 *     road.kml -> virtual_lane_lines.kml / virtual_lane_lines.txt。
 */
bool GenerateVirtualLaneFiles(
    const string& inputKmlFile,
    const string& outputKmlFile,
    const string& outputTxtFile)
{
    vector<vector<TrailPointLL> > roadLines = trailParseKMLLines(inputKmlFile);

    if (roadLines.empty())
    {
        cerr << "未从 KML 中解析到道路中心线，无法生成虚拟车道" << endl;
        return false;
    }

    vector<vector<vector<TrailPointLL> > > allVirtualLines;
    allVirtualLines.resize(roadLines.size());

    for (size_t roadIdx = 0; roadIdx < roadLines.size(); ++roadIdx)
    {
        for (int laneIdx = 0; laneIdx < TRAIL_LANE_OFFSET_COUNT; ++laneIdx)
        {
            double offset = GetVirtualLaneOffset(true, laneIdx);
            vector<TrailPointLL> laneLine =
                trailBuildOffsetLine(roadLines[roadIdx], offset);

            if (!laneLine.empty())
                allVirtualLines[roadIdx].push_back(laneLine);
        }
    }

    bool okKml = trailWriteVirtualLaneKML(outputKmlFile, allVirtualLines);
    bool okTxt = trailWriteVirtualLaneTXT(outputTxtFile, allVirtualLines);

    if (okKml && okTxt)
    {
        cout << "虚拟车道 KML 输出完成: " << outputKmlFile << endl;
        cout << "虚拟车道 TXT 输出完成: " << outputTxtFile << endl;
    }

    return okKml && okTxt;
}

/**
 * 函数功能：
 *     从 RoadNetwork 直接生成虚拟车道路网，输出结构与输入结构保持一致。
 *
 * 输入：
 *     inputRoad：原始道路中心线路网。
 *     useVirtualLane：是否生成虚拟车道。
 *
 * 输出：
 *     RoadNetwork：包含虚拟车道折线的路网；hasVirtualLanes 标记为 true。
 *
 * 调用位置：
 *     main.cpp 的 PrepareRoadWithVirtualLane 在读取 road.kml 后调用。
 *
 * 数据流作用：
 *     原始 RoadNetwork -> 虚拟车道 RoadNetwork -> File.cpp 写出 KML/TXT，HMM 可直接使用虚拟车道 KML。
 */
RoadNetwork GenerateVirtualLaneNetwork(const RoadNetwork& inputRoad, bool useVirtualLane)
{
    RoadNetwork output;
    output.sourcePath = inputRoad.sourcePath;
    output.hasVirtualLanes = useVirtualLane;

    if (!useVirtualLane) {
        output = inputRoad;
        output.hasVirtualLanes = false;
        return output;
    }

    int outId = 0;
    for (const auto& srcLine : inputRoad.polylines) {
        if (srcLine.points.size() < 2) {
            continue;
        }

        vector<TrailPointLL> src;
        src.reserve(srcLine.points.size());
        for (const auto& node : srcLine.points) {
            TrailPointLL p;
            p.lat = node.BLH[0];
            p.lon = node.BLH[1];
            p.h = node.BLH[2];
            src.push_back(p);
        }

        for (int laneIdx = 0; laneIdx < TRAIL_LANE_OFFSET_COUNT; ++laneIdx) {
            double offset = GetVirtualLaneOffset(true, laneIdx);
            vector<TrailPointLL> lane = trailBuildOffsetLine(src, offset);
            if (lane.empty()) {
                continue;
            }

            RoadPolyline outLine;
            outLine.id = outId++;
            outLine.laneId = laneIdx;
            outLine.laneOffset = offset;
            outLine.name = "road_" + std::to_string(srcLine.id)
                + "_lane_" + std::to_string(laneIdx);
            outLine.points.reserve(lane.size());

            for (const auto& p : lane) {
                RoadNode node;
                node.BLH[0] = p.lat;
                node.BLH[1] = p.lon;
                node.BLH[2] = p.h;
                outLine.points.push_back(node);
            }

            output.polylines.push_back(outLine);
        }
    }

    std::cout << "虚拟车道生成完成，线数量: "
        << output.polylines.size() << std::endl;

    return output;
}


