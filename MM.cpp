#include "Head.h"
#include <fstream>
#include <regex>
#include <sstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <queue>
#include <map>
#include <set>
#include <string>
#include <deque>
#include <unordered_map>
#include <functional>
#include <iterator>

using namespace std;


/*
 * 纯 MM 增强版（状态机 + Viterbi 序列优化）
 *
 * 设计思想：
 * 1. 保留状态机、候选收集、虚拟车道、速度辅助；
 * 2. 不再逐帧直接输出当前 best；
 * 3. 先为每一帧生成候选集，再按分链方式做整段 Viterbi；
 * 4. 最终输出最优路径。
 *
 */

static const double R_EARTH = 6378137.0;

// ---------------------------- 状态机 ----------------------------

enum MMState
{
    MM_INIT = 0,
    MM_ON_ROAD = 1,
    MM_AT_INTERSECTION = 2,
    MM_RECOVERY = 3
};

// ---------------------------- 参数 ----------------------------

// 候选搜索
static const double MAX_SNAP_DIST = 50.0;
static const double INIT_SEARCH_RADIUS = 60.0;
static const double ONROAD_SEARCH_RADIUS = 70.0;
static const double INTERSECTION_SEARCH_RADIUS = 85.0;
static const double STOP_SEARCH_RADIUS = 18.0;

// 路口判定
static const double INTERSECTION_NEAR_DIST = 30.0;

// 初始化与恢复
static const int BAD_LIMIT = 3;

// 每帧候选保留数
static const int TOPK_PER_EPOCH = 8;

// 发射代价权重：正常道路
static const double W_OBS_ONROAD = 1.0;
static const double W_DR_ONROAD = 0.35;
static const double W_CONT_ONROAD = 0.12;
static const double W_DIR_ONROAD = 40.0;
static const double W_VEL_DIR_ONROAD = 18.0;

// 发射代价权重：路口
static const double W_OBS_INT = 0.85;
static const double W_DR_INT = 0.12;
static const double W_CONT_INT = 0.10;
static const double W_DIR_INT = 28.0;
static const double W_VEL_DIR_INT = 55.0;

// 发射代价权重：初始化 / 恢复
static const double W_OBS_INIT = 1.0;
static const double W_DR_INIT = 0.10;
static const double W_CONT_INIT = 0.04;
static const double W_DIR_INIT = 18.0;
static const double W_VEL_DIR_INIT = 10.0;

// 转移代价参数
static const double TRANS_ROUTE_BETA = 15.0;
static const double TRANS_PRED_BETA = 12.0;
static const double TRANS_LANE_W = 0.40;
static const double PEN_NON_ADJ_TRANS = 15.0;
static const double PEN_CLUSTER_TRANS = 8.0;
static const double PEN_REVERSE_TRANS = 10.0;

// 坏帧判定
static const double BAD_OBS_DIST = 22.0;
static const double BAD_DIR_DIFF = 65.0 * PI / 180.0;
static const double BAD_COST_TH = 900.0;

// 静止保持
static const double STOP_HOLD_DIST = 6.0;
static const double STOP_HOLD_PEN = 3.0;
static const double STOP_SEG_SWITCH_PEN = 80.0;

// 图节点合并阈值
static const double NODE_MERGE_THRESH = 3.0;

// 分链阈值，单位：秒
static const double CHAIN_RESET_GAP_SEC = 5.0;

// =====================================================
// 功能开关
// =====================================================

// 是否启用虚拟车道：
// true  时 MM 从 trail.cpp 读取虚拟车道偏移，并调用 trail.cpp 输出 virtual_lane_lines.kml/txt；
// false 时只使用道路中心线 offset = 0.0。
static bool g_MMUseVirtualLaneOffset = true;

// SD 模式位掩码，与 SD.cpp 中定义保持一致。
// 0  : SD 不影响 MM
// 1  : SD 预测残差影响搜索中心
// 2  : SD 预测残差调整候选发射代价
// 4  : SD 预测步长进入 Viterbi 转移代价
// 8  : SD 速度方向约束候选道路方向
// 15 : 全部启用
static const int SPEED_SD_MODE_OFF = 0;
static const int SPEED_SD_MODE_CENTER = 1 << 0;
static const int SPEED_SD_MODE_COST = 1 << 1;
static const int SPEED_SD_MODE_TRANS = 1 << 2;
static const int SPEED_SD_MODE_HEADING = 1 << 3;
static const int SPEED_SD_MODE_FULL =
SPEED_SD_MODE_CENTER |
SPEED_SD_MODE_COST |
SPEED_SD_MODE_TRANS |
SPEED_SD_MODE_HEADING;


// SD 总开关 + 模式选择，由 MapMarching() 包装接口在运行时设置。
static bool g_MMEnableSpeedSD = false;
static int g_MMSpeedSDMode = SPEED_SD_MODE_FULL;


// ---------------------------- 基础工具函数 ----------------------------

static double sqr(double x)
{
    return x * x;
}

static double clampDouble(double x, double lo, double hi)
{
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

static double angleDiff(double a, double b)
{
    double d = a - b;

    while (d > PI) d -= 2.0 * PI;
    while (d < -PI) d += 2.0 * PI;

    return std::fabs(d);
}

static double roadHeadingDiffBidirectional(double heading_a, double heading_b)
{
    double d1 = angleDiff(heading_a, heading_b);
    double d2 = angleDiff(heading_a, heading_b + PI);

    return std::min(d1, d2);
}

static void ll2xy(
    double lat_deg,
    double lon_deg,
    const LocalRef& ref,
    double& x,
    double& y)
{
    double lat_rad = lat_deg * PI / 180.0;
    double lon_rad = lon_deg * PI / 180.0;

    double dLat = lat_rad - ref.lat0_rad;
    double dLon = lon_rad - ref.lon0_rad;

    y = R_EARTH * dLat;
    x = R_EARTH * dLon * ref.cosLat0;
}

static void xy2ll(
    double x,
    double y,
    const LocalRef& ref,
    double& lat_deg,
    double& lon_deg)
{
    double lat_rad = y / R_EARTH + ref.lat0_rad;
    double lon_rad = x / (R_EARTH * ref.cosLat0) + ref.lon0_rad;

    lat_deg = lat_rad * 180.0 / PI;
    lon_deg = lon_rad * 180.0 / PI;
}


// ---------------------------- KML 解析与路网构建 ----------------------------

static bool isParallelAndClose(
    const RoadSegmentXY& a,
    const RoadSegmentXY& b,
    double headingThreshRad,
    double distThresh)
{
    double dh = angleDiff(a.heading, b.heading);
    if (dh > headingThreshRad) return false;

    double mx1 = 0.5 * (a.x1 + a.x2);
    double my1 = 0.5 * (a.y1 + a.y2);
    double mx2 = 0.5 * (b.x1 + b.x2);
    double my2 = 0.5 * (b.y1 + b.y2);

    double dx = mx1 - mx2;
    double dy = my1 - my2;
    double d = std::sqrt(dx * dx + dy * dy);

    return d <= distThresh;
}

static double minEndpointDist2(
    const RoadSegmentXY& a,
    const RoadSegmentXY& b)
{
    double d1 = sqr(a.x2 - b.x1) + sqr(a.y2 - b.y1);
    double d2 = sqr(a.x1 - b.x2) + sqr(a.y1 - b.y2);
    double d3 = sqr(a.x1 - b.x1) + sqr(a.y1 - b.y1);
    double d4 = sqr(a.x2 - b.x2) + sqr(a.y2 - b.y2);

    return std::min(std::min(d1, d2), std::min(d3, d4));
}

static std::vector<RoadSegmentXY> parseRoadSegmentsFromKML(
    const std::string& kml_file,
    LocalRef& localRef,
    std::vector<std::vector<int> >& adjacency)
{
    std::vector<RoadSegmentXY> segments;

    std::ifstream file(kml_file.c_str());
    if (!file.is_open())
    {
        std::cerr << "无法打开 KML 文件: " << kml_file << std::endl;
        return segments;
    }

    std::string content(
        (std::istreambuf_iterator<char>(file)),
        std::istreambuf_iterator<char>());

    file.close();

    std::regex re("<coordinates>([^<]+)</coordinates>");
    std::smatch m;
    std::string::const_iterator it = content.cbegin();

    std::vector<std::pair<double, double> > allPoints;

    while (std::regex_search(it, content.cend(), m, re))
    {
        std::stringstream ss(m[1].str());
        std::string coord;
        std::vector<std::pair<double, double> > pts;

        while (ss >> coord)
        {
            size_t c1 = coord.find(',');
            size_t c2 = coord.find(',', c1 + 1);

            if (c1 == std::string::npos || c2 == std::string::npos)
                continue;

            double lon = std::stod(coord.substr(0, c1));
            double lat = std::stod(coord.substr(c1 + 1, c2 - c1 - 1));

            pts.push_back(std::make_pair(lat, lon));
            allPoints.push_back(std::make_pair(lat, lon));
        }

        if (pts.size() >= 2)
        {
            for (size_t i = 1; i < pts.size(); ++i)
            {
                RoadSegmentXY seg;
                seg.lat1 = pts[i - 1].first;
                seg.lon1 = pts[i - 1].second;
                seg.lat2 = pts[i].first;
                seg.lon2 = pts[i].second;
                seg.x1 = seg.y1 = seg.x2 = seg.y2 = 0.0;
                seg.id = -1;
                seg.clusterId = -1;
                seg.heading = 0.0;

                segments.push_back(seg);
            }
        }

        it = m.suffix().first;
    }

    if (segments.empty())
    {
        std::cerr << "KML 中未解析到道路线段" << std::endl;
        return segments;
    }

    double sumLat = 0.0;
    double sumLon = 0.0;

    for (size_t i = 0; i < allPoints.size(); ++i)
    {
        sumLat += allPoints[i].first;
        sumLon += allPoints[i].second;
    }

    double lat0 = sumLat / (double)allPoints.size();
    double lon0 = sumLon / (double)allPoints.size();

    localRef.lat0_rad = lat0 * PI / 180.0;
    localRef.lon0_rad = lon0 * PI / 180.0;
    localRef.cosLat0 = std::cos(localRef.lat0_rad);

    for (size_t i = 0; i < segments.size(); ++i)
    {
        RoadSegmentXY& seg = segments[i];

        ll2xy(seg.lat1, seg.lon1, localRef, seg.x1, seg.y1);
        ll2xy(seg.lat2, seg.lon2, localRef, seg.x2, seg.y2);

        seg.id = (int)i;

        double dx = seg.x2 - seg.x1;
        double dy = seg.y2 - seg.y1;
        seg.heading = std::atan2(dy, dx);
    }

    size_t N = segments.size();
    adjacency.assign(N, std::vector<int>());

    const double CONNECT_THRESH = 8.0;
    double thresh2 = CONNECT_THRESH * CONNECT_THRESH;

    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = i + 1; j < N; ++j)
        {
            if (minEndpointDist2(segments[i], segments[j]) <= thresh2)
            {
                adjacency[i].push_back((int)j);
                adjacency[j].push_back((int)i);
            }
        }
    }

    const double HEADING_THRESH = 30.0 * PI / 180.0;
    const double CLUSTER_DIST = 25.0;

    int clusterCount = 0;

    for (size_t i = 0; i < N; ++i)
    {
        if (segments[i].clusterId != -1)
            continue;

        int cid = clusterCount++;
        segments[i].clusterId = cid;

        std::queue<size_t> q;
        q.push(i);

        while (!q.empty())
        {
            size_t cur = q.front();
            q.pop();

            for (size_t j = 0; j < N; ++j)
            {
                if (segments[j].clusterId != -1)
                    continue;

                if (isParallelAndClose(
                    segments[cur],
                    segments[j],
                    HEADING_THRESH,
                    CLUSTER_DIST))
                {
                    segments[j].clusterId = cid;
                    q.push(j);
                }
            }
        }
    }

    std::cout << "KML 解析完成，道路线段数 = " << segments.size()
        << "，道路簇数 = " << clusterCount << std::endl;

    return segments;
}


// ---------------------------- 路网图 ----------------------------

struct GraphEdge
{
    int to = -1;
    double w = 0.0;
};

struct RoadSegmentExt
{
    RoadSegmentXY seg;
    int u = -1;
    int v = -1;
    double len = 0.0;
};

static int findOrCreateNode(
    double x,
    double y,
    std::vector<std::pair<double, double> >& nodeXY)
{
    for (int i = 0; i < (int)nodeXY.size(); ++i)
    {
        double dx = x - nodeXY[i].first;
        double dy = y - nodeXY[i].second;

        if (dx * dx + dy * dy <= NODE_MERGE_THRESH * NODE_MERGE_THRESH)
            return i;
    }

    nodeXY.push_back(std::make_pair(x, y));
    return (int)nodeXY.size() - 1;
}

static void buildRoadGraph(
    const std::vector<RoadSegmentXY>& rawSegs,
    std::vector<RoadSegmentExt>& segs,
    std::vector<std::vector<GraphEdge> >& graph,
    std::vector<std::pair<double, double> >& nodeXY)
{
    segs.clear();
    nodeXY.clear();

    for (size_t i = 0; i < rawSegs.size(); ++i)
    {
        RoadSegmentExt e;
        e.seg = rawSegs[i];

        double dx = e.seg.x2 - e.seg.x1;
        double dy = e.seg.y2 - e.seg.y1;

        e.len = std::sqrt(dx * dx + dy * dy);

        e.u = findOrCreateNode(e.seg.x1, e.seg.y1, nodeXY);
        e.v = findOrCreateNode(e.seg.x2, e.seg.y2, nodeXY);

        segs.push_back(e);
    }

    graph.assign(nodeXY.size(), std::vector<GraphEdge>());

    for (size_t i = 0; i < segs.size(); ++i)
    {
        const RoadSegmentExt& e = segs[i];

        if (e.u == e.v || e.len < 1e-8)
            continue;

        graph[e.u].push_back({ e.v, e.len });
        graph[e.v].push_back({ e.u, e.len });
    }
}

static std::vector<double> runDijkstra(
    int src,
    const std::vector<std::vector<GraphEdge> >& graph)
{
    const double INF = std::numeric_limits<double>::infinity();

    std::vector<double> dist(graph.size(), INF);

    typedef std::pair<double, int> PDI;
    std::priority_queue<PDI, std::vector<PDI>, std::greater<PDI> > pq;

    dist[src] = 0.0;
    pq.push(std::make_pair(0.0, src));

    while (!pq.empty())
    {
        PDI cur = pq.top();
        pq.pop();

        double d = cur.first;
        int u = cur.second;

        if (d > dist[u])
            continue;

        for (size_t i = 0; i < graph[u].size(); ++i)
        {
            int v = graph[u][i].to;
            double nd = d + graph[u][i].w;

            if (nd < dist[v])
            {
                dist[v] = nd;
                pq.push(std::make_pair(nd, v));
            }
        }
    }

    return dist;
}

static double shortestNodeDistance(
    int src,
    int dst,
    const std::vector<std::vector<GraphEdge> >& graph,
    std::unordered_map<int, std::vector<double> >& cache)
{
    if (src == dst)
        return 0.0;

    if (src < 0 || dst < 0)
        return std::numeric_limits<double>::infinity();

    std::unordered_map<int, std::vector<double> >::iterator it = cache.find(src);

    if (it == cache.end())
    {
        std::vector<double> dist = runDijkstra(src, graph);
        cache[src] = dist;
        it = cache.find(src);
    }

    if (dst >= (int)it->second.size())
        return std::numeric_limits<double>::infinity();

    return it->second[dst];
}


// ---------------------------- 候选结构 ----------------------------

struct SeqCandidate
{
    int seg_id = -1;
    int cluster_id = -1;

    double x = 0.0;
    double y = 0.0;

    double xc = 0.0;
    double yc = 0.0;
    double tproj = 0.0;
    double seg_len = 0.0;
    double s = 0.0;

    double lane_offset = 0.0;
    double d_obs = 0.0;
    double d_dr = 0.0;
    double dir_diff = 0.0;

    double emit_cost = 0.0;
};

struct ProvisionalState
{
    bool valid = false;
    int seg_id = -1;
    int cluster_id = -1;
    double x = 0.0;
    double y = 0.0;
    double s = 0.0;
    double lane_offset = 0.0;
};

struct EpochMeta
{
    bool valid = false;
    long long timestamp = 0;

    double obsx = 0.0;
    double obsy = 0.0;
    double obs_step = 0.0;
    double pred_step = 0.0;

    MMState state = MM_INIT;
};


// ---------------------------- 空间搜索 ----------------------------

static bool projectPointToSegment(
    double x,
    double y,
    const RoadSegmentXY& seg,
    double& xp,
    double& yp,
    double& t,
    double& segLen)
{
    double dx = seg.x2 - seg.x1;
    double dy = seg.y2 - seg.y1;
    double len2 = dx * dx + dy * dy;

    if (len2 < 1e-8)
    {
        segLen = 0.0;
        return false;
    }

    t = ((x - seg.x1) * dx + (y - seg.y1) * dy) / len2;
    t = clampDouble(t, 0.0, 1.0);

    xp = seg.x1 + t * dx;
    yp = seg.y1 + t * dy;

    segLen = std::sqrt(len2);

    return true;
}

static bool projectPointToSegment(
    double x,
    double y,
    const RoadSegmentXY& seg,
    double& xp,
    double& yp,
    double& t)
{
    double segLen = 0.0;

    return projectPointToSegment(
        x,
        y,
        seg,
        xp,
        yp,
        t,
        segLen);
}

static void collectAllNearSegments(
    const std::vector<RoadSegmentXY>& segments,
    double x,
    double y,
    double radius,
    std::vector<int>& outIds)
{
    outIds.clear();

    double r2 = radius * radius;

    for (size_t i = 0; i < segments.size(); ++i)
    {
        double xp = 0.0;
        double yp = 0.0;
        double t = 0.0;

        if (!projectPointToSegment(x, y, segments[i], xp, yp, t))
            continue;

        double d2 = sqr(xp - x) + sqr(yp - y);

        if (d2 <= r2)
            outIds.push_back((int)i);
    }
}

static bool isNearIntersection(
    const std::vector<RoadSegmentXY>& segments,
    const std::vector<std::vector<int> >& adjacency,
    int segId,
    double x,
    double y)
{
    if (segId < 0 || segId >= (int)segments.size())
        return false;

    if (segId >= (int)adjacency.size())
        return false;

    const RoadSegmentXY& seg = segments[segId];

    double d1 = std::sqrt(sqr(x - seg.x1) + sqr(y - seg.y1));
    double d2 = std::sqrt(sqr(x - seg.x2) + sqr(y - seg.y2));

    if ((d1 <= INTERSECTION_NEAR_DIST || d2 <= INTERSECTION_NEAR_DIST) &&
        !adjacency[segId].empty())
        return true;

    return false;
}

static void collectCandidateSegments(
    const std::vector<RoadSegmentXY>& segments,
    const std::vector<std::vector<int> >& adjacency,
    MMState mm_state,
    const ProvisionalState& prev,
    double x_center,
    double y_center,
    bool is_static_epoch,
    std::vector<int>& outIds)
{
    outIds.clear();

    double radius = INIT_SEARCH_RADIUS;

    if (mm_state == MM_ON_ROAD)
        radius = ONROAD_SEARCH_RADIUS;
    else if (mm_state == MM_AT_INTERSECTION)
        radius = INTERSECTION_SEARCH_RADIUS;
    else if (mm_state == MM_RECOVERY)
        radius = INIT_SEARCH_RADIUS;

    if (is_static_epoch)
        radius = std::min(radius, STOP_SEARCH_RADIUS);

    std::set<int> idset;

    std::vector<int> nearAll;
    collectAllNearSegments(segments, x_center, y_center, radius, nearAll);

    for (size_t i = 0; i < nearAll.size(); ++i)
        idset.insert(nearAll[i]);

    if (prev.valid && mm_state != MM_INIT && mm_state != MM_RECOVERY)
    {
        if (prev.seg_id >= 0 && prev.seg_id < (int)segments.size())
            idset.insert(prev.seg_id);

        if (prev.seg_id >= 0 && prev.seg_id < (int)adjacency.size())
        {
            for (size_t i = 0; i < adjacency[prev.seg_id].size(); ++i)
                idset.insert(adjacency[prev.seg_id][i]);
        }
    }

    for (std::set<int>::iterator it = idset.begin(); it != idset.end(); ++it)
        outIds.push_back(*it);

    if (outIds.empty())
    {
        for (size_t i = 0; i < segments.size(); ++i)
            outIds.push_back(segments[i].id);
    }
}


// ---------------------------- 候选生成 ----------------------------

static std::vector<SeqCandidate> buildCandidatesForEpoch(
    const CoorData& c,
    const std::vector<RoadSegmentXY>& segments,
    const std::vector<int>& candidateSegIds,
    MMState mm_state,
    bool has_prev,
    const ProvisionalState& prevState,
    double x_obs,
    double y_obs,
    const SDMMResult& sd_epoch)
{
    std::vector<SeqCandidate> cands;

    double W_OBS = W_OBS_INIT;
    double W_DR = W_DR_INIT;
    double W_CONT = W_CONT_INIT;
    double W_DIR = W_DIR_INIT;
    double W_VEL_DIR = W_VEL_DIR_INIT;

    if (mm_state == MM_ON_ROAD)
    {
        W_OBS = W_OBS_ONROAD;
        W_DR = W_DR_ONROAD;
        W_CONT = W_CONT_ONROAD;
        W_DIR = W_DIR_ONROAD;
        W_VEL_DIR = W_VEL_DIR_ONROAD;
    }
    else if (mm_state == MM_AT_INTERSECTION)
    {
        W_OBS = W_OBS_INT;
        W_DR = W_DR_INT;
        W_CONT = W_CONT_INT;
        W_DIR = W_DIR_INT;
        W_VEL_DIR = W_VEL_DIR_INT;
    }

    SDMMRequest sd_cand_req;
    sd_cand_req.stage = SD_MM_STAGE_CANDIDATE;
    sd_cand_req.use_sd = g_MMEnableSpeedSD;
    sd_cand_req.mode = g_MMSpeedSDMode;
    sd_cand_req.has_prev = has_prev;
    sd_cand_req.is_static_epoch = sd_epoch.is_static_epoch;
    sd_cand_req.use_vel_heading = sd_epoch.use_vel_heading;
    sd_cand_req.x_obs = x_obs;
    sd_cand_req.y_obs = y_obs;
    sd_cand_req.x_dr = sd_epoch.x_dr;
    sd_cand_req.y_dr = sd_epoch.y_dr;
    sd_cand_req.vE_meas = sd_epoch.vE;
    sd_cand_req.vN_meas = sd_epoch.vN;
    sd_cand_req.W_OBS = W_OBS;
    sd_cand_req.W_DR = W_DR;
    sd_cand_req.W_CONT = W_CONT;
    sd_cand_req.W_VEL_DIR = W_VEL_DIR;

    SDMMResult sd_cand_base;
    SD_MMEvaluate(&c, sd_cand_req, sd_cand_base);

    W_OBS = sd_cand_base.W_OBS;
    W_DR = sd_cand_base.W_DR;
    W_CONT = sd_cand_base.W_CONT;
    W_VEL_DIR = sd_cand_base.W_VEL_DIR;

    for (size_t cs = 0; cs < candidateSegIds.size(); ++cs)
    {
        int sid = candidateSegIds[cs];

        if (sid < 0 || sid >= (int)segments.size())
            continue;

        const RoadSegmentXY& seg = segments[sid];

        double xc = 0.0;
        double yc = 0.0;
        double tproj = 0.0;
        double segLen = 0.0;

        if (!projectPointToSegment(x_obs, y_obs, seg, xc, yc, tproj, segLen))
            continue;

        double vx = seg.x2 - seg.x1;
        double vy = seg.y2 - seg.y1;
        double len = std::sqrt(vx * vx + vy * vy);

        if (len < 1e-8)
            continue;

        double nx = -vy / len;
        double ny = vx / len;

        int laneCount = GetVirtualLaneOffsetCount(g_MMUseVirtualLaneOffset);

        for (int k = 0; k < laneCount; ++k)
        {
            double lane_offset = GetVirtualLaneOffset(g_MMUseVirtualLaneOffset, k);

            double xp = xc + lane_offset * nx;
            double yp = yc + lane_offset * ny;

            double d2_obs = sqr(xp - x_obs) + sqr(yp - y_obs);

            if (d2_obs > sqr(MAX_SNAP_DIST))
                continue;

            double d2_dr =
                sqr(xp - sd_epoch.x_dr) +
                sqr(yp - sd_epoch.y_dr);

            double d2_cont = 0.0;

            if (has_prev && prevState.valid)
                d2_cont = sqr(xp - prevState.x) + sqr(yp - prevState.y);

            double dirDiff = 0.0;
            double dirCost = 0.0;

            if (has_prev && prevState.valid && !sd_epoch.is_static_epoch)
            {
                double dxm = xp - prevState.x;
                double dym = yp - prevState.y;
                double norm2 = dxm * dxm + dym * dym;

                if (norm2 > 1.0)
                {
                    double veh_heading = std::atan2(dym, dxm);

                    dirDiff = roadHeadingDiffBidirectional(
                        veh_heading,
                        seg.heading);

                    dirCost = (1.0 - std::cos(dirDiff)) * W_DIR;
                }
            }

            sd_cand_req.road_heading = seg.heading;
            sd_cand_req.W_OBS = W_OBS;
            sd_cand_req.W_DR = W_DR;
            sd_cand_req.W_CONT = W_CONT;
            sd_cand_req.W_VEL_DIR = W_VEL_DIR;

            SDMMResult sd_cand_each;
            SD_MMEvaluate(&c, sd_cand_req, sd_cand_each);

            double velDirCost = sd_cand_each.vel_dir_cost;

            double stopHoldCost = 0.0;

            if (sd_epoch.is_static_epoch && has_prev && prevState.valid)
            {
                double d2_hold =
                    sqr(xp - prevState.x) +
                    sqr(yp - prevState.y);

                stopHoldCost += STOP_HOLD_PEN * d2_hold;

                if (std::sqrt(d2_hold) > STOP_HOLD_DIST)
                    stopHoldCost += 60.0;

                if (prevState.seg_id >= 0 && seg.id != prevState.seg_id)
                    stopHoldCost += STOP_SEG_SWITCH_PEN;
            }

            SeqCandidate cand;
            cand.seg_id = seg.id;
            cand.cluster_id = seg.clusterId;
            cand.x = xp;
            cand.y = yp;
            cand.xc = xc;
            cand.yc = yc;
            cand.tproj = tproj;
            cand.seg_len = segLen;
            cand.s = tproj * segLen;
            cand.lane_offset = lane_offset;
            cand.d_obs = std::sqrt(d2_obs);
            cand.d_dr = std::sqrt(d2_dr);
            cand.dir_diff = dirDiff;

            cand.emit_cost =
                W_OBS * d2_obs +
                W_DR * d2_dr +
                W_CONT * d2_cont +
                dirCost +
                velDirCost +
                stopHoldCost;

            cands.push_back(cand);
        }
    }

    std::sort(
        cands.begin(),
        cands.end(),
        [](const SeqCandidate& a, const SeqCandidate& b)
        {
            return a.emit_cost < b.emit_cost;
        });

    if ((int)cands.size() > TOPK_PER_EPOCH)
        cands.resize(TOPK_PER_EPOCH);

    return cands;
}


// ---------------------------- 转移代价 / Viterbi ----------------------------

static double routeDistanceBetweenCandidates(
    const SeqCandidate& a,
    const SeqCandidate& b,
    const std::vector<RoadSegmentExt>& segs,
    const std::vector<std::vector<GraphEdge> >& graph,
    std::unordered_map<int, std::vector<double> >& cache)
{
    if (a.seg_id < 0 || b.seg_id < 0)
        return std::numeric_limits<double>::infinity();

    if (a.seg_id >= (int)segs.size() || b.seg_id >= (int)segs.size())
        return std::numeric_limits<double>::infinity();

    const RoadSegmentExt& sa = segs[a.seg_id];
    const RoadSegmentExt& sb = segs[b.seg_id];

    if (a.seg_id == b.seg_id)
        return std::fabs(b.s - a.s);

    double a_to_u = a.s;
    double a_to_v = sa.len - a.s;

    double b_from_u = b.s;
    double b_from_v = sb.len - b.s;

    double d_uu = shortestNodeDistance(sa.u, sb.u, graph, cache);
    double d_uv = shortestNodeDistance(sa.u, sb.v, graph, cache);
    double d_vu = shortestNodeDistance(sa.v, sb.u, graph, cache);
    double d_vv = shortestNodeDistance(sa.v, sb.v, graph, cache);

    double ans = std::numeric_limits<double>::infinity();

    if (std::isfinite(d_uu))
        ans = std::min(ans, a_to_u + d_uu + b_from_u);

    if (std::isfinite(d_uv))
        ans = std::min(ans, a_to_u + d_uv + b_from_v);

    if (std::isfinite(d_vu))
        ans = std::min(ans, a_to_v + d_vu + b_from_u);

    if (std::isfinite(d_vv))
        ans = std::min(ans, a_to_v + d_vv + b_from_v);

    return ans;
}

static bool isConnected(
    int a,
    int b,
    const std::vector<std::vector<int> >& adjacency)
{
    if (a == b)
        return true;

    if (a < 0 || b < 0)
        return false;

    if (a >= (int)adjacency.size())
        return false;

    for (size_t i = 0; i < adjacency[a].size(); ++i)
    {
        if (adjacency[a][i] == b)
            return true;
    }

    return false;
}

static double computeTransitionCost(
    const SeqCandidate& prevC,
    const SeqCandidate& curC,
    double obs_step,
    double pred_step,
    const std::vector<std::vector<int> >& adjacency,
    const std::vector<RoadSegmentExt>& segs,
    const std::vector<std::vector<GraphEdge> >& graph,
    std::unordered_map<int, std::vector<double> >& cache)
{
    double cost = 0.0;

    double routeDist = routeDistanceBetweenCandidates(
        prevC,
        curC,
        segs,
        graph,
        cache);

    if (!std::isfinite(routeDist))
        return 1e9;

    cost += std::fabs(routeDist - obs_step) / TRANS_ROUTE_BETA;

    SDMMRequest sd_trans_req;
    sd_trans_req.stage = SD_MM_STAGE_TRANSITION;
    sd_trans_req.use_sd = g_MMEnableSpeedSD;
    sd_trans_req.mode = g_MMSpeedSDMode;
    sd_trans_req.has_pred_step = true;
    sd_trans_req.pred_step = pred_step;
    sd_trans_req.route_dist = routeDist;
    sd_trans_req.trans_pred_beta = TRANS_PRED_BETA;

    SDMMResult sd_trans_out;
    SD_MMEvaluate(nullptr, sd_trans_req, sd_trans_out);

    cost += sd_trans_out.trans_pred_cost;

    if (g_MMUseVirtualLaneOffset)
        cost += TRANS_LANE_W * std::fabs(curC.lane_offset - prevC.lane_offset);

    if (prevC.seg_id == curC.seg_id)
    {
        if (pred_step > 1.5 && curC.s + 2.0 < prevC.s)
            cost += PEN_REVERSE_TRANS;
    }
    else
    {
        if (!isConnected(prevC.seg_id, curC.seg_id, adjacency))
            cost += PEN_NON_ADJ_TRANS;

        if (prevC.cluster_id != curC.cluster_id)
            cost += PEN_CLUSTER_TRANS;
    }

    return cost;
}

static std::vector<int> runChainViterbi(
    const std::vector<std::vector<SeqCandidate> >& allCands,
    const std::vector<EpochMeta>& meta,
    int L,
    int R,
    const std::vector<std::vector<int> >& adjacency,
    const std::vector<RoadSegmentExt>& segs,
    const std::vector<std::vector<GraphEdge> >& graph)
{
    int T = R - L + 1;

    std::vector<std::vector<double> > dp(T);
    std::vector<std::vector<int> > parent(T);

    for (int t = 0; t < T; ++t)
    {
        int idx = L + t;
        int K = (int)allCands[idx].size();

        dp[t].assign(K, 1e18);
        parent[t].assign(K, -1);
    }

    for (int k = 0; k < (int)allCands[L].size(); ++k)
        dp[0][k] = allCands[L][k].emit_cost;

    std::unordered_map<int, std::vector<double> > distCache;

    for (int t = 1; t < T; ++t)
    {
        int idxPrev = L + t - 1;
        int idxCur = L + t;

        for (int j = 0; j < (int)allCands[idxCur].size(); ++j)
        {
            const SeqCandidate& curC = allCands[idxCur][j];

            double best = 1e18;
            int bestp = -1;

            for (int i = 0; i < (int)allCands[idxPrev].size(); ++i)
            {
                if (dp[t - 1][i] >= 1e17)
                    continue;

                const SeqCandidate& prevC = allCands[idxPrev][i];

                double transCost = computeTransitionCost(
                    prevC,
                    curC,
                    meta[idxCur].obs_step,
                    meta[idxCur].pred_step,
                    adjacency,
                    segs,
                    graph,
                    distCache);

                double val =
                    dp[t - 1][i] +
                    transCost +
                    curC.emit_cost;

                if (val < best)
                {
                    best = val;
                    bestp = i;
                }
            }

            dp[t][j] = best;
            parent[t][j] = bestp;
        }
    }

    std::vector<int> chosen(T, -1);

    double best = 1e18;
    int bestk = -1;

    for (int k = 0; k < (int)dp[T - 1].size(); ++k)
    {
        if (dp[T - 1][k] < best)
        {
            best = dp[T - 1][k];
            bestk = k;
        }
    }

    chosen[T - 1] = bestk;

    for (int t = T - 1; t >= 1; --t)
    {
        if (chosen[t] < 0)
            break;

        chosen[t - 1] = parent[t][chosen[t]];
    }

    return chosen;
}


// ---------------------------- 状态机辅助 ----------------------------

static MMState updateState(
    MMState state,
    const SeqCandidate& best,
    bool nearIntersection,
    int& badCounter)
{
    bool bad =
        best.d_obs > BAD_OBS_DIST ||
        best.dir_diff > BAD_DIR_DIFF ||
        best.emit_cost > BAD_COST_TH;

    if (bad)
        badCounter++;
    else
        badCounter = 0;

    if (state == MM_INIT)
    {
        if (!bad)
            return nearIntersection ? MM_AT_INTERSECTION : MM_ON_ROAD;

        return MM_INIT;
    }

    if (state == MM_ON_ROAD)
    {
        if (badCounter >= BAD_LIMIT)
            return MM_RECOVERY;

        if (nearIntersection)
            return MM_AT_INTERSECTION;

        return MM_ON_ROAD;
    }

    if (state == MM_AT_INTERSECTION)
    {
        if (badCounter >= BAD_LIMIT)
            return MM_RECOVERY;

        if (!nearIntersection)
            return MM_ON_ROAD;

        return MM_AT_INTERSECTION;
    }

    if (state == MM_RECOVERY)
    {
        if (!bad)
            return nearIntersection ? MM_AT_INTERSECTION : MM_ON_ROAD;

        return MM_RECOVERY;
    }

    return MM_INIT;
}


// ---------------------------- 输出 ----------------------------

static void writeMatchedTxt(
    const std::string& out_file,
    const std::vector<CoorData>& coords,
    const std::vector<SeqCandidate>& chosen,
    const LocalRef& ref)
{
    std::ofstream fout(out_file.c_str());

    if (!fout.is_open())
    {
        std::cerr << "无法写入匹配结果文件: " << out_file << std::endl;
        return;
    }

    fout << std::fixed << std::setprecision(10);

    for (size_t i = 0; i < coords.size(); ++i)
    {
        double lat = coords[i].BLH[0];
        double lon = coords[i].BLH[1];
        double h = coords[i].BLH[2];

        if (i < chosen.size() && chosen[i].seg_id >= 0)
        {
            xy2ll(chosen[i].x, chosen[i].y, ref, lat, lon);
        }

        // ComputeCoorTimestampSec() 返回秒；GPS_matched 输出使用毫秒。
        long long ts = ComputeCoorTimestampSec(coords[i]) * 1000LL;

        fout << "GPS_matched "
            << ts << " "
            << lat << " "
            << lon << " "
            << h << " "
            << 0.5 << "\n";
    }

    fout.close();
}

static void writeMatchedKML(
    const std::string& out_file,
    const std::vector<CoorData>& coords,
    const std::vector<SeqCandidate>& chosen,
    const LocalRef& ref)
{
    std::string kml_out = out_file;

    size_t pos = kml_out.find_last_of('.');

    if (pos != std::string::npos)
        kml_out = kml_out.substr(0, pos) + ".kml";
    else
        kml_out += ".kml";

    std::ofstream fout(kml_out.c_str());

    if (!fout.is_open())
    {
        std::cerr << "无法写入 KML 文件: " << kml_out << std::endl;
        return;
    }

    fout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    fout << "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n";
    fout << "<Document>\n";
    fout << "<name>MapMarching Matched Trajectory</name>\n";
    fout << "<Style id=\"matchedLine\"><LineStyle><color>ff00ff00</color><width>4</width></LineStyle></Style>\n";
    fout << "<Placemark>\n";
    fout << "<name>matched trajectory</name>\n";
    fout << "<styleUrl>#matchedLine</styleUrl>\n";
    fout << "<LineString><tessellate>1</tessellate><coordinates>\n";

    fout << std::fixed << std::setprecision(10);

    for (size_t i = 0; i < coords.size(); ++i)
    {
        double lat = coords[i].BLH[0];
        double lon = coords[i].BLH[1];
        double h = coords[i].BLH[2];

        if (i < chosen.size() && chosen[i].seg_id >= 0)
        {
            xy2ll(chosen[i].x, chosen[i].y, ref, lat, lon);
        }

        fout << lon << "," << lat << "," << h << "\n";
    }

    fout << "</coordinates></LineString>\n";
    fout << "</Placemark>\n";
    fout << "</Document>\n";
    fout << "</kml>\n";

    fout.close();

    std::cout << "匹配 KML 输出完成: " << kml_out << std::endl;
}


// ---------------------------- 主函数 ----------------------------

/**
 * 函数功能：
 *     状态机地图匹配核心流程。
 *     负责解析路网、构建候选点、运行状态机、执行 Viterbi 序列优化，并写出匹配结果。
 *
 * 输入：
 *     coords：GNSS 原始观测数组，函数会根据 SD/MM 结果写回质量字段。
 *     kml_file：道路 KML 文件路径。
 *     out_file：匹配结果 TXT 输出路径；函数会同步写出同名 KML。
 *
 * 输出：
 *     无返回值；结果写入 out_file 和同名 .kml 文件。
 *
 * 调用位置：
 *     MapMarching 对外入口内部调用。
 *
 * 数据流作用：
 *     CoorData + Road KML -> out\mm_matched_points.txt/kml。
 */
void RunStateMachineMapMatching(
    std::vector<CoorData>& coords,
    const std::string& kml_file,
    const std::string& out_file)
{
    if (coords.empty())
    {
        std::cerr << "MapMatch 输入坐标为空" << std::endl;
        return;
    }

    LocalRef ref;
    std::vector<std::vector<int> > adjacency;

    std::vector<RoadSegmentXY> segments =
        parseRoadSegmentsFromKML(kml_file, ref, adjacency);

    if (segments.empty())
    {
        std::cerr << "道路为空，无法地图匹配" << std::endl;
        return;
    }

    std::vector<RoadSegmentExt> segExt;
    std::vector<std::vector<GraphEdge> > graph;
    std::vector<std::pair<double, double> > nodeXY;

    buildRoadGraph(segments, segExt, graph, nodeXY);

    std::vector<std::vector<SeqCandidate> > allCands(coords.size());
    std::vector<EpochMeta> meta(coords.size());

    MMState mm_state = MM_INIT;
    ProvisionalState prevState;

    int badCounter = 0;

    double vx_last = 0.0;
    double vy_last = 0.0;
    bool has_vel = false;

    bool has_prev = false;

    for (size_t idx = 0; idx < coords.size(); ++idx)
    {
        CoorData& c = coords[idx];

        double x_obs = 0.0;
        double y_obs = 0.0;

        ll2xy(c.BLH[0], c.BLH[1], ref, x_obs, y_obs);

        // ComputeCoorTimestampSec() 返回秒，因此这里不除以 1000。
        double dt = 1.0;

        if (idx > 0)
        {
            long long t1 = ComputeCoorTimestampSec(coords[idx - 1]);
            long long t2 = ComputeCoorTimestampSec(coords[idx]);

            double dts = (double)(t2 - t1);

            if (std::isfinite(dts) && dts > 0.05 && dts < 10.0)
                dt = dts;
        }

        SDMMRequest sd_epoch_req;
        sd_epoch_req.stage = SD_MM_STAGE_EPOCH;
        sd_epoch_req.use_sd = g_MMEnableSpeedSD;
        sd_epoch_req.mode = g_MMSpeedSDMode;
        sd_epoch_req.has_prev = has_prev;
        sd_epoch_req.has_vel = has_vel;
        sd_epoch_req.prev_x = prevState.x;
        sd_epoch_req.prev_y = prevState.y;
        sd_epoch_req.vx_last = vx_last;
        sd_epoch_req.vy_last = vy_last;
        sd_epoch_req.dt = dt;
        sd_epoch_req.x_obs = x_obs;
        sd_epoch_req.y_obs = y_obs;

        SDMMResult sd_epoch;
        SD_MMEvaluate(&c, sd_epoch_req, sd_epoch);

        bool has_speed_meas = sd_epoch.has_speed;
        double speed_meas = sd_epoch.speed;
        bool is_static_epoch = sd_epoch.is_static_epoch;

        double x_center = sd_epoch.x_center;
        double y_center = sd_epoch.y_center;

        c.trusted = sd_epoch.trusted;
        c.trust_reason = sd_epoch.trust_reason;
        c.quality_score = sd_epoch.quality_score;

        std::vector<int> candidateSegIds;

        collectCandidateSegments(
            segments,
            adjacency,
            mm_state,
            prevState,
            x_center,
            y_center,
            is_static_epoch,
            candidateSegIds);

        std::vector<SeqCandidate> curCands = buildCandidatesForEpoch(
            c,
            segments,
            candidateSegIds,
            mm_state,
            has_prev,
            prevState,
            x_obs,
            y_obs,
            sd_epoch);

        if (curCands.empty())
        {
            std::vector<int> extra;

            collectAllNearSegments(
                segments,
                x_obs,
                y_obs,
                INIT_SEARCH_RADIUS,
                extra);

            curCands = buildCandidatesForEpoch(
                c,
                segments,
                extra,
                MM_INIT,
                has_prev,
                prevState,
                x_obs,
                y_obs,
                sd_epoch);
        }

        allCands[idx] = curCands;

        meta[idx].valid = !curCands.empty();
        meta[idx].timestamp = ComputeCoorTimestampSec(c);
        meta[idx].obsx = x_obs;
        meta[idx].obsy = y_obs;
        meta[idx].state = mm_state;

        if (idx == 0)
        {
            meta[idx].obs_step = 0.0;

            SDMMRequest sd_trans_req;
            sd_trans_req.stage = SD_MM_STAGE_TRANSITION;
            sd_trans_req.use_sd = g_MMEnableSpeedSD;
            sd_trans_req.mode = g_MMSpeedSDMode;
            sd_trans_req.has_speed_meas = has_speed_meas;
            sd_trans_req.speed_meas = speed_meas;
            sd_trans_req.dt = dt;
            sd_trans_req.fallback_step = 0.0;

            SDMMResult sd_trans_out;
            SD_MMEvaluate(nullptr, sd_trans_req, sd_trans_out);

            meta[idx].pred_step = sd_trans_out.pred_step;
        }
        else
        {
            meta[idx].obs_step = std::sqrt(
                sqr(x_obs - meta[idx - 1].obsx) +
                sqr(y_obs - meta[idx - 1].obsy));

            SDMMRequest sd_trans_req;
            sd_trans_req.stage = SD_MM_STAGE_TRANSITION;
            sd_trans_req.use_sd = g_MMEnableSpeedSD;
            sd_trans_req.mode = g_MMSpeedSDMode;
            sd_trans_req.has_speed_meas = has_speed_meas;
            sd_trans_req.speed_meas = speed_meas;
            sd_trans_req.dt = dt;
            sd_trans_req.fallback_step = meta[idx].obs_step;

            SDMMResult sd_trans_out;
            SD_MMEvaluate(nullptr, sd_trans_req, sd_trans_out);

            meta[idx].pred_step = sd_trans_out.pred_step;
        }

        if (!curCands.empty())
        {
            const SeqCandidate& best = curCands[0];

            bool nearInt = isNearIntersection(
                segments,
                adjacency,
                best.seg_id,
                best.x,
                best.y);

            mm_state = updateState(
                mm_state,
                best,
                nearInt,
                badCounter);

            double old_prev_x = prevState.x;
            double old_prev_y = prevState.y;

            ProvisionalState bestProv;
            bestProv.valid = true;
            bestProv.seg_id = best.seg_id;
            bestProv.cluster_id = best.cluster_id;
            bestProv.x = best.x;
            bestProv.y = best.y;
            bestProv.s = best.s;
            bestProv.lane_offset = best.lane_offset;

            SDMMRequest sd_update_req;
            sd_update_req.stage = SD_MM_STAGE_UPDATE;
            sd_update_req.use_sd = g_MMEnableSpeedSD;
            sd_update_req.mode = g_MMSpeedSDMode;
            sd_update_req.has_speed_meas = has_speed_meas;
            sd_update_req.has_prev = has_prev;
            sd_update_req.vE_meas = sd_epoch.vE;
            sd_update_req.vN_meas = sd_epoch.vN;
            sd_update_req.best_x = bestProv.x;
            sd_update_req.best_y = bestProv.y;
            sd_update_req.old_prev_x = old_prev_x;
            sd_update_req.old_prev_y = old_prev_y;
            sd_update_req.dt = dt;
            sd_update_req.vx_last = vx_last;
            sd_update_req.vy_last = vy_last;
            sd_update_req.has_vel = has_vel;

            SDMMResult sd_update_out;
            SD_MMEvaluate(nullptr, sd_update_req, sd_update_out);

            vx_last = sd_update_out.vx_last;
            vy_last = sd_update_out.vy_last;
            has_vel = sd_update_out.has_vel;

            prevState = bestProv;
            has_prev = true;
        }
        else
        {
            badCounter++;

            if (badCounter >= BAD_LIMIT)
                mm_state = MM_RECOVERY;
        }
    }

    std::vector<SeqCandidate> chosen(coords.size());

    int L = -1;

    for (int i = 0; i <= (int)coords.size(); ++i)
    {
        bool breakHere = false;

        if (i == (int)coords.size())
        {
            breakHere = true;
        }
        else if (!meta[i].valid || allCands[i].empty())
        {
            breakHere = true;
        }
        else if (i > 0 && meta[i - 1].valid)
        {
            // meta.timestamp 是秒，因此这里不除以 1000。
            double gapSec =
                (double)(meta[i].timestamp - meta[i - 1].timestamp);

            if (gapSec > CHAIN_RESET_GAP_SEC)
                breakHere = true;
        }

        if (!breakHere)
        {
            if (L < 0)
                L = i;

            continue;
        }

        if (L >= 0)
        {
            int R = i - 1;

            if (L == R)
            {
                chosen[L] = allCands[L][0];
            }
            else
            {
                std::vector<int> path = runChainViterbi(
                    allCands,
                    meta,
                    L,
                    R,
                    adjacency,
                    segExt,
                    graph);

                for (int t = 0; t < (int)path.size(); ++t)
                {
                    int gi = L + t;
                    int pk = path[t];

                    if (pk >= 0 && pk < (int)allCands[gi].size())
                        chosen[gi] = allCands[gi][pk];
                    else
                        chosen[gi] = allCands[gi][0];
                }
            }

            L = -1;
        }

        if (i < (int)coords.size() && (!meta[i].valid || allCands[i].empty()))
        {
            SeqCandidate raw;
            raw.seg_id = -1;
            raw.x = meta[i].obsx;
            raw.y = meta[i].obsy;

            chosen[i] = raw;
        }
    }

    writeMatchedTxt(out_file, coords, chosen, ref);
    writeMatchedKML(out_file, coords, chosen, ref);

    std::cout << "MapMarching Viterbi 匹配完成，输出: " << out_file << std::endl;
}

/**
 * 函数功能：
 *     状态机 MM 对外入口，负责调度核心匹配流程并返回统一 CarPoint 结果。
 *
 * 输入：
 *     road：路网数据，sourcePath 指向实际用于匹配的 KML。
 *     gnssData：GNSS 原始观测数组。
 *     enableSD：是否启用速度检测辅助。
 *     sdMode：SD 模式位掩码。
 *     outputFile：匹配结果 TXT 输出路径。
 *
 * 输出：
 *     std::vector<CarPoint>：匹配后的车辆状态点数组。
 *
 * 调用位置：
 *     main.cpp 的模式 2：状态机 MM 模式。
 *
 * 数据流作用：
 *     RoadNetwork + CoorData -> MM 输出文件 -> CarPoint 结果。
 */
std::vector<CarPoint> MapMarching(
    const RoadNetwork& road,
    const std::vector<CoorData>& gnssData,
    bool enableSD,
    int sdMode,
    const std::string& outputFile)
{
    std::vector<CoorData> coords = gnssData;
    if (coords.empty()) {
        std::cerr << "MapMarching 输入 GNSS 为空。" << std::endl;
        return {};
    }

    std::string kmlFile = road.sourcePath;
    if (kmlFile.empty()) {
        kmlFile = "out\\mm_input_road.kml";
        if (!WriteRoadNetworkToKml(kmlFile, road)) {
            return {};
        }
    }

    bool oldVirtualLaneOffset = g_MMUseVirtualLaneOffset;
    bool oldEnableSD = g_MMEnableSpeedSD;
    int oldSDMode = g_MMSpeedSDMode;

    g_MMUseVirtualLaneOffset = !road.hasVirtualLanes;
    g_MMEnableSpeedSD = enableSD;
    g_MMSpeedSDMode = sdMode;

    RunStateMachineMapMatching(coords, kmlFile, outputFile);

    g_MMUseVirtualLaneOffset = oldVirtualLaneOffset;
    g_MMEnableSpeedSD = oldEnableSD;
    g_MMSpeedSDMode = oldSDMode;

    std::vector<CarPoint> result = ReadCarTrajectory(outputFile);

    const size_t n = std::min(result.size(), coords.size());
    for (size_t i = 0; i < n; ++i) {
        result[i].V_ENU[0] = coords[i].velE;
        result[i].V_ENU[1] = coords[i].velN;
        result[i].V_ENU[2] = coords[i].velU;
        result[i].trusted = coords[i].trusted;
        result[i].trust_reason = coords[i].trust_reason;
        result[i].quality_score = coords[i].quality_score;
    }

    return result;
}

