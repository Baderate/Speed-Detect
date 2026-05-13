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
#include <unordered_map>
#include <string>
#include <functional>

using namespace std;

// ============================================================
// 传统 HMM 地图匹配（经典基线版 + SD 开关）
//
// 核心：
// 1. 隐状态：每个历元的道路投影候选点
// 2. 发射概率：点到候选点的距离高斯模型
// 3. 转移概率：路网最短路径长度 与 观测直线距离 的差
// 4. 求解：Viterbi
//
// enableSD = false：
//     严格传统 HMM，不读取 trusted / quality_score
//
// enableSD = true：
//     读取 SD.cpp 写入的 trusted / quality_score
//     对异常点：
//        1. 放大发射 sigma，降低 GNSS 对 HMM 的拉力
//        2. 放大转移 beta，放松异常点附近的运动约束

// ============================================================

static const double R_EARTH_HMM = 6378137.0;

// ------------------------- 参数 -------------------------

// 最大候选吸附距离（m）
static const double HMM_MAX_CAND_DIST = 35.0;

// 每个历元最多保留的候选数
static const int HMM_MAX_CAND_NUM = 8;

// 正常 HMM 发射概率 sigma
static const double HMM_EMIT_SIGMA_M = 8.0;

// SD 异常点发射 sigma
// 越大，异常 GNSS 点对道路候选的约束越弱
static const double HMM_EMIT_SIGMA_SD_BAD_M = 22.0;

// quality_score 下限，防止 sigma 无限放大
static const double HMM_SD_MIN_QUALITY = 0.15;

// 低于该质量分数，也视为弱观测
static const double HMM_SD_WEAK_QUALITY_TH = 0.35;

// 正常转移 beta
static const double HMM_BETA = 15.0;

// SD 异常点附近转移 beta
// 越大，异常点附近转移约束越弱
static const double HMM_BETA_SD_BAD = 30.0;

// SD 异常点附近，限制观测步长对转移项的影响
static const double HMM_SD_OBS_STEP_CLAMP_M = 35.0;

// 长时间间断后重新起链（s）
static const double HMM_RESET_GAP_SEC = 5.0;

// KML 端点合并阈值（m）
static const double HMM_NODE_MERGE_THRESH = 3.0;

// 若没有候选，兜底仍保留最近道路
static const bool HMM_ALWAYS_KEEP_NEAREST = true;

// 全局运行时开关，由主函数参数 enableSD 控制
static bool g_HMM_ENABLE_SD = false;


// ------------------------- 工具函数 -------------------------

inline double sqr_hmm(double x)
{
    return x * x;
}

inline double clamp_hmm(double x, double lo, double hi)
{
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

static std::string makeSiblingPath_hmm(
    const std::string& baseFile,
    const std::string& newName
)
{
    size_t p1 = baseFile.find_last_of("/\\");
    if (p1 == std::string::npos) {
        return newName;
    }

    return baseFile.substr(0, p1 + 1) + newName;
}

static std::string makeDebugPathFromOutput_hmm(
    const std::string& out_file,
    const std::string& suffix
)
{
    size_t p1 = out_file.find_last_of("/\\");
    std::string dir = "";
    std::string name = out_file;

    if (p1 != std::string::npos) {
        dir = out_file.substr(0, p1 + 1);
        name = out_file.substr(p1 + 1);
    }

    size_t dot = name.find_last_of('.');
    std::string stem = name;

    if (dot != std::string::npos) {
        stem = name.substr(0, dot);
    }

    return dir + stem + suffix;
}

inline void ll2xy_hmm(
    double lat_deg,
    double lon_deg,
    const LocalRef& ref,
    double& x,
    double& y
)
{
    double lat_rad = lat_deg * PI / 180.0;
    double lon_rad = lon_deg * PI / 180.0;

    double dLat = lat_rad - ref.lat0_rad;
    double dLon = lon_rad - ref.lon0_rad;

    y = R_EARTH_HMM * dLat;
    x = R_EARTH_HMM * dLon * ref.cosLat0;
}

inline void xy2ll_hmm(
    double x,
    double y,
    const LocalRef& ref,
    double& lat_deg,
    double& lon_deg
)
{
    double lat_rad = y / R_EARTH_HMM + ref.lat0_rad;
    double lon_rad = x / (R_EARTH_HMM * ref.cosLat0) + ref.lon0_rad;

    lat_deg = lat_rad * 180.0 / PI;
    lon_deg = lon_rad * 180.0 / PI;
}

static double pointToPointDist(
    double x1,
    double y1,
    double x2,
    double y2
)
{
    return std::sqrt(
        sqr_hmm(x1 - x2) +
        sqr_hmm(y1 - y2)
    );
}

static void projectPointToSegment_hmm(
    double px,
    double py,
    const RoadSegmentXY& seg,
    double& xc,
    double& yc,
    double& tproj,
    double& segLen
)
{
    double vx = seg.x2 - seg.x1;
    double vy = seg.y2 - seg.y1;

    double segLen2 = vx * vx + vy * vy;

    if (segLen2 < 1e-12)
    {
        xc = seg.x1;
        yc = seg.y1;
        tproj = 0.0;
        segLen = 0.0;
        return;
    }

    segLen = std::sqrt(segLen2);

    double wx = px - seg.x1;
    double wy = py - seg.y1;

    tproj = (wx * vx + wy * vy) / segLen2;

    if (tproj < 0.0) {
        tproj = 0.0;
    }
    else if (tproj > 1.0) {
        tproj = 1.0;
    }

    xc = seg.x1 + tproj * vx;
    yc = seg.y1 + tproj * vy;
}


// ------------------------- SD 辅助函数 -------------------------

static double hmmGetQuality_hmm(const CoorData& c)
{
    if (!g_HMM_ENABLE_SD) {
        return 1.0;
    }

    double q = c.quality_score;

    if (!std::isfinite(q)) {
        q = 1.0;
    }

    if (q < HMM_SD_MIN_QUALITY) {
        q = HMM_SD_MIN_QUALITY;
    }

    if (q > 1.0) {
        q = 1.0;
    }

    return q;
}

static bool hmmIsWeakObs_hmm(const CoorData& c)
{
    if (!g_HMM_ENABLE_SD) {
        return false;
    }

    if (!c.trusted) {
        return true;
    }

    return hmmGetQuality_hmm(c) < HMM_SD_WEAK_QUALITY_TH;
}

static double hmmEmissionSigma_hmm(const CoorData& c)
{
    if (!g_HMM_ENABLE_SD) {
        return HMM_EMIT_SIGMA_M;
    }

    double q = hmmGetQuality_hmm(c);

    double sigma = HMM_EMIT_SIGMA_M;

    if (hmmIsWeakObs_hmm(c)) {
        sigma = HMM_EMIT_SIGMA_SD_BAD_M;
    }

    // q 越低，sigma 越大，GNSS 对 HMM 的拉力越弱
    sigma = sigma / std::sqrt(q);

    return sigma;
}

static double hmmEmissionCost_hmm(
    double obsDistM,
    const CoorData& coor
)
{
    if (!std::isfinite(obsDistM)) {
        return std::numeric_limits<double>::infinity();
    }

    if (obsDistM < 0.0) {
        obsDistM = -obsDistM;
    }

    double sigma = hmmEmissionSigma_hmm(coor);

    return 0.5 * sqr_hmm(obsDistM / sigma);
}


// ------------------------- 内部结构 -------------------------

struct GraphEdge
{
    int to = -1;
    double w = 0.0;
};

struct SegExt
{
    RoadSegmentXY seg;
    int u = -1;
    int v = -1;
    double len = 0.0;
};

struct HMMCandidate
{
    int segIdx = -1;

    double x = 0.0;
    double y = 0.0;

    // 候选点在线段上的弧长位置
    double s = 0.0;

    double obsDist = 1e9;
    double emitCost = 1e9;
};

struct ChainRange
{
    int L = 0;
    int R = 0;
};


// ------------------------- KML 解析 -------------------------

static std::vector<RoadSegmentXY> parseRoadSegmentsFromKML_hmm(
    const std::string& kml_file,
    LocalRef& localRef
)
{
    std::vector<RoadSegmentXY> segments;

    std::ifstream file(kml_file.c_str());

    if (!file.is_open())
    {
        std::cerr << "无法打开 KML 文件: "
            << kml_file
            << std::endl;
        return segments;
    }

    std::string content(
        (std::istreambuf_iterator<char>(file)),
        std::istreambuf_iterator<char>()
    );

    file.close();

    std::regex re("<coordinates>([^<]+)</coordinates>");
    std::smatch m;
    std::string::const_iterator it = content.cbegin();

    std::vector<std::pair<double, double>> allPoints;

    while (std::regex_search(it, content.cend(), m, re))
    {
        std::stringstream ss(m[1].str());
        std::string coord;
        std::vector<std::pair<double, double>> pts;

        while (ss >> coord)
        {
            size_t c1 = coord.find(',');
            size_t c2 = coord.find(',', c1 + 1);

            if (c1 == std::string::npos ||
                c2 == std::string::npos) {
                continue;
            }

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

                seg.x1 = 0.0;
                seg.y1 = 0.0;
                seg.x2 = 0.0;
                seg.y2 = 0.0;

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

    double lat0 = sumLat / static_cast<double>(allPoints.size());
    double lon0 = sumLon / static_cast<double>(allPoints.size());

    localRef.lat0_rad = lat0 * PI / 180.0;
    localRef.lon0_rad = lon0 * PI / 180.0;
    localRef.cosLat0 = std::cos(localRef.lat0_rad);

    for (size_t i = 0; i < segments.size(); ++i)
    {
        RoadSegmentXY& seg = segments[i];

        ll2xy_hmm(
            seg.lat1,
            seg.lon1,
            localRef,
            seg.x1,
            seg.y1
        );

        ll2xy_hmm(
            seg.lat2,
            seg.lon2,
            localRef,
            seg.x2,
            seg.y2
        );

        seg.id = static_cast<int>(i);

        double dx = seg.x2 - seg.x1;
        double dy = seg.y2 - seg.y1;

        seg.heading = std::atan2(dy, dx);
    }

    std::cout << "传统 HMM：KML 解析得到道路线段数 = "
        << segments.size()
        << std::endl;

    return segments;
}


// ------------------------- 道路图构建 -------------------------

static int findOrCreateNode_hmm(
    double x,
    double y,
    std::vector<std::pair<double, double>>& nodeXY
)
{
    for (int i = 0; i < static_cast<int>(nodeXY.size()); ++i)
    {
        double dx = x - nodeXY[i].first;
        double dy = y - nodeXY[i].second;

        if (dx * dx + dy * dy <=
            HMM_NODE_MERGE_THRESH * HMM_NODE_MERGE_THRESH) {
            return i;
        }
    }

    nodeXY.push_back(std::make_pair(x, y));

    return static_cast<int>(nodeXY.size()) - 1;
}

static void buildRoadGraph_hmm(
    const std::vector<RoadSegmentXY>& rawSegs,
    std::vector<SegExt>& segs,
    std::vector<std::vector<GraphEdge>>& graph,
    std::vector<std::pair<double, double>>& nodeXY
)
{
    segs.clear();
    nodeXY.clear();

    for (size_t i = 0; i < rawSegs.size(); ++i)
    {
        SegExt e;

        e.seg = rawSegs[i];

        double dx = e.seg.x2 - e.seg.x1;
        double dy = e.seg.y2 - e.seg.y1;

        e.len = std::sqrt(dx * dx + dy * dy);

        e.u = findOrCreateNode_hmm(e.seg.x1, e.seg.y1, nodeXY);
        e.v = findOrCreateNode_hmm(e.seg.x2, e.seg.y2, nodeXY);

        segs.push_back(e);
    }

    graph.assign(nodeXY.size(), std::vector<GraphEdge>());

    for (size_t i = 0; i < segs.size(); ++i)
    {
        const SegExt& e = segs[i];

        if (e.u == e.v || e.len < 1e-8) {
            continue;
        }

        graph[e.u].push_back({ e.v, e.len });
        graph[e.v].push_back({ e.u, e.len });
    }

    std::cout << "传统 HMM：道路图节点数 = "
        << nodeXY.size()
        << "，边数 = "
        << segs.size()
        << std::endl;
}


// ------------------------- 候选生成 -------------------------

static std::vector<HMMCandidate> buildCandidatesForEpoch_hmm(
    const CoorData& coor,
    const std::vector<SegExt>& segs,
    const LocalRef& localRef
)
{
    std::vector<HMMCandidate> cands;

    double obsx = 0.0;
    double obsy = 0.0;

    ll2xy_hmm(
        coor.BLH[0],
        coor.BLH[1],
        localRef,
        obsx,
        obsy
    );

    double bestGlobalDist = std::numeric_limits<double>::infinity();
    HMMCandidate bestGlobalCand;

    for (int i = 0; i < static_cast<int>(segs.size()); ++i)
    {
        const RoadSegmentXY& seg = segs[i].seg;

        double xc = 0.0;
        double yc = 0.0;
        double tproj = 0.0;
        double segLen = 0.0;

        projectPointToSegment_hmm(
            obsx,
            obsy,
            seg,
            xc,
            yc,
            tproj,
            segLen
        );

        double d = pointToPointDist(
            obsx,
            obsy,
            xc,
            yc
        );

        HMMCandidate cand;

        cand.segIdx = i;
        cand.x = xc;
        cand.y = yc;
        cand.s = tproj * segLen;
        cand.obsDist = d;

        // 发射代价：
        // enableSD = false 时，用固定 sigma；
        // enableSD = true 时，根据 trusted / quality_score 调整 sigma。
        cand.emitCost = hmmEmissionCost_hmm(d, coor);

        if (d < bestGlobalDist)
        {
            bestGlobalDist = d;
            bestGlobalCand = cand;
        }

        if (d <= HMM_MAX_CAND_DIST)
        {
            cands.push_back(cand);
        }
    }

    if (cands.empty() && HMM_ALWAYS_KEEP_NEAREST)
    {
        cands.push_back(bestGlobalCand);
    }

    std::sort(
        cands.begin(),
        cands.end(),
        [](const HMMCandidate& a, const HMMCandidate& b)
        {
            return a.obsDist < b.obsDist;
        }
    );

    if (static_cast<int>(cands.size()) > HMM_MAX_CAND_NUM)
    {
        cands.resize(HMM_MAX_CAND_NUM);
    }

    return cands;
}


// ------------------------- 最短路缓存 -------------------------

static std::vector<double> runDijkstra_hmm(
    int src,
    const std::vector<std::vector<GraphEdge>>& graph
)
{
    const double INF = std::numeric_limits<double>::infinity();

    std::vector<double> dist(graph.size(), INF);

    using PDI = std::pair<double, int>;

    std::priority_queue<
        PDI,
        std::vector<PDI>,
        std::greater<PDI>
    > pq;

    dist[src] = 0.0;
    pq.push(std::make_pair(0.0, src));

    while (!pq.empty())
    {
        PDI cur = pq.top();
        pq.pop();

        double d = cur.first;
        int u = cur.second;

        if (d > dist[u]) {
            continue;
        }

        for (size_t k = 0; k < graph[u].size(); ++k)
        {
            int v = graph[u][k].to;
            double nd = d + graph[u][k].w;

            if (nd < dist[v])
            {
                dist[v] = nd;
                pq.push(std::make_pair(nd, v));
            }
        }
    }

    return dist;
}

static double shortestNodeDistance_hmm(
    int src,
    int dst,
    const std::vector<std::vector<GraphEdge>>& graph,
    std::unordered_map<int, std::vector<double>>& cache
)
{
    if (src == dst) {
        return 0.0;
    }

    auto it = cache.find(src);

    if (it == cache.end())
    {
        std::vector<double> d = runDijkstra_hmm(src, graph);
        it = cache.emplace(src, std::move(d)).first;
    }

    const std::vector<double>& dist = it->second;

    if (dst < 0 || dst >= static_cast<int>(dist.size()))
    {
        return std::numeric_limits<double>::infinity();
    }

    return dist[dst];
}


// ------------------------- 路网距离与转移代价 -------------------------

static double routeDistanceBetweenCandidates_hmm(
    const HMMCandidate& a,
    const HMMCandidate& b,
    const std::vector<SegExt>& segs,
    const std::vector<std::vector<GraphEdge>>& graph,
    std::unordered_map<int, std::vector<double>>& cache
)
{
    const SegExt& sa = segs[a.segIdx];
    const SegExt& sb = segs[b.segIdx];

    if (a.segIdx == b.segIdx)
    {
        return std::fabs(b.s - a.s);
    }

    double a_to_u = a.s;
    double a_to_v = sa.len - a.s;

    double b_from_u = b.s;
    double b_from_v = sb.len - b.s;

    double d_uu = shortestNodeDistance_hmm(sa.u, sb.u, graph, cache);
    double d_uv = shortestNodeDistance_hmm(sa.u, sb.v, graph, cache);
    double d_vu = shortestNodeDistance_hmm(sa.v, sb.u, graph, cache);
    double d_vv = shortestNodeDistance_hmm(sa.v, sb.v, graph, cache);

    double ans = std::numeric_limits<double>::infinity();

    if (std::isfinite(d_uu)) {
        ans = std::min(ans, a_to_u + d_uu + b_from_u);
    }

    if (std::isfinite(d_uv)) {
        ans = std::min(ans, a_to_u + d_uv + b_from_v);
    }

    if (std::isfinite(d_vu)) {
        ans = std::min(ans, a_to_v + d_vu + b_from_u);
    }

    if (std::isfinite(d_vv)) {
        ans = std::min(ans, a_to_v + d_vv + b_from_v);
    }

    return ans;
}

static double transitionCost_hmm(
    const HMMCandidate& prevCand,
    const HMMCandidate& curCand,
    double obsDistStep,
    const CoorData& prevObs,
    const CoorData& curObs,
    const std::vector<SegExt>& segs,
    const std::vector<std::vector<GraphEdge>>& graph,
    std::unordered_map<int, std::vector<double>>& cache
)
{
    double routeDist = routeDistanceBetweenCandidates_hmm(
        prevCand,
        curCand,
        segs,
        graph,
        cache
    );

    if (!std::isfinite(routeDist)) {
        return 1e6;
    }

    double beta = HMM_BETA;
    double step = obsDistStep;

    // SD 开启时，如果前后任一点是弱观测，则放松转移约束
    if (g_HMM_ENABLE_SD &&
        (hmmIsWeakObs_hmm(prevObs) || hmmIsWeakObs_hmm(curObs)))
    {
        beta = HMM_BETA_SD_BAD;
        step = std::min(step, HMM_SD_OBS_STEP_CLAMP_M);
    }

    return std::fabs(routeDist - step) / beta;
}


// ------------------------- 分链 -------------------------

static std::vector<ChainRange> splitIntoChains(
    const std::vector<CoorData>& coords
)
{
    std::vector<ChainRange> chains;

    if (coords.empty()) {
        return chains;
    }

    int L = 0;

    for (int i = 1; i < static_cast<int>(coords.size()); ++i)
    {
        long long t1 = ComputeCoorTimestampSec(coords[i - 1]);
        long long t2 = ComputeCoorTimestampSec(coords[i]);

        bool breakChain = false;

        if (t1 < 0 || t2 < 0)
        {
            breakChain = true;
        }
        else
        {
            double dt = static_cast<double>(t2 - t1);

            if (dt > HMM_RESET_GAP_SEC) {
                breakChain = true;
            }

            if (dt < 0.0) {
                breakChain = true;
            }
        }

        if (breakChain)
        {
            chains.push_back({ L, i - 1 });
            L = i;
        }
    }

    chains.push_back({ L, static_cast<int>(coords.size()) - 1 });

    return chains;
}


// ------------------------- 主函数 -------------------------

/**
 * 函数功能：
 *     经典 HMM 地图匹配核心流程。
 *     负责解析路网、为每个历元生成候选点、构建发射/转移代价，并用 Viterbi 求最优路径。
 *
 * 输入：
 *     coords：GNSS 原始观测数组；enableSD=true 时会读取其中的 trusted 和 quality_score。
 *     kml_file：道路或虚拟车道 KML 文件路径。
 *     out_file：匹配结果 TXT 输出路径。
 *     enableSD：是否使用 SD 质量信息调整 HMM 概率模型。
 *
 * 输出：
 *     无返回值；结果写入 out_file，调试信息写入同目录 *_debug.txt。
 *
 * 调用位置：
 *     HiddenMarkovMapMatching 对外入口内部调用。
 *
 * 数据流作用：
 *     CoorData + Road KML -> out\hmm_matched_points.txt + debug 文件。
 */
void RunClassicHmmMapMatching(
    std::vector<CoorData>& coords,
    const std::string& kml_file,
    const std::string& out_file,
    bool enableSD
)
{
    g_HMM_ENABLE_SD = enableSD;

    std::cout << "========== 传统 HMM 地图匹配 ==========" << std::endl;
    std::cout << "传统 HMM：enableSD = "
        << (g_HMM_ENABLE_SD ? "true" : "false")
        << std::endl;

    if (coords.empty())
    {
        std::cerr << "输入轨迹为空，无法进行传统 HMM 地图匹配"
            << std::endl;
        return;
    }

    LocalRef localRef;

    std::vector<RoadSegmentXY> rawSegs =
        parseRoadSegmentsFromKML_hmm(
            kml_file,
            localRef
        );

    if (rawSegs.empty())
    {
        std::cerr << "传统 HMM：道路解析失败" << std::endl;
        return;
    }

    std::vector<SegExt> segs;
    std::vector<std::vector<GraphEdge>> graph;
    std::vector<std::pair<double, double>> nodeXY;

    buildRoadGraph_hmm(
        rawSegs,
        segs,
        graph,
        nodeXY
    );

    int T = static_cast<int>(coords.size());

    std::vector<std::pair<double, double>> obsXY(T);

    for (int t = 0; t < T; ++t)
    {
        double x = 0.0;
        double y = 0.0;

        ll2xy_hmm(
            coords[t].BLH[0],
            coords[t].BLH[1],
            localRef,
            x,
            y
        );

        obsXY[t] = std::make_pair(x, y);
    }

    std::vector<std::vector<HMMCandidate>> allCands(T);

    for (int t = 0; t < T; ++t)
    {
        allCands[t] = buildCandidatesForEpoch_hmm(
            coords[t],
            segs,
            localRef
        );

        if (allCands[t].empty())
        {
            std::cerr << "传统 HMM：历元 "
                << t
                << " 没有候选，终止"
                << std::endl;
            return;
        }
    }

    std::vector<int> bestIdx(T, -1);

    std::vector<ChainRange> chains = splitIntoChains(coords);

    for (size_t c = 0; c < chains.size(); ++c)
    {
        int L = chains[c].L;
        int R = chains[c].R;
        int N = R - L + 1;

        std::vector<std::vector<double>> dp(N);
        std::vector<std::vector<int>> parent(N);

        for (int k = 0; k < N; ++k)
        {
            int t = L + k;
            int nc = static_cast<int>(allCands[t].size());

            dp[k].assign(
                nc,
                std::numeric_limits<double>::infinity()
            );

            parent[k].assign(nc, -1);
        }

        // 初始化
        for (int j = 0; j < static_cast<int>(allCands[L].size()); ++j)
        {
            dp[0][j] = allCands[L][j].emitCost;
        }

        std::unordered_map<int, std::vector<double>> distCache;

        // Viterbi 递推
        for (int k = 1; k < N; ++k)
        {
            int tPrev = L + k - 1;
            int tCur = L + k;

            double obsStep = pointToPointDist(
                obsXY[tPrev].first,
                obsXY[tPrev].second,
                obsXY[tCur].first,
                obsXY[tCur].second
            );

            obsStep = clamp_hmm(obsStep, 0.0, 200.0);

            for (int j = 0; j < static_cast<int>(allCands[tCur].size()); ++j)
            {
                double bestCost =
                    std::numeric_limits<double>::infinity();

                int bestPrev = -1;

                for (int i = 0; i < static_cast<int>(allCands[tPrev].size()); ++i)
                {
                    double trans = transitionCost_hmm(
                        allCands[tPrev][i],
                        allCands[tCur][j],
                        obsStep,
                        coords[tPrev],
                        coords[tCur],
                        segs,
                        graph,
                        distCache
                    );

                    double total =
                        dp[k - 1][i]
                        + trans
                        + allCands[tCur][j].emitCost;

                    if (total < bestCost)
                    {
                        bestCost = total;
                        bestPrev = i;
                    }
                }

                dp[k][j] = bestCost;
                parent[k][j] = bestPrev;
            }
        }

        // 取最后一帧最优状态
        double bestFinal =
            std::numeric_limits<double>::infinity();

        int bestLast = -1;

        for (int j = 0; j < static_cast<int>(allCands[R].size()); ++j)
        {
            if (dp[N - 1][j] < bestFinal)
            {
                bestFinal = dp[N - 1][j];
                bestLast = j;
            }
        }

        if (bestLast < 0)
        {
            std::cerr << "传统 HMM：链回溯失败" << std::endl;
            return;
        }

        bestIdx[R] = bestLast;

        // 回溯
        for (int k = N - 1; k >= 1; --k)
        {
            int t = L + k;
            int tp = L + k - 1;

            bestIdx[tp] = parent[k][bestIdx[t]];

            if (bestIdx[tp] < 0) {
                bestIdx[tp] = 0;
            }
        }
    }

    std::ofstream fout(out_file.c_str());

    if (!fout.is_open())
    {
        std::cerr << "无法打开输出文件: "
            << out_file
            << std::endl;
        return;
    }

    std::string debug_file =
        makeDebugPathFromOutput_hmm(
            out_file,
            "_debug.txt"
        );

    std::ofstream fdbg(debug_file.c_str());

    fout << "timestamp(ms, UTC)\tlatitude(deg)\tlongitude(deg)\theight(m)\n";
    fout << std::fixed << std::setprecision(8);

    if (fdbg.is_open())
    {
        fdbg << "timestamp(ms)\t"
            << "lat\tlon\t"
            << "seg_id\t"
            << "obs_dist(m)\t"
            << "emit_cost\t"
            << "emit_sigma\t"
            << "input_trusted\t"
            << "input_quality\t"
            << "weak_obs\n";

        fdbg << std::fixed << std::setprecision(8);
    }

    for (int t = 0; t < T; ++t)
    {
        if (bestIdx[t] < 0 ||
            bestIdx[t] >= static_cast<int>(allCands[t].size()))
        {
            std::cerr << "传统 HMM：最终索引异常，t = "
                << t
                << std::endl;
            continue;
        }

        const HMMCandidate& best = allCands[t][bestIdx[t]];

        double lat = 0.0;
        double lon = 0.0;

        xy2ll_hmm(
            best.x,
            best.y,
            localRef,
            lat,
            lon
        );

        long long t_utc = ComputeCoorTimestampSec(coords[t]);
        long long t_ms = (t_utc < 0) ? 0LL : t_utc * 1000LL;

        double q = std::exp(-best.emitCost);
        q = clamp_hmm(q, 0.05, 1.0);

        // HMM SD OFF 时，可以写回 HMM 自己的质量信息；
        // HMM SD ON 时，保留 SD.cpp 给的 trusted / quality_score，
        // 避免把 SD 标记覆盖掉。
        if (!g_HMM_ENABLE_SD)
        {
            coords[t].quality_score = q;
            coords[t].trusted = (best.obsDist < HMM_MAX_CAND_DIST);
            coords[t].trust_reason = coords[t].trusted ? 0 : 1;
        }

        fout << t_ms << "\t"
            << lat << "\t"
            << lon << "\t"
            << coords[t].BLH[2] << "\n";

        if (fdbg.is_open())
        {
            double sigma = hmmEmissionSigma_hmm(coords[t]);
            bool weak = hmmIsWeakObs_hmm(coords[t]);

            fdbg << t_ms << "\t"
                << lat << "\t"
                << lon << "\t"
                << segs[best.segIdx].seg.id << "\t"
                << best.obsDist << "\t"
                << best.emitCost << "\t"
                << sigma << "\t"
                << (coords[t].trusted ? 1 : 0) << "\t"
                << coords[t].quality_score << "\t"
                << (weak ? 1 : 0) << "\n";
        }
    }

    fout.close();

    if (fdbg.is_open()) {
        fdbg.close();
    }

    std::cout << "传统 HMM 匹配完成，输出文件: "
        << out_file
        << std::endl;

    std::cout << "传统 HMM 调试文件: "
        << debug_file
        << std::endl;
}

/**
 * 函数功能：
 *     HMM 对外入口，可选先执行 SD 质量检测，再运行 HMM 匹配并返回统一 CarPoint 结果。
 *
 * 输入：
 *     road：路网数据，sourcePath 指向实际用于 HMM 的 KML。
 *     gnssData：GNSS 原始观测数组。
 *     enableSD：是否先运行 SD 并把质量评分传给 HMM。
 *     outputFile：匹配结果 TXT 输出路径。
 *
 * 输出：
 *     std::vector<CarPoint>：HMM 匹配后的车辆状态点数组。
 *
 * 调用位置：
 *     main.cpp 的模式 3：HMM 模式。
 *
 * 数据流作用：
 *     RoadNetwork + CoorData -> HMM 输出文件 -> CarPoint 结果。
 */
std::vector<CarPoint> HiddenMarkovMapMatching(
    const RoadNetwork& road,
    const std::vector<CoorData>& gnssData,
    bool enableSD,
    const std::string& outputFile)
{
    std::vector<CoorData> coords = gnssData;
    if (coords.empty()) {
        std::cerr << "HMM 输入 GNSS 为空。" << std::endl;
        return {};
    }

    std::string kmlFile = road.sourcePath;
    if (kmlFile.empty()) {
        kmlFile = "out\\hmm_input_road.kml";
        if (!WriteRoadNetworkToKml(kmlFile, road)) {
            return {};
        }
    }

    if (enableSD) {
        std::vector<CarPoint> sdPoints = ConvertCoorToCarPoints(coords);
        RunSpeedDetection(sdPoints, SDMode::SpeedPositionPrediction);

        const size_t n = std::min(sdPoints.size(), coords.size());
        for (size_t i = 0; i < n; ++i) {
            coords[i].trusted = sdPoints[i].trusted;
            coords[i].trust_reason = sdPoints[i].trust_reason;
            coords[i].quality_score = sdPoints[i].quality_score;
        }
    }

    RunClassicHmmMapMatching(coords, kmlFile, outputFile, enableSD);

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

