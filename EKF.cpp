#include "Head.h"

namespace
{
constexpr int EKF_STATE_SIZE = 4;

double ClampDouble(double value, double lo, double hi)
{
    if (value < lo) return lo;
    if (value > hi) return hi;
    return value;
}

long long AbsLongLong(long long value)
{
    return value >= 0 ? value : -value;
}

LocalRef BuildLocalRefFromGnss(const std::vector<CoorData>& gnssData)
{
    double sumLat = 0.0;
    double sumLon = 0.0;
    int count = 0;

    for (const auto& point : gnssData) {
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
        if (std::fabs(ref.cosLat0) < 1e-9) {
            ref.cosLat0 = 1.0;
        }
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

void XYToBlh(double x, double y, const LocalRef& ref, double& latDeg, double& lonDeg)
{
    const double latRad = y / earth_R + ref.lat0_rad;
    const double lonRad = x / (earth_R * ref.cosLat0) + ref.lon0_rad;
    latDeg = latRad * 180.0 / PI;
    lonDeg = lonRad * 180.0 / PI;
}

double QualityNoiseScale(bool trusted, double qualityScore, const EKFOptions& options)
{
    double quality = qualityScore;
    if (!std::isfinite(quality)) {
        quality = 1.0;
    }

    quality = ClampDouble(quality, options.minQuality, 1.0);
    double scale = 1.0 / (quality * quality);
    if (!trusted) {
        scale *= options.untrustedNoiseScale;
    }
    return scale;
}

bool HasValidGnssPosition(const CoorData& point)
{
    return std::isfinite(point.BLH[0]) && std::isfinite(point.BLH[1]);
}

bool HasValidMatchedPosition(const CarPoint& point)
{
    return std::isfinite(point.BLH[0]) && std::isfinite(point.BLH[1]);
}

bool HasValidGnssVelocity(const CoorData& point)
{
    return std::isfinite(point.velE) && std::isfinite(point.velN);
}

long long GetGnssTimestampMs(const CoorData& point)
{
    return GpsTimeToUtcTimestampMs(point.week, point.ws);
}

int FindNearestMatchedIndex(
    const std::vector<CarPoint>& matchedPoints,
    long long timestamp,
    size_t& cursor,
    long long maxDiffMs)
{
    if (matchedPoints.empty()) {
        return -1;
    }

    while (cursor + 1 < matchedPoints.size()
        && matchedPoints[cursor + 1].timestamp <= timestamp) {
        ++cursor;
    }

    int best = static_cast<int>(cursor);
    long long bestDiff = AbsLongLong(matchedPoints[cursor].timestamp - timestamp);

    if (cursor + 1 < matchedPoints.size()) {
        long long nextDiff = AbsLongLong(matchedPoints[cursor + 1].timestamp - timestamp);
        if (nextDiff < bestDiff) {
            best = static_cast<int>(cursor + 1);
            bestDiff = nextDiff;
        }
    }

    if (bestDiff > maxDiffMs) {
        return -1;
    }
    return best;
}

struct MatchedObservation
{
    const CarPoint* point = nullptr;
    int index = -1;
    double x = NaNValue();
    double y = NaNValue();
    bool trusted = true;
    double qualityScore = 1.0;
};

double SafeQualityScore(double value)
{
    return std::isfinite(value) ? ClampDouble(value, 0.0, 1.0) : 1.0;
}

bool PointToXY(const CarPoint& point, const LocalRef& ref, double& x, double& y)
{
    if (!HasValidMatchedPosition(point)) {
        return false;
    }

    BlhToXY(point.BLH[0], point.BLH[1], ref, x, y);
    return std::isfinite(x) && std::isfinite(y);
}

bool TryInterpolatedMatchedObservation(
    const std::vector<CarPoint>& matchedPoints,
    long long timestamp,
    size_t cursor,
    const LocalRef& ref,
    MatchedObservation& observation)
{
    if (matchedPoints.empty() || cursor + 1 >= matchedPoints.size()) {
        return false;
    }

    const CarPoint& left = matchedPoints[cursor];
    const CarPoint& right = matchedPoints[cursor + 1];
    if (timestamp < left.timestamp || timestamp > right.timestamp
        || right.timestamp <= left.timestamp) {
        return false;
    }

    double leftX = 0.0;
    double leftY = 0.0;
    double rightX = 0.0;
    double rightY = 0.0;
    if (!PointToXY(left, ref, leftX, leftY) || !PointToXY(right, ref, rightX, rightY)) {
        return false;
    }

    const double ratio =
        static_cast<double>(timestamp - left.timestamp)
        / static_cast<double>(right.timestamp - left.timestamp);
    if (!std::isfinite(ratio) || ratio < 0.0 || ratio > 1.0) {
        return false;
    }

    observation.x = leftX + (rightX - leftX) * ratio;
    observation.y = leftY + (rightY - leftY) * ratio;

    const double leftQuality = SafeQualityScore(left.quality_score);
    const double rightQuality = SafeQualityScore(right.quality_score);
    observation.qualityScore = leftQuality + (rightQuality - leftQuality) * ratio;
    observation.trusted = left.trusted && right.trusted;
    return true;
}

bool FindMatchedObservation(
    const std::vector<CarPoint>& matchedPoints,
    long long timestamp,
    size_t& cursor,
    long long maxDiffMs,
    const LocalRef& ref,
    MatchedObservation& observation)
{
    const int nearestIndex = FindNearestMatchedIndex(
        matchedPoints,
        timestamp,
        cursor,
        maxDiffMs);
    if (nearestIndex < 0) {
        return false;
    }

    observation.index = nearestIndex;
    observation.point = &matchedPoints[nearestIndex];
    observation.trusted = observation.point->trusted;
    observation.qualityScore = SafeQualityScore(observation.point->quality_score);

    if (TryInterpolatedMatchedObservation(
        matchedPoints,
        timestamp,
        cursor,
        ref,
        observation)) {
        return std::isfinite(observation.x) && std::isfinite(observation.y);
    }

    return PointToXY(*observation.point, ref, observation.x, observation.y);
}

void Symmetrize(double p[EKF_STATE_SIZE][EKF_STATE_SIZE])
{
    for (int r = 0; r < EKF_STATE_SIZE; ++r) {
        for (int c = r + 1; c < EKF_STATE_SIZE; ++c) {
            const double v = 0.5 * (p[r][c] + p[c][r]);
            p[r][c] = v;
            p[c][r] = v;
        }
    }
}

void PredictEKF(
    double state[EKF_STATE_SIZE],
    double p[EKF_STATE_SIZE][EKF_STATE_SIZE],
    double dt,
    const EKFOptions& options)
{
    if (!std::isfinite(dt) || dt <= 0.0) {
        return;
    }

    dt = ClampDouble(dt, 0.001, 10.0);

    state[0] += state[2] * dt;
    state[1] += state[3] * dt;

    double f[EKF_STATE_SIZE][EKF_STATE_SIZE] =
    {
        { 1.0, 0.0, dt,  0.0 },
        { 0.0, 1.0, 0.0, dt  },
        { 0.0, 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 0.0, 1.0 }
    };

    double fp[EKF_STATE_SIZE][EKF_STATE_SIZE] = {};
    for (int r = 0; r < EKF_STATE_SIZE; ++r) {
        for (int c = 0; c < EKF_STATE_SIZE; ++c) {
            for (int k = 0; k < EKF_STATE_SIZE; ++k) {
                fp[r][c] += f[r][k] * p[k][c];
            }
        }
    }

    double predicted[EKF_STATE_SIZE][EKF_STATE_SIZE] = {};
    for (int r = 0; r < EKF_STATE_SIZE; ++r) {
        for (int c = 0; c < EKF_STATE_SIZE; ++c) {
            for (int k = 0; k < EKF_STATE_SIZE; ++k) {
                predicted[r][c] += fp[r][k] * f[c][k];
            }
        }
    }

    const double qPos = options.processPositionNoiseM * options.processPositionNoiseM;
    const double qVel = options.processVelocityNoiseMps * options.processVelocityNoiseMps;

    predicted[0][0] += qPos * dt * dt;
    predicted[1][1] += qPos * dt * dt;
    predicted[2][2] += qVel * dt;
    predicted[3][3] += qVel * dt;

    for (int r = 0; r < EKF_STATE_SIZE; ++r) {
        for (int c = 0; c < EKF_STATE_SIZE; ++c) {
            p[r][c] = predicted[r][c];
        }
    }
    Symmetrize(p);
}

void UpdateDirectState(
    double state[EKF_STATE_SIZE],
    double p[EKF_STATE_SIZE][EKF_STATE_SIZE],
    int stateIndex,
    double observation,
    double variance)
{
    if (!std::isfinite(observation) || !std::isfinite(variance) || variance <= 0.0) {
        return;
    }

    const double innovation = observation - state[stateIndex];
    const double s = p[stateIndex][stateIndex] + variance;
    if (!std::isfinite(s) || s <= 1e-12) {
        return;
    }

    double k[EKF_STATE_SIZE] = {};
    for (int i = 0; i < EKF_STATE_SIZE; ++i) {
        k[i] = p[i][stateIndex] / s;
    }

    double oldP[EKF_STATE_SIZE][EKF_STATE_SIZE] = {};
    for (int r = 0; r < EKF_STATE_SIZE; ++r) {
        for (int c = 0; c < EKF_STATE_SIZE; ++c) {
            oldP[r][c] = p[r][c];
        }
    }

    for (int i = 0; i < EKF_STATE_SIZE; ++i) {
        state[i] += k[i] * innovation;
    }

    for (int r = 0; r < EKF_STATE_SIZE; ++r) {
        for (int c = 0; c < EKF_STATE_SIZE; ++c) {
            p[r][c] = oldP[r][c] - k[r] * oldP[stateIndex][c];
        }
    }
    Symmetrize(p);
}

CarPoint BuildOutputPoint(
    const CoorData& gnss,
    const CarPoint* matched,
    const double state[EKF_STATE_SIZE],
    const LocalRef& ref,
    bool usedRoadObservation)
{
    CarPoint out;
    out.timestamp = GetGnssTimestampMs(gnss);
    out.x = state[0];
    out.y = state[1];
    out.V_ENU[0] = state[2];
    out.V_ENU[1] = state[3];
    out.V_ENU[2] = gnss.velU;

    XYToBlh(out.x, out.y, ref, out.BLH[0], out.BLH[1]);
    out.BLH[2] = gnss.BLH[2];
    if (!std::isfinite(out.BLH[2]) && matched != nullptr) {
        out.BLH[2] = matched->BLH[2];
    }

    double gnssQuality = std::isfinite(gnss.quality_score) ? gnss.quality_score : 1.0;
    gnssQuality = ClampDouble(gnssQuality, 0.0, 1.0);
    double fusedQuality = gnssQuality;

    if (matched != nullptr && usedRoadObservation) {
        double roadQuality = std::isfinite(matched->quality_score) ? matched->quality_score : 1.0;
        roadQuality = ClampDouble(roadQuality, 0.0, 1.0);
        fusedQuality = ClampDouble(0.6 * gnssQuality + 0.4 * roadQuality, 0.0, 1.0);
    }

    out.quality_score = fusedQuality;
    out.trusted = fusedQuality >= 0.35 && (gnss.trusted || usedRoadObservation);
    out.trust_reason = out.trusted ? 0 : 40;
    return out;
}
}

/**
 * 函数功能：
 *     EKF 融合主流程。
 *     状态量为 [x, y, velE, velN]，GNSS 位置/速度作为主观测，
 *     MM 或 HMM 结果作为道路伪观测，并根据 SD 的 trusted / quality_score 自适应调整观测噪声。
 *
 * 输入：
 *     gnssData：GNSS 原始观测序列。
 *     matchedPoints：MM 或 HMM 输出的车辆状态点序列。
 *     options：EKF 参数配置。
 *
 * 输出：
 *     std::vector<CarPoint>：融合后的车辆状态点序列。
 *
 * 调用位置：
 *     main.cpp 的模式 4：EKF 融合。
 *
 * 数据流作用：
 *     CoorData + MM/HMM CarPoint -> EKF CarPoint。
 */
std::vector<CarPoint> RunEKFFusion(
    const std::vector<CoorData>& gnssData,
    const std::vector<CarPoint>& matchedPoints,
    const EKFOptions& options)
{
    std::vector<CarPoint> result;
    if (gnssData.empty()) {
        std::cerr << "EKF 输入 GNSS 为空。" << std::endl;
        return result;
    }

    LocalRef ref = BuildLocalRefFromGnss(gnssData);
    result.reserve(gnssData.size());

    double state[EKF_STATE_SIZE] = {};
    double p[EKF_STATE_SIZE][EKF_STATE_SIZE] = {};
    bool initialized = false;
    long long prevTimestamp = 0;
    size_t matchedCursor = 0;

    for (const auto& gnss : gnssData) {
        if (!HasValidGnssPosition(gnss)) {
            continue;
        }

        const long long timestamp = GetGnssTimestampMs(gnss);
        double gnssX = 0.0;
        double gnssY = 0.0;
        BlhToXY(gnss.BLH[0], gnss.BLH[1], ref, gnssX, gnssY);

        MatchedObservation matchedObservation;
        const bool hasMatchedObservation = FindMatchedObservation(
            matchedPoints,
            timestamp,
            matchedCursor,
            options.maxPseudoTimeDiffMs,
            ref,
            matchedObservation);
        const CarPoint* matched = nullptr;
        double roadX = NaNValue();
        double roadY = NaNValue();
        if (hasMatchedObservation) {
            matched = matchedObservation.point;
            roadX = matchedObservation.x;
            roadY = matchedObservation.y;
        }

        if (!initialized) {
            state[0] = gnssX;
            state[1] = gnssY;
            state[2] = std::isfinite(gnss.velE) ? gnss.velE : 0.0;
            state[3] = std::isfinite(gnss.velN) ? gnss.velN : 0.0;

            p[0][0] = options.initialPositionStdM * options.initialPositionStdM;
            p[1][1] = options.initialPositionStdM * options.initialPositionStdM;
            p[2][2] = options.initialVelocityStdMps * options.initialVelocityStdMps;
            p[3][3] = options.initialVelocityStdMps * options.initialVelocityStdMps;

            initialized = true;
            prevTimestamp = timestamp;
        }
        else {
            double dt = (timestamp - prevTimestamp) / 1000.0;
            if (!std::isfinite(dt) || dt <= 0.0 || dt > 10.0) {
                dt = 1.0;
            }
            PredictEKF(state, p, dt, options);
            prevTimestamp = timestamp;
        }

        const double gnssScale = QualityNoiseScale(
            gnss.trusted,
            gnss.quality_score,
            options);
        const double gnssPosVar = options.gnssPositionNoiseM * options.gnssPositionNoiseM * gnssScale;

        UpdateDirectState(state, p, 0, gnssX, gnssPosVar);
        UpdateDirectState(state, p, 1, gnssY, gnssPosVar);

        if (options.useGnssVelocity && HasValidGnssVelocity(gnss)) {
            const double gnssVelVar =
                options.gnssVelocityNoiseMps * options.gnssVelocityNoiseMps * gnssScale;
            UpdateDirectState(state, p, 2, gnss.velE, gnssVelVar);
            UpdateDirectState(state, p, 3, gnss.velN, gnssVelVar);
        }

        bool usedRoadObservation = false;
        if (options.useRoadPseudo && matched != nullptr
            && std::isfinite(roadX) && std::isfinite(roadY)) {
            const double roadScale = QualityNoiseScale(
                matchedObservation.trusted,
                matchedObservation.qualityScore,
                options);
            const double roadVar =
                options.roadPositionNoiseM * options.roadPositionNoiseM * roadScale;

            UpdateDirectState(state, p, 0, roadX, roadVar);
            UpdateDirectState(state, p, 1, roadY, roadVar);
            usedRoadObservation = true;
        }

        result.push_back(BuildOutputPoint(gnss, matched, state, ref, usedRoadObservation));
    }

    std::cout << "EKF 融合完成：输入 GNSS 点数 = " << gnssData.size()
        << "，输出融合点数 = " << result.size()
        << "，道路伪观测点数 = " << matchedPoints.size()
        << std::endl;

    return result;
}
