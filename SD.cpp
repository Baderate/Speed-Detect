#include "Head.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>
#include <algorithm>
#include <string>


// =====================================================
// MM/HMM 内部使用的 SD 辅助逻辑
// =====================================================
// SD_MM_MODE_* 是状态机地图匹配使用的位掩码。
//   0  : SD 不影响 MM。
//   1  : SD 预测结果可调整搜索中心。
//   2  : SD 预测残差可调整候选点发射代价。
//   4  : SD 预测步长可进入维特比转移代价。
//   8  : SD 速度方向可约束候选道路方向。
//   15 : 启用完整 SD 辅助。
static const int SD_MM_MODE_OFF = 0;
static const int SD_MM_MODE_CENTER = 1 << 0;
static const int SD_MM_MODE_COST = 1 << 1;
static const int SD_MM_MODE_TRANS = 1 << 2;
static const int SD_MM_MODE_HEADING = 1 << 3;
static const int SD_MM_MODE_FULL =
SD_MM_MODE_CENTER | SD_MM_MODE_COST | SD_MM_MODE_TRANS | SD_MM_MODE_HEADING;

static const double SD_MM_SPEED_STATIC_TH = 0.35;
static const double SD_MM_SPEED_HEADING_VALID_TH = 0.80;
static const double SD_MM_CENTER_BLEND_MAX = 0.30;
static const double SD_MM_STOP_CENTER_BLEND = 0.85;
static const double SD_MM_QUALITY_MIN_SOFT = 0.20;
static const double SD_MM_QUALITY_MIN_HARD = 0.08;
static const int SD_MM_SIGMA_WARMUP_MIN = 8;
static const double SD_MM_SIGMA_FLOOR_M = 1.0;
static const double SD_MM_DEFAULT_SOFT_TH = 5.0;
static const double SD_MM_DEFAULT_HARD_TH = 12.0;
static const double SD_MM_DEFAULT_EXTREME_TH = 20.0;
static const double SD_MM_TRANS_COST_BASE_GAIN = 0.45;
static const double SD_MM_TRANS_COST_RESIDUAL_GAIN = 0.35;

struct SDMMResidualStats
{
    int count = 0;
    double mean = 0.0;
    double m2 = 0.0;
};

struct SDMMResidualJudge
{
    double mean = 0.0;
    double sigma = 0.0;
    double soft_th = SD_MM_DEFAULT_SOFT_TH;
    double hard_th = SD_MM_DEFAULT_HARD_TH;
    double extreme_th = SD_MM_DEFAULT_EXTREME_TH;
    bool use_sigma = false;
};

static SDMMResidualStats g_sd_mm_residual_stats;

static double sd_mm_clamp(double x, double lo, double hi)
{
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

static double sd_mm_smoothstep01(double x)
{
    x = sd_mm_clamp(x, 0.0, 1.0);
    return x * x * (3.0 - 2.0 * x);
}

static void SD_MMResetResidualStats()
{
    g_sd_mm_residual_stats = SDMMResidualStats();
}

static void SD_MMUpdateResidualStats(double residual_m)
{
    if (!std::isfinite(residual_m) || residual_m < 0.0)
        return;

    g_sd_mm_residual_stats.count++;

    double delta = residual_m - g_sd_mm_residual_stats.mean;
    g_sd_mm_residual_stats.mean += delta / g_sd_mm_residual_stats.count;
    double delta2 = residual_m - g_sd_mm_residual_stats.mean;
    g_sd_mm_residual_stats.m2 += delta * delta2;
}

static SDMMResidualJudge SD_MMBuildResidualJudge()
{
    SDMMResidualJudge judge;

    if (g_sd_mm_residual_stats.count < SD_MM_SIGMA_WARMUP_MIN)
        return judge;

    double variance = 0.0;
    if (g_sd_mm_residual_stats.count > 1)
        variance = g_sd_mm_residual_stats.m2 / (g_sd_mm_residual_stats.count - 1);

    variance = std::max(variance, 0.0);
    double sigma = std::max(std::sqrt(variance), SD_MM_SIGMA_FLOOR_M);

    judge.mean = g_sd_mm_residual_stats.mean;
    judge.sigma = sigma;
    judge.soft_th = judge.mean + 2.0 * sigma;
    judge.hard_th = judge.mean + 3.0 * sigma;
    judge.extreme_th = judge.mean + 4.0 * sigma;
    judge.use_sigma = true;
    return judge;
}

static bool SD_MMShouldLearnResidual(double residual_m, const SDMMResidualJudge& judge)
{
    if (!std::isfinite(residual_m) || residual_m < 0.0)
        return false;

    if (!judge.use_sigma)
        return residual_m <= SD_MM_DEFAULT_EXTREME_TH;

    return residual_m <= judge.extreme_th;
}

static bool sd_mm_has_mode(int mode, int bit)
{
    return (mode & bit) != 0;
}

static bool SD_MMIsEnabled(bool use_sd, int mode)
{
    return use_sd && mode != SD_MM_MODE_OFF;
}

static bool SD_MMUseCenter(bool use_sd, int mode)
{
    return SD_MMIsEnabled(use_sd, mode) && sd_mm_has_mode(mode, SD_MM_MODE_CENTER);
}

static bool SD_MMUseCost(bool use_sd, int mode)
{
    return SD_MMIsEnabled(use_sd, mode) && sd_mm_has_mode(mode, SD_MM_MODE_COST);
}

static bool SD_MMUseTrans(bool use_sd, int mode)
{
    return SD_MMIsEnabled(use_sd, mode) && sd_mm_has_mode(mode, SD_MM_MODE_TRANS);
}

static bool SD_MMUseHeading(bool use_sd, int mode)
{
    return SD_MMIsEnabled(use_sd, mode) && sd_mm_has_mode(mode, SD_MM_MODE_HEADING);
}

static double SD_MMResidualToQuality(double r, double soft_th, double hard_th, double q_floor)
{
    if (r <= soft_th) return 1.0;
    if (r >= hard_th) return q_floor;

    double alpha = (r - soft_th) / std::max(hard_th - soft_th, 1e-6);
    alpha = sd_mm_smoothstep01(alpha);
    double q = 1.0 - alpha * (1.0 - q_floor);
    if (q < q_floor) q = q_floor;
    if (q > 1.0) q = 1.0;
    return q;
}

static double SD_MMResidualSeverity(double residual_m, double soft_th, double extreme_th)
{
    if (residual_m <= soft_th) return 0.0;
    if (residual_m >= extreme_th) return 1.0;

    double alpha = (residual_m - soft_th) /
        std::max(extreme_th - soft_th, 1e-6);
    return sd_mm_smoothstep01(alpha);
}

static bool SD_MMGetValidSpeedEN(const CoorData& c,
    double& vE, double& vN, double& speed)
{
    vE = 0.0;
    vN = 0.0;
    speed = 0.0;

    if (std::isfinite(c.phaE) && std::isfinite(c.phaN))
    {
        double s = std::sqrt(c.phaE * c.phaE + c.phaN * c.phaN);
        if (std::isfinite(s) && s >= 0.05 && s <= 60.0)
        {
            vE = c.phaE;
            vN = c.phaN;
            speed = s;
            return true;
        }
    }

    if (std::isfinite(c.velE) && std::isfinite(c.velN))
    {
        double s = std::sqrt(c.velE * c.velE + c.velN * c.velN);
        if (std::isfinite(s) && s >= 0.05 && s <= 60.0)
        {
            vE = c.velE;
            vN = c.velN;
            speed = s;
            return true;
        }
    }

    return false;
}

static bool SD_MMIsStaticEpoch(bool use_sd, int mode,
    bool has_prev, bool has_speed_meas, double speed_meas)
{
    return SD_MMIsEnabled(use_sd, mode) &&
        has_prev && has_speed_meas && speed_meas <= SD_MM_SPEED_STATIC_TH;
}

static bool SD_MMUseVelocityHeadingForEpoch(bool use_sd, int mode,
    bool has_speed_meas, double speed_meas)
{
    return SD_MMUseHeading(use_sd, mode) &&
        has_speed_meas && speed_meas >= SD_MM_SPEED_HEADING_VALID_TH;
}

static void SD_MMPredictXY(
    bool is_static_epoch,
    bool has_speed_meas,
    bool has_vel,
    double prev_x,
    double prev_y,
    double vx_last,
    double vy_last,
    double vE_meas,
    double vN_meas,
    double dt,
    double x_obs,
    double y_obs,
    double& x_dr,
    double& y_dr)
{
    x_dr = x_obs;
    y_dr = y_obs;

    if (is_static_epoch)
    {
        x_dr = prev_x;
        y_dr = prev_y;
    }
    else if (has_speed_meas)
    {
        double vE_pred = vE_meas;
        double vN_pred = vN_meas;

        if (has_vel)
        {
            vE_pred = 0.65 * vE_meas + 0.35 * vx_last;
            vN_pred = 0.65 * vN_meas + 0.35 * vy_last;
        }

        x_dr = prev_x + vE_pred * dt;
        y_dr = prev_y + vN_pred * dt;
    }
    else if (has_vel)
    {
        x_dr = prev_x + vx_last * dt;
        y_dr = prev_y + vy_last * dt;
    }
    else
    {
        x_dr = prev_x;
        y_dr = prev_y;
    }
}

static void SD_MMApplyCenterByResidual(
    bool is_static_epoch,
    double x_obs,
    double y_obs,
    double x_dr,
    double y_dr,
    double prev_x,
    double prev_y,
    bool use_sigma,
    double soft_th,
    double hard_th,
    double extreme_th,
    double& x_center,
    double& y_center,
    bool& trusted,
    int& trust_reason,
    double& quality_score)
{
    double pred_res = std::sqrt((x_obs - x_dr) * (x_obs - x_dr) +
        (y_obs - y_dr) * (y_obs - y_dr));
    double q_speed = SD_MMResidualToQuality(
        pred_res,
        soft_th,
        hard_th,
        SD_MM_QUALITY_MIN_HARD);
    double severity = SD_MMResidualSeverity(pred_res, soft_th, extreme_th);

    if (is_static_epoch)
    {
        double alpha = SD_MM_STOP_CENTER_BLEND;
        x_center = (1.0 - alpha) * x_obs + alpha * prev_x;
        y_center = (1.0 - alpha) * y_obs + alpha * prev_y;

        trusted = true;
        trust_reason = 6;
        quality_score = q_speed;
    }
    else
    {
        double alpha = SD_MM_CENTER_BLEND_MAX * severity;
        x_center = (1.0 - alpha) * x_obs + alpha * x_dr;
        y_center = (1.0 - alpha) * y_obs + alpha * y_dr;

        quality_score = q_speed;

        if (pred_res <= soft_th)
        {
            trusted = true;
            trust_reason = use_sigma ? 10 : 0;
            quality_score = 1.0;
        }
        else if (pred_res <= hard_th)
        {
            trusted = true;
            trust_reason = use_sigma ? 11 : 1;
        }
        else if (pred_res <= extreme_th)
        {
            trusted = true;
            trust_reason = use_sigma ? 12 : 2;
        }
        else
        {
            trusted = false;
            trust_reason = use_sigma ? 13 : 3;
            quality_score = SD_MM_QUALITY_MIN_HARD;
        }
    }
}

static void SD_MMAdjustCandidateWeights(
    bool has_prev,
    double x_obs,
    double y_obs,
    double x_dr,
    double y_dr,
    double soft_th,
    double hard_th,
    double extreme_th,
    double& W_OBS,
    double& W_DR,
    double& W_CONT)
{
    double pred_res = has_prev
        ? std::sqrt((x_obs - x_dr) * (x_obs - x_dr) + (y_obs - y_dr) * (y_obs - y_dr))
        : 0.0;

    double q_eff = SD_MMResidualToQuality(
        pred_res,
        soft_th,
        hard_th,
        SD_MM_QUALITY_MIN_SOFT);
    double severity = SD_MMResidualSeverity(pred_res, soft_th, extreme_th);

    W_OBS *= std::max(0.30, 0.15 + 0.85 * q_eff);
    W_DR *= (1.0 + 0.80 * severity + 0.60 * (1.0 - q_eff));
    W_CONT *= (1.0 + 0.45 * severity + 0.45 * (1.0 - q_eff));
}

static double sd_mm_angle_diff(double a, double b)
{
    double d = a - b;
    while (d > PI) d -= 2.0 * PI;
    while (d < -PI) d += 2.0 * PI;
    return std::fabs(d);
}

static double sd_mm_road_heading_diff_bidirectional(double heading_a, double heading_b)
{
    double d1 = sd_mm_angle_diff(heading_a, heading_b);
    double d2 = sd_mm_angle_diff(heading_a, heading_b + PI);
    return std::min(d1, d2);
}

static double SD_MMVelocityHeadingCost(
    bool use_vel_heading,
    bool is_static_epoch,
    double vE_meas,
    double vN_meas,
    double road_heading,
    double W_VEL_DIR)
{
    if (!use_vel_heading || is_static_epoch)
        return 0.0;

    double vel_heading = std::atan2(vN_meas, vE_meas);
    double vel_dir_diff = sd_mm_road_heading_diff_bidirectional(vel_heading, road_heading);
    return (1.0 - std::cos(vel_dir_diff)) * W_VEL_DIR;
}

static double SD_MMPredictedStep(
    bool use_trans,
    bool has_speed_meas,
    double speed_meas,
    double dt,
    double fallback_step)
{
    if (use_trans && has_speed_meas)
        return speed_meas * std::max(dt, 0.2);
    return fallback_step;
}

static double SD_MMTransitionPredictionCost(
    double routeDist,
    double pred_step,
    double trans_pred_beta)
{
    double residual = std::fabs(routeDist - pred_step);
    double severity = SD_MMResidualSeverity(residual, SD_MM_DEFAULT_SOFT_TH, SD_MM_DEFAULT_EXTREME_TH);
    double gain = SD_MM_TRANS_COST_BASE_GAIN + SD_MM_TRANS_COST_RESIDUAL_GAIN * (1.0 - severity);
    return gain * residual / std::max(trans_pred_beta, 1e-6);
}

static void SD_MMUpdateSmoothedVelocity(
    bool has_speed_meas,
    bool has_prev,
    double vE_meas,
    double vN_meas,
    double best_x,
    double best_y,
    double old_prev_x,
    double old_prev_y,
    double dt,
    double& vx_last,
    double& vy_last,
    bool& has_vel)
{
    if (has_speed_meas)
    {
        if (!has_vel)
        {
            vx_last = vE_meas;
            vy_last = vN_meas;
            has_vel = true;
        }
        else
        {
            const double alpha = 0.7;
            vx_last = alpha * vE_meas + (1.0 - alpha) * vx_last;
            vy_last = alpha * vN_meas + (1.0 - alpha) * vy_last;
        }
    }
    else if (has_prev)
    {
        double dt_eff = std::max(dt, 0.2);
        vx_last = (best_x - old_prev_x) / dt_eff;
        vy_last = (best_y - old_prev_y) / dt_eff;
        has_vel = true;
    }
}


/**
 * 函数功能：
 *     为 MM/HMM 提供速度辅助评价，根据不同阶段返回搜索中心、候选代价或转移代价修正。
 *
 * 输入：
 *     c：当前 GNSS 点；只有历元阶段需要，其他阶段可传 nullptr。
 *     req：SD 请求结构，包含阶段、速度、候选、上一帧状态和权重等上下文。
 *     out：SD 输出结构，由函数写入速度、质量评分、候选权重修正和转移代价。
 *
 * 输出：
 *     bool：SD 启用并执行评价时返回 true；SD 关闭时返回 false。
 *
 * 调用位置：
 *     1. MM.cpp 的历元阶段修正搜索中心。
 *     2. MM.cpp 的候选阶段修正发射代价和速度方向代价。
 *     3. MM.cpp 的转移阶段修正 Viterbi 转移代价。
 *     4. MM.cpp 的更新阶段平滑速度状态。
 *
 * 数据流作用：
 *     CoorData + MM 上下文 -> SDMMResult -> MM/HMM 质量辅助。
 */
bool SD_MMEvaluate(const CoorData* c, const SDMMRequest& req, SDMMResult& out)
{
    out = SDMMResult();

    out.enabled = SD_MMIsEnabled(req.use_sd, req.mode);
    out.x_dr = req.x_obs;
    out.y_dr = req.y_obs;
    out.x_center = req.x_obs;
    out.y_center = req.y_obs;
    out.W_OBS = req.W_OBS;
    out.W_DR = req.W_DR;
    out.W_CONT = req.W_CONT;
    out.W_VEL_DIR = req.W_VEL_DIR;
    out.pred_step = req.has_pred_step ? req.pred_step : req.fallback_step;
    out.vx_last = req.vx_last;
    out.vy_last = req.vy_last;
    out.has_vel = req.has_vel;

    if (!out.enabled)
    {
        if (req.stage == SD_MM_STAGE_CANDIDATE)
        {
            // 保持纯 MM 行为：SD 关闭时，不让航位推算项和速度方向项进入发射代价。
            out.W_DR = 0.0;
            out.W_VEL_DIR = 0.0;
        }
        return false;
    }

    if (req.stage == SD_MM_STAGE_EPOCH)
    {
        if (!req.has_prev)
            SD_MMResetResidualStats();

        if (c != nullptr)
            out.has_speed = SD_MMGetValidSpeedEN(*c, out.vE, out.vN, out.speed);

        out.is_static_epoch = SD_MMIsStaticEpoch(
            req.use_sd, req.mode, req.has_prev, out.has_speed, out.speed);

        out.use_vel_heading = SD_MMUseVelocityHeadingForEpoch(
            req.use_sd, req.mode, out.has_speed, out.speed);

        if (req.has_prev)
        {
            SD_MMPredictXY(
                out.is_static_epoch, out.has_speed, req.has_vel,
                req.prev_x, req.prev_y, req.vx_last, req.vy_last,
                out.vE, out.vN, req.dt,
                req.x_obs, req.y_obs,
                out.x_dr, out.y_dr);
        }

        SDMMResidualJudge judge = SD_MMBuildResidualJudge();
        double pred_res = std::sqrt(
            (req.x_obs - out.x_dr) * (req.x_obs - out.x_dr) +
            (req.y_obs - out.y_dr) * (req.y_obs - out.y_dr));

        if (SD_MMUseCenter(req.use_sd, req.mode) && req.has_prev)
        {
            SD_MMApplyCenterByResidual(
                out.is_static_epoch,
                req.x_obs, req.y_obs,
                out.x_dr, out.y_dr,
                req.prev_x, req.prev_y,
                judge.use_sigma,
                judge.soft_th,
                judge.hard_th,
                judge.extreme_th,
                out.x_center, out.y_center,
                out.trusted, out.trust_reason, out.quality_score);
        }
        else
        {
            out.x_center = req.x_obs;
            out.y_center = req.y_obs;
            out.trusted = true;
            out.trust_reason = 0;
            out.quality_score = 1.0;
        }

        if (req.has_prev && SD_MMShouldLearnResidual(pred_res, judge))
            SD_MMUpdateResidualStats(pred_res);

        return true;
    }

    if (req.stage == SD_MM_STAGE_CANDIDATE)
    {
        if (!SD_MMUseHeading(req.use_sd, req.mode))
            out.W_VEL_DIR = 0.0;

        if (SD_MMUseCost(req.use_sd, req.mode))
        {
            SDMMResidualJudge judge = SD_MMBuildResidualJudge();
            SD_MMAdjustCandidateWeights(
                req.has_prev,
                req.x_obs, req.y_obs,
                req.x_dr, req.y_dr,
                judge.soft_th,
                judge.hard_th,
                judge.extreme_th,
                out.W_OBS, out.W_DR, out.W_CONT);
        }

        out.vel_dir_cost = SD_MMVelocityHeadingCost(
            SD_MMUseHeading(req.use_sd, req.mode) && req.use_vel_heading,
            req.is_static_epoch,
            req.vE_meas,
            req.vN_meas,
            req.road_heading,
            out.W_VEL_DIR);

        return true;
    }

    if (req.stage == SD_MM_STAGE_TRANSITION)
    {
        if (req.has_pred_step)
        {
            out.pred_step = req.pred_step;
        }
        else
        {
            out.pred_step = SD_MMPredictedStep(
                SD_MMUseTrans(req.use_sd, req.mode),
                req.has_speed_meas,
                req.speed_meas,
                req.dt,
                req.fallback_step);
        }

        if (SD_MMUseTrans(req.use_sd, req.mode) && std::isfinite(req.route_dist) && req.route_dist >= 0.0)
        {
            out.trans_pred_cost = SD_MMTransitionPredictionCost(
                req.route_dist,
                out.pred_step,
                req.trans_pred_beta);
        }

        return true;
    }

    if (req.stage == SD_MM_STAGE_UPDATE)
    {
        SD_MMUpdateSmoothedVelocity(
            req.has_speed_meas,
            req.has_prev,
            req.vE_meas,
            req.vN_meas,
            req.best_x,
            req.best_y,
            req.old_prev_x,
            req.old_prev_y,
            req.dt,
            out.vx_last,
            out.vy_last,
            out.has_vel);

        return true;
    }

    return true;
}

/**
 * 函数功能：
 *     兼容旧接口的 SD 验证流程。实际验证逻辑统一转入 EvaluateSpeedDetection，
 *     保证所有 SD 验证都先经过 RunSpeedDetection 主流程。
 */
void RunSDModule_RecallOnly(
    const std::vector<CoorData>& coords,
    const std::vector<TxtData>& truthdata,
    double sd_threshold_m,
    double truth_abnormal_threshold_m,
    const std::string& debug_output_file)
{
    (void)sd_threshold_m;
    EvaluateSpeedDetection(
        coords,
        truthdata,
        SDMode::SpeedPositionPrediction,
        truth_abnormal_threshold_m,
        debug_output_file);
}

// ======================== 独立 SD 主流程辅助函数 ========================
static const char* GetSDModeName(SDMode mode)
{
    switch (mode) {
    case SDMode::Off:
        return "关闭";
    case SDMode::ThreeSigma:
        return "3sigma统计检测";
    case SDMode::SpeedPositionPrediction:
        return "速度位置预测检测";
    default:
        return "未知模式";
    }
}

static void ResetSpeedDetectionResult(std::vector<CarPoint>& points)
{
    for (auto& p : points) {
        p.trusted = true;
        p.trust_reason = 0;
        p.quality_score = 1.0;
    }
}

static bool GetSDDeltaTimeSec(const CarPoint& prev, const CarPoint& curr, double& dt)
{
    dt = (curr.timestamp - prev.timestamp) / 1000.0;
    return std::isfinite(dt) && dt > 0.05 && dt <= 10.0;
}

static bool GetSDPointDistanceM(const CarPoint& prev, const CarPoint& curr, double& dist)
{
    if (std::isfinite(prev.x) && std::isfinite(prev.y)
        && std::isfinite(curr.x) && std::isfinite(curr.y)) {
        const double dx = curr.x - prev.x;
        const double dy = curr.y - prev.y;
        dist = std::sqrt(dx * dx + dy * dy);
        return true;
    }

    if (std::isfinite(prev.BLH[0]) && std::isfinite(prev.BLH[1])
        && std::isfinite(curr.BLH[0]) && std::isfinite(curr.BLH[1])) {
        dist = haversine(prev.BLH[0], prev.BLH[1], curr.BLH[0], curr.BLH[1]);
        return true;
    }

    dist = NaNValue();
    return false;
}

static double BuildSDQualityByResidual(double residual, double softThreshold, double hardThreshold)
{
    if (!std::isfinite(residual) || residual <= softThreshold) {
        return 1.0;
    }

    const double span = std::max(hardThreshold - softThreshold, 1e-6);
    const double t = sd_mm_clamp((residual - softThreshold) / span, 0.0, 1.0);
    return 1.0 - 0.85 * t;
}

static void RunThreeSigmaSpeedDetection(std::vector<CarPoint>& points)
{
    ResetSpeedDetectionResult(points);
    if (points.size() < 3) {
        return;
    }

    std::vector<double> speeds(points.size(), NaNValue());
    double sum = 0.0;
    int count = 0;

    for (size_t i = 1; i < points.size(); ++i) {
        double dt = 0.0;
        double dist = 0.0;
        if (!GetSDDeltaTimeSec(points[i - 1], points[i], dt)
            || !GetSDPointDistanceM(points[i - 1], points[i], dist)) {
            continue;
        }

        speeds[i] = dist / dt;
        sum += speeds[i];
        ++count;
    }

    if (count < 3) {
        return;
    }

    const double mean = sum / count;
    double squareSum = 0.0;
    for (double speed : speeds) {
        if (std::isfinite(speed)) {
            const double d = speed - mean;
            squareSum += d * d;
        }
    }

    const double sigma = std::sqrt(squareSum / std::max(count - 1, 1));
    if (!std::isfinite(sigma) || sigma < 1e-6) {
        return;
    }

    const double softThreshold = 2.0 * sigma;
    const double hardThreshold = 3.0 * sigma;

    for (size_t i = 1; i < points.size(); ++i) {
        if (!std::isfinite(speeds[i])) {
            continue;
        }

        const double residual = std::fabs(speeds[i] - mean);
        points[i].quality_score = BuildSDQualityByResidual(
            residual,
            softThreshold,
            hardThreshold);
        points[i].trusted = residual <= hardThreshold;
        points[i].trust_reason = points[i].trusted ? 0 : 10;
    }
}

static void RunSpeedPositionPredictionDetection(std::vector<CarPoint>& points)
{
    ResetSpeedDetectionResult(points);
    if (points.size() < 2) {
        return;
    }

    const double softThreshold = 5.0;
    const double hardThreshold = 12.0;

    for (size_t i = 1; i < points.size(); ++i) {
        const CarPoint& prev = points[i - 1];
        CarPoint& curr = points[i];

        double dt = 0.0;
        if (!GetSDDeltaTimeSec(prev, curr, dt)
            || !std::isfinite(prev.x) || !std::isfinite(prev.y)
            || !std::isfinite(curr.x) || !std::isfinite(curr.y)
            || !std::isfinite(prev.V_ENU[0]) || !std::isfinite(prev.V_ENU[1])) {
            continue;
        }

        const double predX = prev.x + prev.V_ENU[0] * dt;
        const double predY = prev.y + prev.V_ENU[1] * dt;
        const double dx = curr.x - predX;
        const double dy = curr.y - predY;
        const double residual = std::sqrt(dx * dx + dy * dy);

        curr.quality_score = BuildSDQualityByResidual(
            residual,
            softThreshold,
            hardThreshold);
        curr.trusted = residual <= hardThreshold;
        curr.trust_reason = curr.trusted ? 0 : 20;
    }
}

/**
 * 函数功能：
 *     SD 主流程入口。所有独立 SD 算法都通过 mode 枚举在这里分发，
 *     并直接写回质量评分和可信标记。
 *
 * 输入：
 *     points：车辆状态点数组，函数会原地修改。
 *     mode：SD 模式枚举。
 *
 * 输出：
 *     无返回值；修改 points 中每个点的 trusted、trust_reason、quality_score。
 *
 * 调用位置：
 *     1. EvaluateSpeedDetection 内部调用，用于 SD 效果验证。
 *     2. HMM.cpp 的 HiddenMarkovMapMatching 在 enableSD=true 时调用，用质量评分辅助 HMM。
 *
 * 数据流作用：
 *     CarPoint 原始轨迹 -> 带质量评分的 CarPoint。
 */
void RunSpeedDetection(std::vector<CarPoint>& points, SDMode mode)
{
    if (points.empty()) {
        return;
    }

    switch (mode) {
    case SDMode::Off:
        ResetSpeedDetectionResult(points);
        break;
    case SDMode::ThreeSigma:
        RunThreeSigmaSpeedDetection(points);
        break;
    case SDMode::SpeedPositionPrediction:
        RunSpeedPositionPredictionDetection(points);
        break;
    default:
        std::cerr << "未知 SD 模式，已按关闭处理。" << std::endl;
        ResetSpeedDetectionResult(points);
        break;
    }
}

/**
 * 函数功能：
 *     SD 验证入口。函数内部先调用 RunSpeedDetection 主流程，
 *     再将 SD 前后的 GNSS 点与参考真值按时间对齐，输出召回率、精确率等指标。
 *
 * 输入：
 *     coords：GNSS 原始观测数组。
 *     truthdata：参考真值数组。
 *     mode：SD 模式枚举。
 *     truthAbnormalThresholdM：真值异常阈值，单位 m。
 *     debugOutputFile：逐点调试结果输出路径；为空时不输出。
 *
 * 输出：
 *     无返回值；评估指标输出到控制台，可选逐点明细写入 debugOutputFile。
 *
 * 调用位置：
 *     main.cpp 的模式 1：SD 效果验证。
 *
 * 数据流作用：
 *     CoorData -> CarPoint -> SD 主流程 -> TxtData 对齐评估。
 */
void EvaluateSpeedDetection(
    const std::vector<CoorData>& coords,
    const std::vector<TxtData>& truthdata,
    SDMode mode,
    double truthAbnormalThresholdM,
    const std::string& debugOutputFile)
{
    if (coords.empty() || truthdata.empty()) {
        std::cout << "SD验证跳过：GNSS 点或真值为空。" << std::endl;
        return;
    }

    std::vector<CarPoint> beforePoints = ConvertCoorToCarPoints(coords);
    std::vector<CarPoint> afterPoints = beforePoints;

    RunSpeedDetection(afterPoints, mode);

    size_t truthIdx = 0;
    int matchedCount = 0;
    int truthAbnormalCount = 0;
    int sdAbnormalCount = 0;
    int trustedCount = 0;
    int truePositive = 0;
    int falsePositive = 0;
    int falseNegative = 0;
    int trueNegative = 0;
    int afterMatchedCount = 0;
    double beforeSumError = 0.0;
    double beforeSumSqError = 0.0;
    double beforeMaxError = 0.0;
    double afterSumError = 0.0;
    double afterSumSqError = 0.0;
    double afterMaxError = 0.0;

    std::ofstream debug;
    if (!debugOutputFile.empty()) {
        debug.open(debugOutputFile.c_str());
        if (debug.is_open()) {
            debug << "index\ttimestamp_ms\tbefore_error_m\tafter_error_m\t"
                << "truth_abnormal\tsd_abnormal\ttrusted\tquality_score\n";
            debug << std::fixed << std::setprecision(8);
        }
        else {
            std::cerr << "无法打开 SD 验证调试输出文件: " << debugOutputFile << std::endl;
        }
    }

    for (size_t i = 0; i < afterPoints.size(); ++i) {
        const CarPoint& before = beforePoints[i];
        const CarPoint& after = afterPoints[i];

        auto AbsTimeDiff = [](long long a, long long b) -> long long {
            return a >= b ? a - b : b - a;
        };

        while (truthIdx + 1 < truthdata.size()
            && AbsTimeDiff(truthdata[truthIdx + 1].timestamp, after.timestamp)
            <= AbsTimeDiff(truthdata[truthIdx].timestamp, after.timestamp)) {
            ++truthIdx;
        }

        const TxtData& t = truthdata[truthIdx];
        if (!std::isfinite(t.BLH[0]) || !std::isfinite(t.BLH[1])
            || !std::isfinite(after.BLH[0]) || !std::isfinite(after.BLH[1])) {
            continue;
        }

        const long long timeDiff = AbsTimeDiff(t.timestamp, after.timestamp);
        if (timeDiff > 2000) {
            continue;
        }

        const double beforeErr = haversine(before.BLH[0], before.BLH[1], t.BLH[0], t.BLH[1]);
        const double afterErr = haversine(after.BLH[0], after.BLH[1], t.BLH[0], t.BLH[1]);
        const bool truthAbnormal = afterErr > truthAbnormalThresholdM;
        const bool sdAbnormal = !after.trusted || after.quality_score < 0.35;

        beforeSumError += beforeErr;
        beforeSumSqError += beforeErr * beforeErr;
        beforeMaxError = std::max(beforeMaxError, beforeErr);
        ++matchedCount;

        if (truthAbnormal) {
            ++truthAbnormalCount;
        }
        if (sdAbnormal) {
            ++sdAbnormalCount;
        }
        if (after.trusted) {
            ++trustedCount;
        }

        if (!sdAbnormal) {
            afterSumError += afterErr;
            afterSumSqError += afterErr * afterErr;
            afterMaxError = std::max(afterMaxError, afterErr);
            ++afterMatchedCount;
        }

        if (sdAbnormal && truthAbnormal) {
            ++truePositive;
        }
        else if (sdAbnormal && !truthAbnormal) {
            ++falsePositive;
        }
        else if (!sdAbnormal && truthAbnormal) {
            ++falseNegative;
        }
        else {
            ++trueNegative;
        }

        if (debug.is_open()) {
            debug << i << "\t"
                << after.timestamp << "\t"
                << beforeErr << "\t"
                << afterErr << "\t"
                << (truthAbnormal ? 1 : 0) << "\t"
                << (sdAbnormal ? 1 : 0) << "\t"
                << (after.trusted ? 1 : 0) << "\t"
                << after.quality_score << "\n";
        }
    }

    if (debug.is_open()) {
        debug.close();
    }

    const double beforeAvgError = matchedCount > 0 ? beforeSumError / matchedCount : 0.0;
    const double beforeRmse = matchedCount > 0 ? std::sqrt(beforeSumSqError / matchedCount) : 0.0;
    const double afterAvgError = afterMatchedCount > 0 ? afterSumError / afterMatchedCount : 0.0;
    const double afterRmse = afterMatchedCount > 0 ? std::sqrt(afterSumSqError / afterMatchedCount) : 0.0;
    const double sdAbnormalRatio = matchedCount > 0
        ? static_cast<double>(sdAbnormalCount) / matchedCount * 100.0
        : 0.0;
    const double truthAbnormalRatio = matchedCount > 0
        ? static_cast<double>(truthAbnormalCount) / matchedCount * 100.0
        : 0.0;
    const double trustedRatio = matchedCount > 0
        ? static_cast<double>(trustedCount) / matchedCount * 100.0
        : 0.0;
    const double recall = (truePositive + falseNegative) > 0
        ? static_cast<double>(truePositive) / (truePositive + falseNegative) * 100.0
        : 0.0;
    const double precision = (truePositive + falsePositive) > 0
        ? static_cast<double>(truePositive) / (truePositive + falsePositive) * 100.0
        : 0.0;
    const double f1 = (precision + recall) > 0.0
        ? 2.0 * precision * recall / (precision + recall)
        : 0.0;

    std::cout << "\n========== SD验证 ==========\n";
    std::cout << "SD模式: " << GetSDModeName(mode) << "\n";
    std::cout << "输入点数: " << coords.size() << "\n";
    std::cout << "真值点数: " << truthdata.size() << "\n";
    std::cout << "时间对齐点数: " << matchedCount << "\n";
    std::cout << "真值异常阈值: " << truthAbnormalThresholdM << " m\n";
    std::cout << "SD前平均位置误差: " << beforeAvgError << " m\n";
    std::cout << "SD前RMSE: " << beforeRmse << " m\n";
    std::cout << "SD前最大位置误差: " << beforeMaxError << " m\n";
    std::cout << "SD后保留点数: " << afterMatchedCount << "\n";
    std::cout << "SD后保留点平均位置误差: " << afterAvgError << " m\n";
    std::cout << "SD后保留点RMSE: " << afterRmse << " m\n";
    std::cout << "SD后保留点最大位置误差: " << afterMaxError << " m\n";
    std::cout << "真值异常点数: " << truthAbnormalCount << "\n";
    std::cout << "真值异常点占比: " << truthAbnormalRatio << "%\n";
    std::cout << "SD异常点数: " << sdAbnormalCount << "\n";
    std::cout << "SD异常点占比: " << sdAbnormalRatio << "%\n";
    std::cout << "可信点占比: " << trustedRatio << "%\n";
    std::cout << "TP: " << truePositive << "\n";
    std::cout << "FP: " << falsePositive << "\n";
    std::cout << "FN: " << falseNegative << "\n";
    std::cout << "TN: " << trueNegative << "\n";
    std::cout << "召回率: " << recall << "%\n";
    std::cout << "精确率: " << precision << "%\n";
    std::cout << "F1: " << f1 << "%\n";

    if (!debugOutputFile.empty()) {
        std::cout << "SD验证调试文件: " << debugOutputFile << std::endl;
    }
}

