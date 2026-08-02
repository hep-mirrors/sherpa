#pragma once

#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <iostream>

/**
 * @class MonteCarloWeights
 * @brief A class for storing and tracking weights in Monte Carlo simulations
 * 
 * Provides functionality to:
 * - Store weights from Monte Carlo samples
 * - Calculate weight statistics (mean, variance, min, max)
 * - Track weight distribution and normalization
 * - Detect weight degeneracy issues
 * - Support incremental updates
 */
class MonteCarloWeights {
private:
    std::vector<double> weights_;
    double sum_weights_;
    double sum_squared_weights_;
    bool is_normalized_;
    
    // Weight monitoring
    std::vector<double> mean_history_;
    std::vector<double> jump_ratios_;
    std::vector<double> ess_history_;
    std::vector<double> variance_history_;
    std::vector<double> min_history_;
    std::vector<double> max_history_;
    double previous_mean_;
    bool has_previous_mean_;
    double max_jump_ratio_;
    size_t max_jump_index_;
    
    // Statistical change detection
    struct StatisticSnapshot {
        double mean;
        double variance;
        double std_dev;
        double ess;
        double degeneracy_ratio;
        double min_weight;
        double max_weight;
        size_t num_weights;
    };
    
    std::vector<StatisticSnapshot> stat_snapshots_;
    StatisticSnapshot last_snapshot_;
    bool has_last_snapshot_;
    
    /**
     * @brief Update running statistics when weights change
     */
    void updateStatistics() {
        if (weights_.empty()) {
            sum_weights_ = 0.0;
            sum_squared_weights_ = 0.0;
            return;
        }
        
        sum_weights_ = std::accumulate(weights_.begin(), weights_.end(), 0.0);
        sum_squared_weights_ = std::inner_product(
            weights_.begin(), weights_.end(), 
            weights_.begin(), 0.0
        );
    }
    
    /**
     * @brief Record mean and detect jump from previous mean
     */
    void recordMeanJump() {
        if (weights_.empty()) {
            return;
        }
        
        double current_mean = getMean();
        mean_history_.push_back(current_mean);
        
        if (has_previous_mean_ && previous_mean_ > 0.0) {
            double jump_ratio = current_mean / previous_mean_;
            jump_ratios_.push_back(jump_ratio);
            
            // Track maximum jump
            if (jump_ratio > max_jump_ratio_) {
                max_jump_ratio_ = jump_ratio;
                max_jump_index_ = jump_ratios_.size() - 1;
            }
        }
        
        previous_mean_ = current_mean;
        has_previous_mean_ = true;
    }
    
    /**
     * @brief Create a snapshot of all statistics
     */
    StatisticSnapshot createSnapshot() {
        StatisticSnapshot snap;
        snap.num_weights = weights_.size();
        snap.mean = getMean();
        snap.variance = getVariance();
        snap.std_dev = getStdDev();
        snap.ess = getEffectiveSampleSize();
        snap.degeneracy_ratio = getDegeneracyRatio();
        snap.min_weight = getMin();
        snap.max_weight = getMax();
        return snap;
    }

public:
    /**
     * @brief Default constructor
     */
    MonteCarloWeights() 
        : sum_weights_(0.0), sum_squared_weights_(0.0), is_normalized_(false),
          previous_mean_(0.0), has_previous_mean_(false),
          max_jump_ratio_(1.0), max_jump_index_(0), has_last_snapshot_(false) {}
    
    /**
     * @brief Constructor with initial capacity
     * @param capacity Initial capacity to reserve
     */
    explicit MonteCarloWeights(size_t capacity) 
        : sum_weights_(0.0), sum_squared_weights_(0.0), is_normalized_(false),
          previous_mean_(0.0), has_previous_mean_(false),
          max_jump_ratio_(1.0), max_jump_index_(0), has_last_snapshot_(false) {
        weights_.reserve(capacity);
    }
    
    /**
     * @brief Constructor from a vector of weights
     * @param weights Initial weights
     */
    explicit MonteCarloWeights(const std::vector<double>& weights)
        : weights_(weights), is_normalized_(false),
          previous_mean_(0.0), has_previous_mean_(false),
          max_jump_ratio_(1.0), max_jump_index_(0), has_last_snapshot_(false) {
        updateStatistics();
    }
    
    // ===== Basic Operations =====
    
    /**
     * @brief Add a single weight
     * @param weight The weight value to add
     */
    void addWeight(double weight) {
        // if (weight < 0.0) {
        //     throw std::invalid_argument("Weight cannot be negative");
        // }
        weights_.push_back(weight);
        sum_weights_ += weight;
        sum_squared_weights_ += weight * weight;
        is_normalized_ = false;
    }
    
    /**
     * @brief Add multiple weights
     * @param new_weights Vector of weights to add
     */
    void addWeights(const std::vector<double>& new_weights) {
        for (double w : new_weights) {
            // if (w < 0.0) {
            //     throw std::invalid_argument("Weight cannot be negative");
            // }
            weights_.push_back(w);
            sum_weights_ += w;
            sum_squared_weights_ += w * w;
        }
        is_normalized_ = false;
        recordMeanJump();
    }
    
    /**
     * @brief Clear all weights
     */
    void clear() {
        weights_.clear();
        sum_weights_ = 0.0;
        sum_squared_weights_ = 0.0;
        is_normalized_ = false;
        mean_history_.clear();
        jump_ratios_.clear();
        ess_history_.clear();
        variance_history_.clear();
        min_history_.clear();
        max_history_.clear();
        previous_mean_ = 0.0;
        has_previous_mean_ = false;
        max_jump_ratio_ = 1.0;
        max_jump_index_ = 0;
        stat_snapshots_.clear();
        has_last_snapshot_ = false;
    }
    
    // ===== Normalization =====
    
    /**
     * @brief Normalize weights to sum to 1.0
     */
    void normalize() {
        if (weights_.empty()) {
            throw std::runtime_error("Cannot normalize empty weights");
        }
        
        if (sum_weights_ <= 0.0) {
            throw std::runtime_error("Sum of weights must be positive");
        }
        
        for (double& w : weights_) {
            w /= sum_weights_;
        }
        
        updateStatistics();
        is_normalized_ = true;
    }
    
    /**
     * @brief Check if weights are normalized
     * @return true if normalized, false otherwise
     */
    bool isNormalized() const {
        return is_normalized_;
    }
    
    // ===== Statistics =====
    
    /**
     * @brief Get the number of weights
     * @return Number of weights stored
     */
    size_t size() const {
        return weights_.size();
    }
    
    /**
     * @brief Get the sum of all weights
     * @return Sum of weights
     */
    double getSumWeights() const {
        return sum_weights_;
    }
    
    /**
     * @brief Get the mean weight
     * @return Mean value of weights
     */
    double getMean() const {
        if (weights_.empty()) {
            throw std::runtime_error("Cannot compute mean of empty weights");
        }
        return sum_weights_ / weights_.size();
    }
    
    /**
     * @brief Get the variance of weights
     * @return Variance (population variance)
     */
    double getVariance() const {
        if (weights_.empty()) {
            throw std::runtime_error("Cannot compute variance of empty weights");
        }
        
        double mean = getMean();
        double variance = (sum_squared_weights_ / weights_.size()) - (mean * mean);
        
        // Handle floating point errors
        return std::max(0.0, variance);
    }
    
    /**
     * @brief Get the standard deviation of weights
     * @return Standard deviation
     */
    double getStdDev() const {
        return std::sqrt(getVariance());
    }
    
    /**
     * @brief Get the minimum weight
     * @return Minimum weight value
     */
    double getMin() const {
        if (weights_.empty()) {
            throw std::runtime_error("Cannot get min of empty weights");
        }
        return *std::min_element(weights_.begin(), weights_.end());
    }
    
    /**
     * @brief Get the maximum weight
     * @return Maximum weight value
     */
    double getMax() const {
        if (weights_.empty()) {
            throw std::runtime_error("Cannot get max of empty weights");
        }
        return *std::max_element(weights_.begin(), weights_.end());
    }
    
    /**
     * @brief Get the effective sample size
     * 
     * ESS = (sum of weights)^2 / (sum of squared weights)
     * Useful for detecting weight degeneracy
     * 
     * @return Effective sample size
     */
    double getEffectiveSampleSize() const {
        if (sum_squared_weights_ <= 0.0) {
            throw std::runtime_error("Sum of squared weights must be positive");
        }
        return (sum_weights_ * sum_weights_) / sum_squared_weights_;
    }
    
    /**
     * @brief Get the weight degeneracy ratio
     * 
     * Ratio of ESS to actual number of samples.
     * Values close to 1.0 indicate good diversity.
     * Values close to 0.0 indicate degeneracy.
     * 
     * @return Degeneracy ratio (ESS / N)
     */
    double getDegeneracyRatio() const {
        if (weights_.empty()) {
            throw std::runtime_error("Cannot compute degeneracy ratio of empty weights");
        }
        return getEffectiveSampleSize() / weights_.size();
    }
    
    /**
     * @brief Get the coefficient of variation
     * 
     * CV = std_dev / mean
     * Higher values indicate more variation
     * 
     * @return Coefficient of variation
     */
    double getCoefficientOfVariation() const {
        double mean = getMean();
        if (mean == 0.0) {
            throw std::runtime_error("Mean is zero, cannot compute CV");
        }
        return getStdDev() / mean;
    }
    
/**
     * @brief Remove last weight
     * @throws std::runtime_error if weights list is empty
     */
    void removeLastWeight() {
        if (weights_.empty()) {
            throw std::runtime_error("Cannot remove from empty weights");
        }
        
        double removed = weights_.back();
        weights_.pop_back();
        sum_weights_ -= removed;
        sum_squared_weights_ -= removed * removed;
        is_normalized_ = false;
    }
    

    // ===== Access =====
    
    /**
     * @brief Access a weight by index
     * @param index Index of the weight
     * @return The weight at the given index
     */
    double operator[](size_t index) const {
        if (index >= weights_.size()) {
            throw std::out_of_range("Weight index out of range");
        }
        return weights_[index];
    }
    
    /**
     * @brief Get a const reference to the weights vector
     * @return Const reference to weights
     */
    const std::vector<double>& getWeights() const {
        return weights_;
    }
    
    /**
     * @brief Get a mutable reference to the weights vector
     * 
     * Warning: Direct modification may invalidate statistics
     * Call updateStatistics() after modification if needed
     * 
     * @return Mutable reference to weights
     */
    std::vector<double>& getWeightsMutable() {
        is_normalized_ = false;
        return weights_;
    }
    
    // ===== Diagnostics =====
    
    /**
     * @brief Check if weights show significant degeneracy
     * @param threshold Degeneracy ratio threshold (default 0.1)
     * @return true if degeneracy ratio is below threshold
     */
    bool isDegenerative(double threshold = 0.1) const {
        return (getDegeneracyRatio() < threshold && weights_.size() > 100);
    }
    
    // ===== Weight Jump Monitoring =====
    
    /**
     * @brief Manually record the current mean weight for jump detection
     * 
     * Call this method between batches/iterations to track mean changes
     */
    void recordWeightSnapshot() {
        recordMeanJump();
    }
    
    /**
     * @brief Record a complete snapshot of all statistics for change detection
     * 
     * This captures: mean, variance, ESS, degeneracy ratio, min, max
     * Use this to detect changes in ANY statistical measure, not just mean
     */
    void recordStatisticSnapshot() {
        if (weights_.empty()) {
            throw std::runtime_error("Cannot snapshot empty weights");
        }
        
        StatisticSnapshot snap = createSnapshot();
        stat_snapshots_.push_back(snap);
        
        // Track individual statistics in history
        mean_history_.push_back(snap.mean);
        ess_history_.push_back(snap.ess);
        variance_history_.push_back(snap.variance);
        min_history_.push_back(snap.min_weight);
        max_history_.push_back(snap.max_weight);
        
        has_last_snapshot_ = true;
        last_snapshot_ = snap;
    }
    
    /**
     * @brief Detect significant change in ANY statistic
     * 
     * Compares current snapshot with previous one. Triggers on:
     * - Mean change > threshold
     * - ESS drop > threshold  
     * - Variance increase > threshold
     * - Min/Max range change > threshold
     * 
     * @param metric Which metric to check: "mean", "ess", "variance", "range"
     * @param percent_threshold Percent change threshold (default 10%)
     * @return Pair of (changed, percent_change)
     */
    std::pair<bool, double> detectStatisticChange(
        const std::string& metric = "mean",
        double percent_threshold = 10.0) const {
        
        if (stat_snapshots_.size() < 2) {
            return {false, 0.0};
        }
        
        const StatisticSnapshot& prev = stat_snapshots_[stat_snapshots_.size() - 2];
        const StatisticSnapshot& curr = last_snapshot_;
        
        double percent_change = 0.0;
        
        if (metric == "mean") {
            if (prev.mean != 0.0) {
                percent_change = ((curr.mean - prev.mean) / std::abs(prev.mean)) * 100.0;
            }
        } else if (metric == "ess") {
            if (prev.ess != 0.0) {
                percent_change = ((curr.ess - prev.ess) / prev.ess) * 100.0;
            }
        } else if (metric == "variance") {
            if (prev.variance != 0.0) {
                percent_change = ((curr.variance - prev.variance) / prev.variance) * 100.0;
            }
        } else if (metric == "range") {
            double prev_range = prev.max_weight - prev.min_weight;
            double curr_range = curr.max_weight - curr.min_weight;
            if (prev_range != 0.0) {
                percent_change = ((curr_range - prev_range) / prev_range) * 100.0;
            }
        } else if (metric == "degeneracy") {
            if (prev.degeneracy_ratio != 0.0) {
                percent_change = ((curr.degeneracy_ratio - prev.degeneracy_ratio) / 
                                 prev.degeneracy_ratio) * 100.0;
            }
        }
        
        bool significant = std::abs(percent_change) > percent_threshold;
        return {significant, percent_change};
    }
    
    /**
     * @brief Get comprehensive change report for all statistics
     */
    struct StatisticChangeReport {
        bool mean_changed;
        double mean_percent_change;
        
        bool ess_changed;
        double ess_percent_change;
        
        bool variance_changed;
        double variance_percent_change;
        
        bool range_changed;
        double range_percent_change;
        
        bool degeneracy_changed;
        double degeneracy_percent_change;
        
        // Thresholds used
        double threshold;
        
        // Overall change detected
        bool any_change;
        int num_changes;
    };
    
    /**
     * @brief Generate comprehensive change detection report
     * @param percent_threshold Percent change threshold (default 10%)
     * @return Report structure with all changes
     */
    StatisticChangeReport getStatisticChangeReport(double percent_threshold = 10.0) const {
        StatisticChangeReport report;
        report.threshold = percent_threshold;
        report.num_changes = 0;
        
        if (stat_snapshots_.size() < 2) {
            report.any_change = false;
            return report;
        }
        
        auto [mean_ch, mean_pct] = detectStatisticChange("mean", percent_threshold);
        report.mean_changed = mean_ch;
        report.mean_percent_change = mean_pct;
        if (mean_ch) report.num_changes++;
        
        auto [ess_ch, ess_pct] = detectStatisticChange("ess", percent_threshold);
        report.ess_changed = ess_ch;
        report.ess_percent_change = ess_pct;
        if (ess_ch) report.num_changes++;
        
        auto [var_ch, var_pct] = detectStatisticChange("variance", percent_threshold);
        report.variance_changed = var_ch;
        report.variance_percent_change = var_pct;
        if (var_ch) report.num_changes++;
        
        auto [rng_ch, rng_pct] = detectStatisticChange("range", percent_threshold);
        report.range_changed = rng_ch;
        report.range_percent_change = rng_pct;
        if (rng_ch) report.num_changes++;
        
        auto [deg_ch, deg_pct] = detectStatisticChange("degeneracy", percent_threshold);
        report.degeneracy_changed = deg_ch;
        report.degeneracy_percent_change = deg_pct;
        if (deg_ch) report.num_changes++;
        
        report.any_change = (report.num_changes > 0);
        return report;
    }
    
    /**
     * @brief Get the last recorded statistics snapshot
     */
    const StatisticSnapshot& getLastSnapshot() const {
        if (!has_last_snapshot_) {
            throw std::runtime_error("No snapshot recorded yet");
        }
        return last_snapshot_;
    }
    
    /**
     * @brief Get all statistics snapshots
     */
    const std::vector<StatisticSnapshot>& getAllSnapshots() const {
        return stat_snapshots_;
    }
    
    /**
     * @brief Get ESS history
     */
    const std::vector<double>& getESSHistory() const {
        return ess_history_;
    }
    
    /**
     * @brief Get variance history
     */
    const std::vector<double>& getVarianceHistory() const {
        return variance_history_;
    }
    
    /**
     * @brief Get min weight history
     */
    const std::vector<double>& getMinHistory() const {
        return min_history_;
    }
    
    /**
     * @brief Get max weight history
     */
    const std::vector<double>& getMaxHistory() const {
        return max_history_;
    }
    
    /**
     * @brief Print detailed change detection report
     * @param percent_threshold Threshold for significant change (default 10%)
     */
    void printDetailedChangeReport(double percent_threshold = 10.0) const {
        if (!has_last_snapshot_) {
            std::cout << "No statistics snapshots recorded yet\n";
            return;
        }
        
        if (stat_snapshots_.size() < 2) {
            std::cout << "Need at least 2 snapshots for change detection\n";
            return;
        }
        
        auto report = getStatisticChangeReport(percent_threshold);
        
        std::cout << "=== Statistic Change Detection Report ===\n";
        std::cout << "Threshold: " << percent_threshold << "%\n";
        std::cout << "Snapshots: " << stat_snapshots_.size() << "\n";
        std::cout << "Changes detected: " << report.num_changes << "\n\n";
        
        const StatisticSnapshot& prev = stat_snapshots_[stat_snapshots_.size() - 2];
        const StatisticSnapshot& curr = last_snapshot_;
        
        std::cout << "Previous Snapshot | Current Snapshot | Change\n";
        std::cout << "--------\n";
        
        std::cout << "MEAN: " << prev.mean << " → " << curr.mean;
        std::cout << " | " << (report.mean_percent_change > 0 ? "↑" : "↓") 
                 << " " << std::abs(report.mean_percent_change) << "%";
        if (report.mean_changed) std::cout << " ⚠️ CHANGED";
        std::cout << "\n";
        
        std::cout << "ESS:  " << prev.ess << " → " << curr.ess;
        std::cout << " | " << (report.ess_percent_change > 0 ? "↑" : "↓") 
                 << " " << std::abs(report.ess_percent_change) << "%";
        if (report.ess_changed) std::cout << " ⚠️ CHANGED";
        std::cout << "\n";
        
        std::cout << "VAR:  " << prev.variance << " → " << curr.variance;
        std::cout << " | " << (report.variance_percent_change > 0 ? "↑" : "↓") 
                 << " " << std::abs(report.variance_percent_change) << "%";
        if (report.variance_changed) std::cout << " ⚠️ CHANGED";
        std::cout << "\n";
        
        double prev_range = prev.max_weight - prev.min_weight;
        double curr_range = curr.max_weight - curr.min_weight;
        std::cout << "RANGE: " << prev_range << " → " << curr_range;
        std::cout << " | " << (report.range_percent_change > 0 ? "↑" : "↓") 
                 << " " << std::abs(report.range_percent_change) << "%";
        if (report.range_changed) std::cout << " ⚠️ CHANGED";
        std::cout << "\n";
        
        std::cout << "DEG RATIO: " << prev.degeneracy_ratio << " → " << curr.degeneracy_ratio;
        std::cout << " | " << (report.degeneracy_percent_change > 0 ? "↑" : "↓") 
                 << " " << std::abs(report.degeneracy_percent_change) << "%";
        if (report.degeneracy_changed) std::cout << " ⚠️ CHANGED";
        std::cout << "\n";
        
        std::cout << "\n";
        if (report.any_change) {
            std::cout << "🚨 OVERALL: " << report.num_changes << " statistics changed\n";
        } else {
            std::cout << "✓ OVERALL: No significant changes\n";
        }
        std::cout << "=========================================\n";
    }
    
    /**
     * @brief Print weight statistics to console
     */
    void printStatistics() const {
        if (weights_.empty()) {
            std::cout << "No weights stored\n";
            return;
        }
        
        std::cout << "=== Monte Carlo Weight Statistics ===\n";
        std::cout << "Number of weights: " << weights_.size() << "\n";
        std::cout << "Sum of weights: " << sum_weights_ << "\n";
        std::cout << "Mean weight: " << getMean() << "\n";
        std::cout << "Std deviation: " << getStdDev() << "\n";
        std::cout << "Min weight: " << getMin() << "\n";
        std::cout << "Max weight: " << getMax() << "\n";
        std::cout << "Effective sample size: " << getEffectiveSampleSize() << "\n";
        std::cout << "Degeneracy ratio: " << getDegeneracyRatio() << "\n";
        std::cout << "Normalized: " << (is_normalized_ ? "Yes" : "No") << "\n";
        
        if (isDegenerative()) {
            std::cout << "WARNING: Weights show significant degeneracy!\n";
        }
        
        // Jump monitoring info
        // if (!jump_ratios_.empty()) {
        //     std::cout << "\n--- Weight Jump History ---\n";
        //     std::cout << "Number of snapshots: " << mean_history_.size() << "\n";
        //     std::cout << "Last jump ratio: " << getLastJumpRatio() 
        //               << " (" << (getLastJumpRatio() > 1.0 ? "↑" : "↓") << ")\n";
        //     std::cout << "Max jump ratio: " << getMaxJumpRatio() 
        //               << " (at snapshot " << getMaxJumpIndex() << ")\n";
            
        //     // if (hasSignificantJump()) {
        //     //     std::cout << "⚠️  SIGNIFICANT JUMP DETECTED! Last mean changed by "
        //     //               << (std::abs(getLastJumpRatio() - 1.0) * 100.0) << "%\n";
        //     // }
            
        //     auto [mean_j, max_j] = getJumpStatistics();
        //     std::cout << "Average jump ratio: " << mean_j << "\n";
        // }
        
        std::cout << "====================================\n";
    }
    
    // ===== Utility =====
    
    /**
     * @brief Get the index of the maximum weight
     * @return Index of maximum weight
     */
    size_t getMaxIndex() const {
        if (weights_.empty()) {
            throw std::runtime_error("Cannot get max index of empty weights");
        }
        return std::distance(weights_.begin(), 
                           std::max_element(weights_.begin(), weights_.end()));
    }
    
    /**
     * @brief Manually update statistics after direct weight modifications
     */
    void updateStats() {
        updateStatistics();
    }
    
    /**
     * @brief Get a copy of weights for analysis
     * @return Copy of the weights vector
     */
    std::vector<double> getWeightsCopy() const {
        return weights_;
    }
};