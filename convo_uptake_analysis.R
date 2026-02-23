###############################################################################
# Minimum Observations for Reliable Inferences About Instructional Uptake
#
# Generalizability analysis of the conversational uptake dataset from:
#   Demszky et al. (2021) "Measuring Conversational Uptake" (ACL 2021)
#   Demszky & Hill (2023) "The NCTE Transcripts" (BEA 2023)
#
# Data source: https://github.com/ddemszky/conversational-uptake
#
# Author: Generated for AIED 2026 analysis
# Date: 2026-02-23
###############################################################################

# ── 0. Setup ─────────────────────────────────────────────────────────────────

library(tidyverse)
library(lme4)       # mixed-effects models for ICC
library(boot)       # bootstrapping

# If packages aren't installed:
# install.packages(c("tidyverse", "lme4", "boot"))

# Set seed for reproducibility
set.seed(42)

# ── 1. Load Data ─────────────────────────────────────────────────────────────

# Clone the repo first:
#   git clone https://github.com/ddemszky/conversational-uptake.git
# Then set the path below:

data_path <- "conversational-uptake/data/uptake_data.csv"
df <- read_csv(data_path, show_col_types = FALSE)

cat("=== DATASET OVERVIEW ===\n")
cat("Total exchanges:", nrow(df), "\n")
cat("Unique observations (obs_id):", n_distinct(df$obs_id), "\n")
cat("Columns:", paste(names(df), collapse = ", "), "\n\n")

# ── 2. Descriptive Statistics ────────────────────────────────────────────────

cat("=== DESCRIPTIVE STATISTICS ===\n\n")

# Exchanges per observation
exchanges_per_obs <- df %>%
    count(obs_id, name = "n_exchanges")

cat("Exchanges per observation:\n")
print(summary(exchanges_per_obs$n_exchanges))
cat("\n")

# Distribution table
exchange_bins <- exchanges_per_obs %>%
    mutate(bin = case_when(
        n_exchanges == 1 ~ "1",
        n_exchanges == 2 ~ "2",
        n_exchanges == 3 ~ "3",
        n_exchanges == 4 ~ "4",
        n_exchanges <= 9 ~ "5-9",
        n_exchanges <= 19 ~ "10-19",
        n_exchanges <= 49 ~ "20-49",
        TRUE ~ "50+"
    )) %>%
    count(bin) %>%
    mutate(pct = round(100 * n / sum(n), 1))

cat("Distribution of exchanges per observation:\n")
print(exchange_bins, n = Inf)
cat("\n")

# Inter-rater agreement
cat("Inter-rater agreement (uptake_num_agree):\n")
df %>%
    count(uptake_num_agree) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    print()
cat("\n")

# Uptake z-score summary
cat("Uptake z-score summary:\n")
df %>%
    filter(!is.na(uptake_zscore)) %>%
    summarise(
        n = n(),
        mean = round(mean(uptake_zscore), 4),
        sd = round(sd(uptake_zscore), 4),
        min = round(min(uptake_zscore), 4),
        q25 = round(quantile(uptake_zscore, 0.25), 4),
        median = round(median(uptake_zscore), 4),
        q75 = round(quantile(uptake_zscore, 0.75), 4),
        max = round(max(uptake_zscore), 4)
    ) %>%
    print()
cat("\n")

# ── 3. ICC(1) via ANOVA ─────────────────────────────────────────────────────
#
# One-way random-effects ICC(1) treating obs_id as the grouping variable.
# This decomposes variance into between-observation and within-observation
# components, analogous to a G-study in Generalizability Theory.

cat("=== ICC(1) ANALYSIS — ANOVA METHOD ===\n\n")

# Filter to observations with >=2 exchanges and valid z-scores
df_valid <- df %>% filter(!is.na(uptake_zscore))

multi_obs_ids <- df_valid %>%
    count(obs_id) %>%
    filter(n >= 2) %>%
    pull(obs_id)

df_multi <- df_valid %>% filter(obs_id %in% multi_obs_ids)

# One-way ANOVA
aov_result <- aov(uptake_zscore ~ factor(obs_id), data = df_multi)
aov_summary <- summary(aov_result)

SSB <- aov_summary[[1]]["factor(obs_id)", "Sum Sq"]
SSW <- aov_summary[[1]]["Residuals", "Sum Sq"]
dfB <- aov_summary[[1]]["factor(obs_id)", "Df"]
dfW <- aov_summary[[1]]["Residuals", "Df"]
MSB <- SSB / dfB
MSW <- SSW / dfW

# Compute n0 (harmonic-adjusted mean group size for unbalanced design)
n_per_group <- df_multi %>% count(obs_id) %>% pull(n)
N <- sum(n_per_group)
k <- length(n_per_group)
n0 <- (N - sum(n_per_group^2) / N) / (k - 1)

# ICC(1)
ICC1 <- (MSB - MSW) / (MSB + (n0 - 1) * MSW)

cat("Number of observations with >=2 exchanges:", k, "\n")
cat("Total exchanges in analysis:", N, "\n")
cat("Harmonic-adjusted mean group size (n0):", round(n0, 2), "\n")
cat("MSB (between observations):", round(MSB, 4), "\n")
cat("MSW (within observations):", round(MSW, 4), "\n")
cat("ICC(1):", round(ICC1, 4), "\n\n")

# Variance decomposition
sigma2_b <- (MSB - MSW) / n0
sigma2_w <- MSW

cat("=== VARIANCE DECOMPOSITION ===\n")
cat("Between-observation variance (sigma^2_b):", round(sigma2_b, 4), "\n")
cat("Within-observation variance (sigma^2_w):", round(sigma2_w, 4), "\n")
cat("% between:", round(100 * sigma2_b / (sigma2_b + sigma2_w), 1), "%\n")
cat("% within:", round(100 * sigma2_w / (sigma2_b + sigma2_w), 1), "%\n\n")

# ── 4. ICC(1) via Mixed-Effects Model (lme4) ────────────────────────────────
#
# Cross-check the ANOVA-based ICC with a random-intercept model.

cat("=== ICC(1) — MIXED-EFFECTS MODEL (lme4) ===\n\n")

mod <- lmer(uptake_zscore ~ 1 + (1 | obs_id), data = df_multi, REML = TRUE)

vc <- as.data.frame(VarCorr(mod))
sigma2_b_lmer <- vc$vcov[vc$grp == "obs_id"]
sigma2_w_lmer <- vc$vcov[vc$grp == "Residual"]
ICC1_lmer <- sigma2_b_lmer / (sigma2_b_lmer + sigma2_w_lmer)

cat("Between-observation variance:", round(sigma2_b_lmer, 4), "\n")
cat("Within-observation variance:", round(sigma2_w_lmer, 4), "\n")
cat("ICC(1) from lmer:", round(ICC1_lmer, 4), "\n")
cat("(Should closely match ANOVA ICC =", round(ICC1, 4), ")\n\n")

# ── 5. Spearman-Brown Reliability Projections (D-Study) ─────────────────────
#
# Given ICC(1) as single-exchange reliability, project the reliability
# of the mean of k exchanges using the Spearman-Brown prophecy formula:
#   Rel(k) = k * ICC / (1 + (k - 1) * ICC)

cat("=== SPEARMAN-BROWN RELIABILITY PROJECTIONS ===\n\n")

spearman_brown <- function(icc, k) {
    (k * icc) / (1 + (k - 1) * icc)
}

k_values <- c(1, 2, 3, 5, 7, 10, 15, 20, 25, 30, 40, 50, 75, 100)
sb_table <- tibble(
    k_exchanges = k_values,
    reliability = round(sapply(k_values, function(k) spearman_brown(ICC1, k)), 3)
)

cat("Projected reliability by number of exchanges:\n")
print(sb_table, n = Inf)
cat("\n")

# Minimum exchanges for target reliabilities
min_exchanges <- function(icc, target_rel) {
    ceiling(target_rel * (1 - icc) / (icc * (1 - target_rel)))
}

targets <- c(0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90)
target_table <- tibble(
    target_reliability = targets,
    min_exchanges_needed = sapply(targets, function(t) min_exchanges(ICC1, t))
)

cat("Minimum exchanges for target reliability thresholds:\n")
print(target_table, n = Inf)
cat("\n")

# ── 6. Sensitivity: Raw Uptake Scale (0–2) ──────────────────────────────────

cat("=== SENSITIVITY CHECK: RAW UPTAKE (0-2 scale) ===\n\n")

df_multi_raw <- df %>%
    filter(!is.na(uptake), obs_id %in% multi_obs_ids)

mod_raw <- lmer(uptake ~ 1 + (1 | obs_id), data = df_multi_raw, REML = TRUE)
vc_raw <- as.data.frame(VarCorr(mod_raw))
sigma2_b_raw <- vc_raw$vcov[vc_raw$grp == "obs_id"]
sigma2_w_raw <- vc_raw$vcov[vc_raw$grp == "Residual"]
ICC1_raw <- sigma2_b_raw / (sigma2_b_raw + sigma2_w_raw)

cat("ICC(1) on raw uptake:", round(ICC1_raw, 4), "\n\n")

raw_target_table <- tibble(
    target_reliability = targets,
    min_exchanges_zscore = sapply(targets, function(t) min_exchanges(ICC1, t)),
    min_exchanges_raw = sapply(targets, function(t) min_exchanges(ICC1_raw, t))
)

cat("Comparison — minimum exchanges by metric:\n")
print(raw_target_table, n = Inf)
cat("\n")

# ── 7. Split-Half Reliability ────────────────────────────────────────────────

cat("=== SPLIT-HALF RELIABILITY ===\n\n")

obs_10plus <- df_valid %>%
    count(obs_id) %>%
    filter(n >= 10) %>%
    pull(obs_id)

cat("Observations with >=10 exchanges:", length(obs_10plus), "\n")

split_half_data <- df_valid %>%
    filter(obs_id %in% obs_10plus) %>%
    group_by(obs_id) %>%
    mutate(
        row_id = row_number(),
        half = if_else(row_id %% 2 == 1, "odd", "even")
    ) %>%
    group_by(obs_id, half) %>%
    summarise(mean_uptake = mean(uptake_zscore), .groups = "drop") %>%
    pivot_wider(names_from = half, values_from = mean_uptake)

r_split <- cor(split_half_data$odd, split_half_data$even, use = "complete.obs")
sb_corrected <- 2 * r_split / (1 + r_split)

cat("Split-half Pearson r:", round(r_split, 3), "\n")
cat("Spearman-Brown corrected:", round(sb_corrected, 3), "\n")
cat("p-value:", round(cor.test(split_half_data$odd, split_half_data$even)$p.value, 4), "\n\n")

# ── 8. Bootstrap Stability Analysis ─────────────────────────────────────────
#
# For observations with >=30 exchanges, repeatedly subsample k exchanges
# and compute correlation between subsample mean and full-observation mean.

cat("=== BOOTSTRAP STABILITY ANALYSIS ===\n\n")

obs_30plus <- df_valid %>%
    count(obs_id) %>%
    filter(n >= 30) %>%
    pull(obs_id)

cat("Observations with >=30 exchanges:", length(obs_30plus), "\n\n")

if (length(obs_30plus) >= 5) {

    big_obs_data <- df_valid %>%
        filter(obs_id %in% obs_30plus) %>%
        group_by(obs_id) %>%
        summarise(
            values = list(uptake_zscore),
            full_mean = mean(uptake_zscore),
            n = n(),
            .groups = "drop"
        )

    n_boot <- 200
    k_test <- c(5, 10, 15, 20, 25, 30)

    boot_results <- map_dfr(k_test, function(k) {
        correlations <- replicate(n_boot, {
            sampled_means <- map_dbl(big_obs_data$values, function(vals) {
                if (length(vals) >= k) {
                    mean(sample(vals, k, replace = FALSE))
                } else {
                    NA_real_
                }
            })
            cor(sampled_means, big_obs_data$full_mean, use = "complete.obs")
        })

        tibble(
            k = k,
            mean_r = round(mean(correlations, na.rm = TRUE), 3),
            sd_r = round(sd(correlations, na.rm = TRUE), 3),
            ci_lower = round(quantile(correlations, 0.025, na.rm = TRUE), 3),
            ci_upper = round(quantile(correlations, 0.975, na.rm = TRUE), 3)
        )
    })

    cat("Bootstrap correlation with full-observation mean:\n")
    print(boot_results, n = Inf)
    cat("\n")

} else {
    cat("Insufficient observations with >=30 exchanges for bootstrap analysis.\n\n")
}

# ── 9. Visualization ────────────────────────────────────────────────────────

cat("=== GENERATING PLOTS ===\n\n")

# Plot 1: Spearman-Brown reliability curve
k_seq <- 1:150
rel_seq <- sapply(k_seq, function(k) spearman_brown(ICC1, k))

p1 <- ggplot(tibble(k = k_seq, reliability = rel_seq), aes(k, reliability)) +
    geom_line(color = "#2E75B6", linewidth = 1.2) +
    geom_hline(yintercept = c(0.60, 0.70, 0.80),
               linetype = "dashed", color = "grey50", linewidth = 0.5) +
    annotate("text", x = 140, y = 0.61, label = "0.60 (screening)",
             size = 3, color = "grey40", hjust = 1) +
    annotate("text", x = 140, y = 0.71, label = "0.70 (research)",
             size = 3, color = "grey40", hjust = 1) +
    annotate("text", x = 140, y = 0.81, label = "0.80 (high-stakes)",
             size = 3, color = "grey40", hjust = 1) +
    geom_point(data = target_table %>%
                   filter(target_reliability %in% c(0.60, 0.70, 0.80)) %>%
                   mutate(rel = target_reliability),
               aes(x = min_exchanges_needed, y = rel),
               color = "#C00000", size = 3) +
    labs(
        title = "Reliability of Mean Uptake Score by Number of Exchanges",
        subtitle = paste0("ICC(1) = ", round(ICC1, 3),
                          " | Based on Demszky et al. conversational uptake dataset"),
        x = "Number of Teacher-Student Exchanges (k)",
        y = "Projected Reliability"
    ) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
    )

ggsave("reliability_curve.png", p1, width = 8, height = 5, dpi = 150)
cat("Saved: reliability_curve.png\n")

# Plot 2: Variance decomposition pie chart
p2 <- ggplot(
    tibble(
        component = c("Between observations\n(teacher/lesson signal)",
                      "Within observations\n(exchange-to-exchange noise)"),
        variance = c(sigma2_b, sigma2_w),
        pct = round(100 * variance / sum(variance), 1)
    ),
    aes(x = "", y = variance, fill = component)
) +
    geom_col(width = 1) +
    coord_polar(theta = "y") +
    geom_text(aes(label = paste0(pct, "%")),
              position = position_stack(vjust = 0.5), size = 5, fontface = "bold") +
    scale_fill_manual(values = c("#2E75B6", "#D5E8F0")) +
    labs(
        title = "Variance Decomposition of Uptake Scores",
        subtitle = paste0("ICC(1) = ", round(ICC1, 3)),
        fill = NULL
    ) +
    theme_void(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom"
    )

ggsave("variance_decomposition.png", p2, width = 7, height = 5, dpi = 150)
cat("Saved: variance_decomposition.png\n")

# Plot 3: Bootstrap stability
if (exists("boot_results")) {
    p3 <- ggplot(boot_results, aes(k, mean_r)) +
        geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
                    fill = "#2E75B6", alpha = 0.2) +
        geom_line(color = "#2E75B6", linewidth = 1.2) +
        geom_point(color = "#2E75B6", size = 3) +
        geom_hline(yintercept = 0.80, linetype = "dashed", color = "grey50") +
        labs(
            title = "Bootstrap Stability: Subsample vs. Full-Observation Agreement",
            subtitle = paste0("Based on ", length(obs_30plus),
                              " observations with >=30 exchanges, 200 bootstrap samples"),
            x = "Subsample Size (k exchanges)",
            y = "Mean Correlation with Full-Observation Mean"
        ) +
        scale_y_continuous(limits = c(0.4, 1), breaks = seq(0.4, 1, 0.1)) +
        theme_minimal(base_size = 12) +
        theme(
            plot.title = element_text(face = "bold"),
            panel.grid.minor = element_blank()
        )

    ggsave("bootstrap_stability.png", p3, width = 8, height = 5, dpi = 150)
    cat("Saved: bootstrap_stability.png\n")
}

# Plot 4: Distribution of exchanges per observation
p4 <- ggplot(exchanges_per_obs, aes(n_exchanges)) +
    geom_histogram(binwidth = 1, fill = "#2E75B6", color = "white",
                   boundary = 0.5) +
    geom_vline(xintercept = c(21, 33), linetype = "dashed",
               color = c("#C00000", "#C00000")) +
    annotate("text", x = 23, y = Inf, label = "k=21\n(rel≥.60)",
             vjust = 1.5, size = 3, color = "#C00000") +
    annotate("text", x = 35, y = Inf, label = "k=33\n(rel≥.70)",
             vjust = 1.5, size = 3, color = "#C00000") +
    scale_x_continuous(breaks = c(1, 5, 10, 20, 33, 50, 75, 100)) +
    labs(
        title = "Distribution of Annotated Exchanges per Observation",
        subtitle = "Most observations have very few annotated exchanges",
        x = "Number of Exchanges",
        y = "Number of Observations"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
    )

ggsave("exchanges_distribution.png", p4, width = 8, height = 5, dpi = 150)
cat("Saved: exchanges_distribution.png\n")

# ── 10. Summary ──────────────────────────────────────────────────────────────

cat("\n")
cat("=" |> strrep(70), "\n")
cat("SUMMARY OF KEY FINDINGS\n")
cat("=" |> strrep(70), "\n\n")

cat("ICC(1) = ", round(ICC1, 3), "\n")
cat("  → Only ", round(100 * sigma2_b / (sigma2_b + sigma2_w), 1),
    "% of variance is between observations\n")
cat("  → ", round(100 * sigma2_w / (sigma2_b + sigma2_w), 1),
    "% is within-observation (exchange-to-exchange) noise\n\n")

cat("MINIMUM EXCHANGES FOR RELIABLE INFERENCES:\n")
cat("  Formative screening (rel ≥ 0.60): ",
    min_exchanges(ICC1, 0.60), " exchanges\n")
cat("  Research standard   (rel ≥ 0.70): ",
    min_exchanges(ICC1, 0.70), " exchanges\n")
cat("  High-stakes         (rel ≥ 0.80): ",
    min_exchanges(ICC1, 0.80), " exchanges\n")
cat("  Summative eval      (rel ≥ 0.90): ",
    min_exchanges(ICC1, 0.90), " exchanges\n\n")

cat("PRACTICAL RECOMMENDATION:\n")
cat("  For a math course research study:\n")
cat("  6-8 sessions × 40 teachers × 1 semester\n")
cat("  with full-transcript automated scoring (~30-60 exchanges/session)\n")
cat("  → ~180-480 exchanges per teacher\n")
cat("  → Teacher-level reliability ≈ 0.70-0.80+\n")

cat("\n=== ANALYSIS COMPLETE ===\n")