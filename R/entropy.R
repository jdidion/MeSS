#!/usr/bin/env Rscript --vanilla

# Simplification and speed-up of SMART-BS-Seq by Liu et al. 
# (https://pypi.python.org/pypi/SMART-BS-Seq).
#
# Rather than starting from multiple .wig files, this script starts from a
# matrix of CpG x sample. The output is a data.frame with the following columns:
#
# * A specificity value [0-1] for each CpG, which Liu interpret as follows:
# < 0.5 == low specificity; 0.5-0.75 == intermediate specificity;
# > 0.75 == high specificity.
# * Median methylation, computed by one-step Tukey biweight.
# * Distance: geometric difference from previous row).
# * DiffEntropy: entropy of difference from previous row).
# * MeSS.*: Methylation Specificity Score for each tissue.
# * Optionally (see below):
#   * p.*: Significance of deviation from median, for each tissue.
#   * padj.*: Adjusted p-value for each tissue.
#
# Minimum beta range threshold: When there are many beta values very close to
# each other (this happens frequently at the extremes), but one or a few betas
# are slightly different, the standard entropy calculation might yield a high
# specificity. While technically correct, these scores cause difficulties for
# interpretation. Liu et al. got around this problem by enforcing a minimum
# difference of 1.0 (1%) between each beta and the median (biweight). Note that
# this is not mentioned in their paper; one must read their code to find it out.
# It is difficult to known why they chose 1.0, but we suspect it was an
# arbitrary decision. This approach leads to a strange distribution of MeSSes,
# with a large spike toward the low end (the value of the spike is determined
# by the sample size; for 21 samples it was between 0.18-0.19). We instead
# recommend determining this value empirically from the distribution of per-CpG
# beta ranges, or by the resolution of the assay (1 / average read depth for
# WGBS). We chose to go with the first option and selected 3.0 (3%) as our
# threshold; YMMV.
#
# P-values: To determine the specificity of methylation in each tissue,
# Liu et al. use the per-sample weights from the Tukey Biweight calculation to
# identify samples that belong to the "background" distribution and then run
# a one-sample t-test for each non-background sample. They do not perform
# multiple test correction. We assert that p-values are problematic for
# determining tissue specificity for the reasons described above. Nevertheless,
# we show that, after removing CpGs with beta value range below the empirically
# determined threshold, there is a high degree of correlation between MeSS and
# p-value. We also show that computing p-values based on betas is problematic
# due to the non-normality of the data, and we show that first converting betas
# to M values leads to higher correlation with MeSS. We also recommend computing
# p-values using a 'leave-one-out' strategy (i.e. treat one sample as the mean
# and use the remaining samples as the distribution). We provide an option
# (disabled by default) to compute p-values.

# TODO: Adapt conditional entropy/information gain (e.g.
# http://bioinformatics.oxfordjournals.org/content/32/20/3107.full.pdf)
# to identify the sites that are the most important determinants of tissue
# specificity.

library(GenomicRanges)
library(plyr)
library(matrixStats)

cpg_specificity <- function(
        pos, meth, min.beta.range=3.0, log.base=NULL, max.cores=22,
        validate=TRUE, pairwise=TRUE, pvalue_function=NULL, ranks=TRUE, ...) {
    # Compute methylation specificity for each CpG.
    # Args:
    #   pos: GenomicRanges object with N CpGs.
    #   meth: matrix of N x M samples. Contains methylation values between 0-1
    #       (missing values are not allowed).
    #   min.beta.range: The minimum difference (in percent) between the max and
    #       min value for each CpG.
    #   max.cores: maximum number of parallel execution threads
    #   validate: whether to pre-process methylation data
    #   pairwise: Whether to calcluate distance and entropy of difference for
    #       consecutive pairs of CpGs.
    #   pvalue_function: Function to use for computing p-values.
    #   ...: Additional arguments to p-value function.
    #
    # Returns N x (6+M) data.frame where the columns are chromosome, start, end,
    # specificity, weight (calculated using Tukey's Biweight function), and
    # whether the beta range is below 'min.beta.range'. If `pairwise==TRUE`,
    # there are two additional columns, entropy of difference from previous CpG,
    # and Euclidian distance from previous CpG. Finally, there are M MeSS
    # scores. If `pvalue_function` is non-NULL, pvalues are computed and both
    # raw and adjusted pvalues are reported (an additional 2M columns). If
    # `ranks==TRUE`, MeSSes are converted to ranks (an additional M columns).
    
    # ensure meth is a matrix of the proper form
    if (validate) {
        if (!is.matrix(meth)) {
            meth <- as.matrix(meth)
        }
        # Remove rows with missing values
        w <- apply(meth, 1, function(x) any(is.na(x)))
        if (any(w)) {
            pos <- pos[!w,]
            meth <- meth[!w,]
        }
        # Make sure values are in the range [0-100]
        if (all(meth <= 1)) {
            meth <- meth * 100
        }
        meth[meth < 0] <- 0
        meth[meth > 100] <- 100
    }
    
    # Only look at autosomes
    chroms <- as.vector(seqnames(pos))
    keep <- chroms %in% paste0('chr', 1:22)
    if (any(!keep)) {
        chroms <- chroms[keep]
        pos <- pos[keep,]
        meth <- meth[keep,]
    }
    
    sample_names <- colnames(meth)
    num_samples <- length(sample_names)
    
    # Here we follow Liu et al and use M as the base for the logarithm.
    # Ultimately, it doesn't matter what the base is; using a smaller value
    # will generate a wider range of entropy values, but normalization at the
    # end should lead to roughly (exactly?) the same values at the end.
    if (is.null(log.base)) {
        log.base <- num_samples
    }
    
    # If we want pairwise info, the data has to be processed per-chromosome.
    # This is somewhat wasteful because chr22 will take much less time than
    # chr1 to process, so the faster threads will be idle. Otherwise we can
    # parallize over the whole matrix to speed things up.
    
    if (pairwise) {
        cores <- min(length(unique(chroms)), max.cores)
    }
    else {
        cores <- max.cores
    }
    
    parallel <- FALSE
    if (cores > 1) {
        doMC::registerDoMC(cores=cores)
        parallel <- TRUE
    }
    
    if (pairwise) {
        # Run across chromsomes
        df <- ldply(split.data.frame(meth, chroms), chunk_cpg_specificity,
            num_samples, min.beta.range, log.base,
            pairwise=TRUE, pvalue_function=pvalue_function, ...,
            .parallel=parallel)[,-1]
    }
    else if (parallel) {
        # Run across entire data frame
        num_cpg <- nrow(meth)
        chunks <- ceiling(num_cpg / cores)
        chunk_size <- ceiling(num_cpg / chunks)
        df <- ldply(
            split.data.frame(meth, rep(1:chunks, each=chunk_size, length.out=num_cpg)),
            chunk_cpg_specificity,
            num_samples, min.beta.range, log.base,
            pairwise=FALSE, pvalue_function=pvalue_function, ...,
            .parallel=TRUE)[,-1]
    }
    else {
        df <- chunk_cpg_specificity(
            meth, num_samples, min.beta.range, log.base, pairwise=FALSE)
    }
    
    rm(meth)
    
    # Normalize MeSS
    mess.idx <- grep('MeSS.', colnames(df))
    mess <- df[,mess.idx]
    max.mess <- apply(mess, 2, max)
    df[,mess.idx] <- 1 - sweep(mess, 2, FUN='/', max.mess)
    rm(mess, mess.idx, max.mess)
    
    if (ranks) {
        # Convert MeSS to ranks
        ranks <- adply(df[,mess.idx], 2, function(x) rank(-x, ties.method='random'))
        ranks <- t(ranks[,-1])
        rownames(ranks) <- NULL
        df <- data.frame(df, rank=ranks)
        rm(ranks)
    }
    
    if (!is.null(pvalue_function)) {
        # Perform multiple-test correction. This probably isn't the right way
        # to do this, but Liu et al don't offer guidance as to the correct
        # procedure.
        p.idx <- grep('p.', colnames(df), fixed=TRUE)
        p <- as.matrix(df[,p.idx])
        padj <- matrix(NA, nrow=nrow(p), ncol=ncol(p))
        colnames(padj) <- sample_names
        padj[!is.na(p)] <- p.adjust(p[!is.na(p)])
        df <- data.frame(df, padj=padj)
    }
    
    # add coordinates to df and return
    data.frame(pos, df)[,-c(4,5)]
}

chunk_cpg_specificity <- function(meth, num_samples, ..., pairwise=TRUE, pvalue_function=NULL) {
    # Calculate specificity scores for a single chromosome.
    # Args:
    #   meth: the CpG x sample methylation matrix
    #   num_samples: Number of samples; usually equals the number of columns
    #       in 'meth'
    #   min.meth.range: Minimum difference between max and min beta for each CpG
    
    # compute weights, Tukey biweights, aned entropies
    result <- calc_entropy_matrix(meth, num_samples, ..., pvalue_function=pvalue_function)

    # note that MeSS requires a global normalization
    df <- data.frame(
        Specificity=1-result$entropy1, TukeyBiweight=result$biweights,
        RangeBelowCutoff=result$thresholded, MeSS=result$Q)

    if (pairwise) {
        # calculate pairwise row distance
        diff <- meth[-nrow(meth),] - meth[-1,]
        rm(meth)
        # compute pairwise distance and entropy
        df$Distance <- c(0, sqrt(rowSums(diff^2)))
        df$DiffEntropy <- c(1, calc_entropy_matrix(diff, num_samples, ..., entropy.only=TRUE))
    }
    
    if (!is.null(pvalue_function)) {
        df$p=result$pvalues
    }
    
    df
}

# This is a faster method of computing the weights and entropy values
# for the entire matrix. See the original row-wise functions below
# for more details.
calc_entropy_matrix <- function(
        meth, num_samples, min.beta.range, log.base, entropy.only=FALSE,
        err=0.0001, const=5, pvalue_function=NULL, ...) {
    # eq 2
    med_diff <- meth - rowMedians(meth, na.rm=TRUE)
    u <- med_diff / ((const * rowMedians(abs(med_diff), na.rm=TRUE)) + err)
    
    # eq 3
    weights <- (1 - u^2)^2
    weights[is.na(weights) | weights < 0 | abs(u) > 1] <- 0
    
    # eq 4
    biweights <- rowSums(weights * meth) / rowSums(weights)
    biweights[is.na(biweights) | is.nan(biweights) | is.infinite(biweights)] <- 0
    
    # Here we deviate from BS-SMART-Seq by using an empirically-determined value
    # for the minimum m_prime rather than 1.0. This improves the distribution of
    # MeSSes and still avoids divide-by-zero and log-of-0 errors. We also flag
    # CpGs for which the range between min and max < min.beta.range.
    
    # eq #5
    m_prime <- abs(meth - biweights)
    m_prime[is.na(m_prime) | m_prime < min.beta.range] <- min.beta.range
    
    # eq #6
    m_rel <- m_prime / rowSums(m_prime) # p
    log_m_rel <- log(m_rel, log.base)
    entropy0 <- -rowSums(m_rel * log_m_rel) # H
    
    # eq #7
    beta_ranges <- rowMaxs(meth, na.rm=TRUE) - rowMins(meth, na.rm=TRUE)
    entropy_weights <- 1.0 - (beta_ranges / 100.0)
    
    # eq #8
    entropy1 <- pmin(pmin(entropy0, 1.0) * entropy_weights, 1.0)
    
    # if only entropy is desired, return the weighted entropy
    if (entropy.only) {
        entropy1
    }
    else {
        # otherwise compute specificity and p-values and return all the data
        result <- list(meth=m_rel, weights=weights, biweights=biweights,
             entropy0=entropy0, entropy1=entropy1,
             thresholded=beta_ranges < min.beta.range,
             Q=entropy1-log_m_rel)
        if (!is.null(pvalue_function)) {
            result$pvalues=pvalue_function(meth, weights, num_samples, ...)
        }
        result
    }
}

# This is the method used by Liu et al. to identify tissues that exhibit
# significant specificity at each t-test.
calc_pvalues_sparse_ttest_matrix <- function(
        meth, weights, num_samples,
        values=c(0, 0.001, 0.005, -0.001, -0.005),
        probs=c(1/6, 1/6, 1/6, 2/6, 1/6), min.pvalue=1E-100,
        parallel=FALSE) {
    background <- meth + matrix(
        sample(values, length(meth), TRUE, probs),
        nrow(meth), ncol(meth))
    test <- weights < 0.5
    background[test] <- NA
    num_bg <- rowSums(test)
    test_rows <- which(
        num_bg >= 2 &
        num_bg < num_samples &
        sqrt(rowVars(background, na.rm=T) / num_samples) > 0)
    test_cells <- which(test, arr.ind=T)
    test_cells <- test_cells[!(test_cells[,1] %in% test_rows),]
    pvalues <- matrix(NA, nrow(meth), ncol(meth))
    pvalues[test_cells] <- aaply(test_cells, 1, function(x) {
        t.test(background[x[1],], mu=meth[x], paired=FALSE)$p.value
    }, .parallel=parallel)
    pvalues[!is.na(pvalues) & pvalues < min.pvalue] <- min.pvalue
    colnames(pvalues) <- names(meth)
    pvalues
}

# This is a variant of the Liu et al method that tests all values against
# the background, where the background consists of all values with corresponding
# weights < 0.5 and excluding the value being tested.
calc_pvalues_all_ttest_matrix <- function(
        meth, weights, num_samples, min.pvalue=1E-100, parallel=FALSE) {
    background <- meth
    colnames(background) <- paste0(colnames(background), '.bg')
    background[weights < 0.5] <- NA
    ix <- seq(num_samples+1, 2*num_samples)
    pvalues <- aaply(cbind(meth=meth, bg=background), 1, function(x) {
        bg <- x[ix]
        if (sum(!is.na(bg)) < 3) {
            return(rep(NA, num_samples))
        }
        sapply(1:num_samples, function(i)
            t.test(bg[-i], mu=x[i], paired=FALSE)$p.value)
    }, .parallel=parallel)
    pvalues[!is.na(pvalues) & pvalues < min.pvalue] <- min.pvalue
    colnames(pvalues) <- colnames(meth)
    pvalues
}

# P-value calculations using M-values

beta_to_M <- function(meth, alpha) {
    meth01 <- meth / 100
    log2((meth01 + alpha) / (1 - meth01 + alpha))
}

calc_pvalues_ttest_matrix_mvalue <- function(meth, ..., alpha=1E-5) {
    calc_pvalues_ttest_matrix(beta_to_M(meth, alpha), ...)
}

calc_pvalues_all_ttest_matrix_mvalue <- function(meth, ..., alpha=1E-5) {
    calc_pvalues_all_ttest_matrix(beta_to_M(meth, alpha), ...)
}

# Command-line interface
parser <- argparse::ArgumentParser()
parser$add_argument("-t", "--threads", type="integer", default=1)
parser$add_argument("-o", "--outfile", default="specificity.RData")
parser$add_argument("-w", "--pairwise", action="store_true", default=FALSE)
parser$add_argument("-k", "--ranks", action="store_true", default=FALSE)
parser$add_argument("-p", "--pvalues", action="store_true", default=FALSE)
parser$add_argument("-d", "--data", default=NULL)
parser$add_argument("-m", "--meth", default=NULL)
parser$add_argument("-c", "--coords", default=NULL)
parser$add_argument("-r", "--min-beta-range", type="double", default=3.0)
parser$add_argument("-V", "--novalidate", action="store_true", default=FALSE)
args <- parser$parse_args()

if (!is.null(args$data)) {
    vars <- load(args$data)
} else {
    vars <- c(load(args$coords), load(args$meth))
}
if (vars[1] != 'pos') {
    pos <- get(vars[1])
    rm(vars[1])
}
if (vars[2] != 'meth') {
    meth <- get(vars[2])
    rm(vars[2])
}
result <- cpg_specificity(
    pos, meth, min.beta.range=args$min_beta_range, max.cores=args$threads,
    validate=!args$novalidate, pairwise=args$pairwise, pvalues=ifelse(
        args$pvalues, calc_pvalues_all_ttest_matrix_mvalue, NULL),
    ranks=args$ranks)
save(result, file=outfile)

###
# These are the slower row-wise functions. I'm keeping them here because
# they make the computation steps a bit more clear than the matrix functions.
###

# Compute the Tukey Biweight for a vector of methylation values.
# This calculation is derived from equations 2-4 in Zhang et al.
# Args:
#    err: small value used to avoid zero values in the denominator
#    const: tuning constant (Liu et al recommend 5)
tukey_biweight <- function(meth, err=0.0001, const=5) {
    # eq #2
    med_meth <- median(meth)  # M_r
    med_diff <- meth - med_meth
    med_med_diff <- median(abs(med_diff))  # S_r
    center_dist <- (const * med_med_diff) + err  # denominator
    
    # eq #3
    weights <- rep(0, length(meth))
    names(weights) <- names(meth)
    u <- med_diff / center_dist
    w <- abs(u) <= 1
    biweight <- 0
    if (any(w)) {
        weights[w] <- (1 - u[w]^2)^2
        if (any(weights > 0)) {
            # eq #4
            biweight <- sum(weights * meth) / sum(weights)
        }
    }
    list(biweight=biweight, weights=weights)
}

# Calcuate entropy for a vector of methylation values.
# This calculation is derived from equations 5-8 of Zhang et al.,
# with the modifications from Liu et al: 1) the log base in eq #6
# is the number of samples rather than 2, and 2) eq #7 is changed to
# w_r = 1 - ((max(m) - min(m)) / 100).
# Args:
#   meth: methylation matrix
#   num_samples: number of columns in the matrix
#   entropy.only: whether to only calculate the entropy value, or to also
#                 calculate p-values
#   err, const: Paremeters for Tukey biweight calculation.
#   values, probs, min.pvalue: Paramemters for p-value calculation (see
#                              below)
#
# P-values are computed by 1) identifying the tissues in which this CpG is
# most likely to have specific methylation, as determined by having
# weight > 0.5 from the Tukey biweight procedure; 2) adding a small amount
# of noise to the background values (sampled from 'values' with
# probabilities 'probs'); and 3) performing one-sample t-tests for each
# non-background tissue.
#
# Returns: If entropy.only=TRUE, returns a list (H0, H1, p), where H0 is
# the unweighted entropy, H1 is the weighted entropy, and p is the relative
# methylation signal for each sample. Otherwise,
# returns a 1-row data.frame of specificity (1-entropy), the Tukey biweight,
# and num_samples P values.
calc_entropy <- function(meth, num_samples) {
    tb <- tukey_biweight(meth)
    # eq #5
    m_prime <- abs(meth - tb$biweight)
    # eq #6
    m_sum <- sum(m_prime)
    if (m_sum == 0) {
        m_rel <- rep(1, num_samples)
    } else {
        m_rel <- m_prime / m_sum # p
    }
    log_m_rel <- log(m_rel, num_samples)
    entropy0 <- -sum(m_rel * log_m_rel) # H
    # eq #7
    weight <- 1.0 - ((max(meth) - min(meth)) / 100.0)
    # eq #8
    entropy1 <- min(min(entropy0, 1.0) * weight, 1.0)
    
    # If we only care about the entropy, return here
    return(list(p=m_rel, H0=entropy0, H1=entropy1, TB=tb, Q=entropy1-log_m_rel))
}

# Outlier tests
# Want to calculate on normalized (Z-score) M values (log2(b/(1-b)))

## Identify outliers iteratively
#outliers_iterative <- function() {
    # Use any of these functions to identify an outlier and remove it
    # iteratively until the p-value is no longer significant
#    dixon.test()
#    grubbs.test()
#    chisq.out.test()
#}

## GESD predicts how many outliers are in the data, and ranks observations by
## probability of outlier status.
#outliers_gesd <- function() {
#    gesd()
#}

#' Generalized Extreme Studentized Deviate (GESD) test.
#'
#' obs: Numeric vector of observations.
#' alpha: Significance level.
#' value.zscore: Whether observations are already z-normalized.
#' r: Upper bound of number of outliers.
#'
#' Value: List with two elements:
#'   total: total number of outliers
#'   ranks: outlier ranks for corresponding observations.
#'
#' From https://github.com/raunakms/GESD
gesd <- function(obs, alpha, value.zscore, r=NA) {
    #### Define and Declare Variables ####
    n <- length(obs)
    if(is.na(r)){r <- floor(n/2)} # by default, set upper bound on number of outliers 'r' to 1/2 sample size
    R <- numeric(length=r) # test statistics for 'r' outliers
    lambda <- numeric(length=r) # critical values for 'r' outliers
    outlier_ind <- numeric(length=r) # removed outlier observation values
    outlier_val <- numeric(length=r) # removed outlier observation values
    m <- 0 # number of outliers
    obs_new <- obs # temporary observation values
    
    #### Find outliers ####
    for(i in 1:r){
    
        #### Compute test statistic ####
        if(value.zscore) {
            z <- abs(obs_new) # If Z-score is alrealy computed
        } else {
            z <- abs(obs_new - mean(obs_new))/sd(obs_new) # Z-scores
        }
        
        max_ind <- which(z==max(z),arr.ind=T)[1] # in case of ties, return first one
        R[i] <- z[max_ind] # max Z-score
        outlier_val[i] <- obs_new[max_ind] # removed outlier observation values
        outlier_ind[i] <- which(obs_new[max_ind] == obs, arr.ind=T)[1] # index of removed outlier observation values
        obs_new <- obs_new[-max_ind] # remove observation that maximizes |x_i - x_mean|
        
        #### Compute critical values ####
        ##n_obs <- length(obs_new) # subset sample size
        p <- 1 - alpha/(2*(n-i+1)) # probability
        t_pv <- qt(p,df=(n-i-1)) # Critical value from Student's t distribution
        lambda[i] <- ((n-i)*t_pv) / (sqrt((n-i-1+t_pv^2)*(n-i+1)))
        
        #### Find exact number of outliers: largest 'i' such that R_i > lambda_i ####
        # print(c(i, R[i], lambda[i]))
        # try ((R[i] > lambda[i]), finally <- print(c(i, R[i], lambda[i]))
        # )
        if(!is.na(R[i]) & !is.na(lambda[i])) { # qt can produce NaNs
            if (R[i] > lambda[i]) {
                m <- i
            }
            # else {
                # break
            # }
            
        }
    }

    vals <- data.frame(NumOutliers=1:r, TestStatistic=R, CriticalValue=lambda)
    # print(vals)
    
    #### Rank outlier observations ####
    outlier_rank <- numeric(length=n)
    if (m > 0) {
        for (i in 1:m) {
            # outlier_rank[outlier_ind[i]] <- i
            outlier_rank[which(obs==outlier_val[i])] <- i
        }
    }
    names(outlier_rank) <- names(obs)
    num_outliers <- sum(outlier_rank != 0) #14 and 25 missing
    res <- c(num_outliers, outlier_rank)
    list(total=num_outliers, ranks=outlier_rank)
}

# Alternative measures of inequality (mostly drawn from econometrics)

# library(IC2)
#
# atkinson.index <- function(meth, eps=0) {
#     tb <- tukey_biweight(meth)
#     calcAtkinson(meth, tb$weights, eps)
# }
#
# generalized.entropy <- function(meth, alpha=1) {
#     tb <- tukey_biweight(meth)
#     calcGEI(meth, tb$weights, alpha)
# }
