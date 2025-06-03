#' Biological antirheumatic agents in rheumatoid arthritis
#'
#' Samples from the joint distribution of the relative effects for each
#' active treatment (abatacept, anti-TNF, rituximab, tocilizumab) vs. placebo
#' estimated in a random effects Bayesian network meta-analysis with a binomial
#' likelihood and logit link. Three chains were run, the first 50,000 samples
#' were discarded as burn-in, and the subsequent 100,000 samples (a total of
#' 300,000 samples) were saved in this data frame. The outcome is American
#' College of Rheumatology 50% improvement (ACR50) response rate at 6 months.
#'
#' @format A data frame with 300,000 rows and 6 columns:
#' \describe{
#'   \item{`Chain`}{The index of the chain.}
#'   \item{`Placebo`}{The sampled relative effects of placebo vs. placebo
#'   (should be 0).}
#'   \item{`Abatacept`}{The sampled relative effects of abatacept vs. placebo.}
#'   \item{`Anti.TNF`}{The sampled relative effects of anti-TNF vs. placebo.}
#'   \item{`Rituximab`}{The sampled relative effects of rituximab vs. placebo.}
#'   \item{`Tocilizumab`}{The sampled relative effects of tocilizumab vs.
#'   placebo.}
#' }
#'
#' @source Salliot C, Finckh A, Katchamart W, Lu Y, Sun Y, Bombardier C, Keystone E. (2011). Indirect comparisons of the efficacy of biological antirheumatic agents in rheumatoid arthritis in patients with an inadequate response to conventional disease-modifying antirheumatic drugs or to an anti-tumour necrosis factor agent: a meta-analysis. Ann Rheum Dis 70(2):266-71. [doi:10.1136/ard.2010.132134]. Adopted from nmadb::readByID(501387)$data.
#'
#' data(dat_Salliott2010)
"dat_Salliott2010"
