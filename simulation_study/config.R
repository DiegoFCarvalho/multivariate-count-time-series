# Tunable parameters for a simulation run.

get_config <- function() {
  list(
    n          = 500,     # Monte Carlo replications
    d          = 3,       # dimension of the count series
    T          = 100,     # series length for fitting and for forecast simulations
    T_burnin   = 200,     # observations simulated and discarded before the recorded DGP path
    K          = 100,     # simulated paths per fitted model
    lags       = 0:3,     # lags for the autocovariance and rank-correlation statistics
    seed       = 1,
    n_workers  = 20,      # parallel workers over replications
    est_timeout_seconds = 480,   # per-fit timeout; a longer fit is reported as "timed_out"
    output_dir = "results/",
    # mcar1 is DGP only (it has no est()); the other DGPs are also Est.Mets.
    dgp_names     = c("minar1_poi", "minar1_negbin", "mdacp1", "mcar1", "svine1"),
    est_met_names = c("minar1_poi", "minar1_negbin", "mdacp1", "svine1")  # fitted against every DGP
  )
}
