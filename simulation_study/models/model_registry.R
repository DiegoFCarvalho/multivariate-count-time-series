# Sources every model file and builds `all_models`; the only file that names individual model files.

source("models/model_interface.R")
source("models/model_svines.R")
source("models/model_minar.R")
source("models/model_mdacp.R")
source("models/model_multiv_count_ar.R")

# Each model is one object used under the same name as DGP and as Est.Met.
# True parameters live in run_simulation.R.
minar1_poi    <- make_minar_model("minar1_poi", p = 1, family = "poisson")
minar1_negbin <- make_minar_model("minar1_negbin", p = 1, family = "negbin")
mdacp1        <- make_mdacp_model("mdacp1", r = 2)

# DGP only: has no est().
mcar1 <- make_mcar_model("mcar1", copula_family = "gaussian")

# The DGP-side theta is a hand-built svine_dist object; the Est.Met side gets its theta from est().
svine1 <- make_svines_model("svine1", p = 2)

all_models <- list(
  minar1_poi    = minar1_poi,
  minar1_negbin = minar1_negbin,
  mdacp1        = mdacp1,
  mcar1         = mcar1,
  svine1        = svine1
)
