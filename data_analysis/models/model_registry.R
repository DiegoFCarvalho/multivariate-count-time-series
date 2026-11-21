# Model objects used in the analysis.
source("models/model_interface.R")
source("models/model_svines.R")
source("models/model_mdacp.R")

svines_exact  <- make_svines_model("svines_exact", p = 2, margins = "exact_ecdf")
svines_scaled <- make_svines_model("svines_scaled", p = 2, margins = "scaled_ecdf")
mdacp1        <- make_mdacp_model("mdacp1", r = 2)

all_models <- list(svines_exact = svines_exact, svines_scaled = svines_scaled, mdacp1 = mdacp1)
