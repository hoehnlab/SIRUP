# docs/build.R - rebuild the generated parts of the docs. Run from the package root:
#
#   Rscript docs/build.R              # everything
#   Rscript docs/build.R topics       # only man/*.Rd -> docs/topics/*.md
#   Rscript docs/build.R vignettes    # only vignettes -> docs/vignettes/*.md
#
# Needs: R, pandoc, the packages in Suggests that the vignettes use, and python3
# with PyYAML. In CI all of that comes from the docs Docker image.
# This file mirrors Dowser's docs/build.R (a thin entry point), but it does not use markr.

steps <- commandArgs(trailingOnly = TRUE)
if (length(steps) == 0) steps <- c("topics", "vignettes", "nav", "extras")

config <- yaml::read_yaml("docs/docs-config.yml")
source("docs/build_topics.R")
source("docs/render_vignettes.R")

if ("topics" %in% steps)     build_topics(config)
if ("vignettes" %in% steps)  render_vignettes(config)
# The text-only steps live in Python so they also run where R is absent.
if ("nav" %in% steps)    stopifnot(system2("python3", c("docs/docs_tools.py", "nav")) == 0)
if ("extras" %in% steps) stopifnot(system2("python3", c("docs/docs_tools.py", "extras")) == 0)
message("Done. Next: python3 docs/docs_tools.py check")
