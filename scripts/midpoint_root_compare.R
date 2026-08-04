library(ape)
library(phangorn)

base <- "data/test_data/strain/test1/prediction/protal051/tree"
tree_files <- c(
  file.path(base, "s__Bacteroides_ovatus.partitioned.nwk"),
  file.path(base, "s__Bacteroides_ovatus.unpartitioned.nwk"),
  file.path(base, "s__Bacteroides_xylanisolvens.partitioned.nwk"),
  file.path(base, "s__Bacteroides_xylanisolvens.unpartitioned.nwk")
)

for (f in tree_files) {
  t <- read.tree(f)
  t_rooted <- midpoint(t)
  out <- sub("\\.nwk$", ".midpoint.nwk", f)
  write.tree(t_rooted, file = out)
  cat("Written:", out, "\n")
}
