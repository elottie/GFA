library(readr)
library(dplyr)

comparison_files <- list.files("ldsc_block", pattern = "set_\\d+_\\d+\\.csv", full.names = TRUE)
print(comparison_files)

# Helper to load trait IDs for a set (from comparison file)
get_trait_ids <- function(comparison_file, setnum) {
  comparisons <- read_csv(comparison_file)
  setnums <- as.numeric(regmatches(basename(comparison_file), gregexpr("\\d+", basename(comparison_file)))[[1]])
  if(setnum == setnums[1]) {
    unique(comparisons$trait1)
  } else {
    unique(comparisons$trait2)
  }
}

make_trait_df <- function(trait_ids, set_name) {
  data.frame(trait = trait_ids, value = runif(length(trait_ids)), set = set_name)
}

set_dfs <- list() # to hold only what's needed

for (i in seq_along(comparison_files)) {
  file <- comparison_files[i]
  setnums <- as.numeric(regmatches(basename(file), gregexpr("\\d+", basename(file)))[[1]])
  print(setnums)
  s1name <- paste0("s", setnums[1])
  s2name <- paste0("s", setnums[2])

  print(paste("s1name",s1name))
  print(paste("s2name",s2name))

  # Only create set dataframe if it doesn't already exist/nor just deleted
  if (!s1name %in% names(set_dfs)) {
    trait_ids_s1 <- get_trait_ids(file, setnums[1])
    set_dfs[[s1name]] <- make_trait_df(trait_ids_s1, s1name)
  }
  if (!s2name %in% names(set_dfs)) {
    trait_ids_s2 <- get_trait_ids(file, setnums[2])
    set_dfs[[s2name]] <- make_trait_df(trait_ids_s2, s2name)
  }

  comparisons <- read_csv(file)
  # head(comparisons) # Remove unless you want to print each time

  cat("\n--- Comparing traits in", s1name, "to traits in", s2name, "---\n")

  for (j in 1:nrow(comparisons)) {
    t1 <- comparisons$trait1[j]
    t2 <- comparisons$trait2[j]
    cat(sprintf("Comparing %s (from %s) to %s (from %s)\n", t1, s1name, t2, s2name))
  }

  # ---- TRIANGLE PROCESSING ----
  # For set_1_2: also process t1, t2
  # For set_1_3: also process t3
  # For set_1_4: also process t4
  triangle_dir <- "ldsc_block" # change if triangles are elsewhere

  if (setnums[2] == 2) { # For set_1_2, do t1 AND t2
    print('in setnums2')
    for (tri_num in 1:2) {
      tri_file <- file.path(triangle_dir, paste0("triangle_", tri_num, ".csv"))
      print(paste('tri_file:',tri_file))
      if (file.exists(tri_file)) {
        print('tri_file exists')
        triangle_comparisons <- read_csv(tri_file)
        cat(sprintf("\nProcessing triangle comparisons in %s:\n", basename(tri_file)))
        for (k in 1:nrow(triangle_comparisons)) {
          tri_t1 <- triangle_comparisons$trait1[k]
          tri_t2 <- triangle_comparisons$trait2[k]
          cat(sprintf("Triangle: Comparing %s to %s\n", tri_t1, tri_t2))
        }
      }
    }
  } else if (setnums[2] >= 3) { # For set_1_3, set_1_4, etc: just t3, t4, ...
    tri_file <- file.path(triangle_dir, paste0("triangle_", setnums[2], ".csv"))
    if (file.exists(tri_file)) {
      triangle_comparisons <- read_csv(tri_file)
      cat(sprintf("\nProcessing triangle comparisons in %s:\n", basename(tri_file)))
      for (k in 1:nrow(triangle_comparisons)) {
        tri_t1 <- triangle_comparisons$trait1[k]
        tri_t2 <- triangle_comparisons$trait2[k]
        cat(sprintf("Triangle: Comparing %s to %s\n", tri_t1, tri_t2))
      }
    }
  }

  # ---- END TRIANGLE PROCESSING ----

  # Drop s2 after comparison if not the last set (you could adapt this logic as needed)
  set_dfs[[s1name]] <- NULL # optional, if you only want to keep s2 (or s3...) for next round
  cat(sprintf("Dropped dataframe for %s\n", s2name))
}

cat("\n--- All comparisons (and triangles) complete! ---\n")
