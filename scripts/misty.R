
suppressPackageStartupMessages({
  library(Seurat)
  library(schard)
  library(mistyR)
  library(tidyverse)
})

## ------------------------------ Parameters -----------------------------------
h5ad       <- "Data/dummy1.h5ad"
sample_id  <- "dummy1"
juxta_thr  <- 200   
para_l     <- 500  
cv_folds   <- 10
seed       <- 42
out_dir    <- "misty"
fig_dir    <- "Figures"
visium_pitch_um <- 100
dir.create(fig_dir, showWarnings = FALSE)

## ------------------------------ Data -----------------------------------------
seurat_obj <- schard::h5ad2seurat_spatial(h5ad)
cn <- colnames(seurat_obj)

## Targets: one binary indicator per flow domain
seurat_obj$flow_cluster <- as.factor(seurat_obj$flow_cluster)
targets_df <- as.data.frame(model.matrix(~ 0 + flow_cluster, data = seurat_obj@meta.data))
colnames(targets_df) <- gsub("^flow_cluster", "", colnames(targets_df))
rownames(targets_df) <- cn

## Predictor: pathogen label (1 = P+, 0 = P-)
pathogen_raw <- as.character(seurat_obj$bacterial_labels)
pathogen_vec <- ifelse(grepl("\\+|1|true", tolower(pathogen_raw)), 1, 0)

## Spot coordinates
coords <- GetTissueCoordinates(seurat_obj)[, c("x", "y")]
xy <- as.matrix(coords)



## ------------------------------ Views ----------------------------------------
juxta_name <- paste0("juxtaview.", juxta_thr)
para_name  <- paste0("paraview.", para_l)

build_views <- function(pvec) {
  # `own` view: the label of the spot itself
  own_df <- data.frame(Pathogen = pvec, row.names = cn)
  # juxta / para views: built from a two-column table (a single column makes
  # add_juxtaview fail), then only the Pathogen columns are kept
  p2 <- data.frame(Pathogen = pvec, NonPathogen = 1 - pvec, row.names = cn)
  pv <- create_initial_view(p2) %>%
    add_juxtaview(positions = coords, neighbor.thr = juxta_thr) %>%
    add_paraview(positions = coords, l = para_l)
  list(
    own   = own_df,
    juxta = pv[[juxta_name]]$data %>% dplyr::select(dplyr::starts_with("Pathogen")),
    para  = pv[[para_name]]$data  %>% dplyr::select(dplyr::starts_with("Pathogen"))
  )
}

## View sets: own only; own + each scale; all views
view_sets <- list(own       = "own",
                  own_juxta = c("own", "juxta"),
                  own_para  = c("own", "para"),
                  all       = c("own", "juxta", "para"))

run_sets <- function(pvec, tag) {
  v <- build_views(pvec)
  purrr::map_dfr(names(view_sets), function(nm) {
    mv <- create_initial_view(targets_df)
    for (u in view_sets[[nm]]) mv <- mv %>% add_views(create_view(u, v[[u]], u))
    folder <- file.path(out_dir, tag, nm, sample_id)
    run_misty(mv, results.folder = folder,
              bypass.intra = TRUE, cv.folds = cv_folds, seed = seed)
    collect_results(folder)$improvements %>%
      filter(measure == "gain.R2") %>%
      transmute(tag = tag, set = nm, target, gain_R2 = value)
  })
}

## ------------------------ Spatial null: reflected maps -----------------------
## The pathogen map is reflected over the tissue and each reflected position is
## assigned the label of the nearest spot, so that the spatial structure of the
## map is preserved but decoupled from the tissue.
nearest_spot <- function(m) apply(m, 1, function(p)
  which.min(rowSums((xy - matrix(p, nrow(xy), 2, byrow = TRUE))^2)))
sx <- sum(range(xy[, 1])); sy <- sum(range(xy[, 2]))
flipx  <- cbind(sx - xy[, 1], xy[, 2])
flipy  <- cbind(xy[, 1], sy - xy[, 2])
rot180 <- cbind(sx - xy[, 1], sy - xy[, 2])

misty_res <- bind_rows(
  run_sets(pathogen_vec,                         "real"),
  run_sets(pathogen_vec[nearest_spot(flipx)],    "null_flipx"),
  run_sets(pathogen_vec[nearest_spot(flipy)],    "null_flipy"),
  run_sets(pathogen_vec[nearest_spot(rot180)],   "null_rot180"))

write_csv(misty_res, file.path(out_dir, "misty_gainR2_all_runs.csv"))

## ------------------------- Variance explained and increments -----------------
wide <- misty_res %>%
  pivot_wider(names_from = set, values_from = gain_R2) %>%
  mutate(own_R2  = own,
         d_local = own_juxta - own,    # increase when adding the local view
         d_wide  = own_para  - own,    # increase when adding the wide view
         d_all   = all       - own)    # increase when adding both views

real_tab <- wide %>% filter(tag == "real") %>%
  select(target, own_R2, d_local, d_wide, d_all)

null_max <- wide %>% filter(tag != "real") %>%
  group_by(target) %>%
  summarise(across(c(d_local, d_wide, d_all), max, .names = "null_max_{.col}"),
            .groups = "drop")

summary_tab <- real_tab %>%
  left_join(null_max, by = "target") %>%
  mutate(ratio_all_vs_null = d_all / null_max_d_all)   # real increase / largest null increase

cat("\n=== Variance explained (own) and increments over `own` (percentage points) ===\n")
print(summary_tab %>% mutate(across(where(is.numeric), ~ round(.x, 2))), n = Inf)
write_csv(summary_tab, file.path(out_dir, "misty_summary.csv"))

## ----------------------------- Direction of the association ------------------
v_real <- build_views(pathogen_vec)
cor_tab <- round(sapply(colnames(targets_df), function(dm) c(
  own   = cor(pathogen_vec,       targets_df[[dm]]),
  juxta = cor(v_real$juxta[[1]],  targets_df[[dm]]),
  para  = cor(v_real$para[[1]],   targets_df[[dm]]))), 3)
cat("\n=== Pearson correlation: domain indicator vs pathogen variable (rows: view) ===\n")
print(cor_tab)
write.csv(cor_tab, file.path(out_dir, "misty_correlations.csv"))

## ----------------------------------- Figure ----------------------------------
delta <- wide %>%
  transmute(tag, target, Local = d_local, Wide = d_wide) %>%
  pivot_longer(c(Local, Wide), names_to = "scale", values_to = "dR2")

juxta_lab <- sprintf("Juxtaview (%.0f um)", juxta_thr )
para_lab  <- sprintf("Paraview (%.0f um)",  para_l)

p <- ggplot() +
  geom_col(data = filter(delta, tag == "real"),
           aes(target, dR2, fill = scale),
           position = position_dodge(0.8), width = 0.7, colour = "black") +
  geom_point(data = filter(delta, tag != "real"),
             aes(target, dR2, group = scale),
             position = position_dodge(0.8), colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c(Local = "#eebb11", Wide = "#991111"),
                    labels = c(Local = juxta_lab, Wide = para_lab)) +
  labs(x = "", y = expression(Delta * R^2 ~ "(percentage points) over the spot label"),
       fill = "Neighbourhood scale") +
  theme_classic() + theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "mistyR_deltaR2.png"), p, dpi = 300, width = 6, height = 5)