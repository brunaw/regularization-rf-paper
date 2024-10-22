# Package loading -------
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(extrafont)

# Setting a few colors 
colfunc <- colorRampPalette(c("#efe1c6", "#3B9AB2"))
cols <- colfunc(9)

# Loading previous results files
# 1. Standard RF
std_rf <- read_rds("auxiliar_results/std_rf.rds")

# 2. Constant lambda_0 with rrf
df_rrf_scaled <- readRDS("auxiliar_results/df_rrf_scaled_new.Rds")

# 3. GRRF with rrf
df_grrf_scaled <- read_rds("auxiliar_results/df_grrf_scaled_new.rds")

# 4. Generalized regularization models 
df_mix <- read.table("auxiliar_results/df_mix.txt") %>% 
  filter(lambdas_df %in% c(0.1, 0.5, 0.9))

# Figure 1. Standard RF Results  -----------------
pdf("plots/regression/scaled_RF.pdf", height = 3.85,
    width = 7, paper='special') 
std_rf$df_rf %>% 
  group_by(mtry_df) %>% 
  dplyr::summarise_if(is.numeric, mean) %>% 
  ungroup() %>% 
  ggplot() +
  geom_point(aes(y = factor(pars), 
                 x = rmse, 
                 fill = factor(mtry_df)), 
             size = 4, shape = 23) + 
  geom_vline(aes(xintercept = mean(rmse)), colour = "red") +
  annotate(geom = "text", y = 1.42, x = 0.493, 
           label = paste0("Mean = ", round(mean(std_rf$df_rf$rmse), 3)),
           colour = "red", size = 6) + 
  scale_x_continuous(
    breaks = scales::pretty_breaks()) + 
  scale_fill_manual(values = cols)+
  theme_minimal(22) +
  labs(fill = expression(italic(mtry)), 
       y = "# Features", 
       x = expression("RMSE"[test])) +
  theme( 
    strip.placement = "outside",
    legend.title = element_text(color = "black", size = 16),
    legend.position = "bottom", 
    text = element_text(family = "Times"))
dev.off()

# Figure 2. Using only lambda_0 with the rrf and varying
# gamma in the GRRF ------

# First: merging the ranger and rrf results and wrangling data
df_wrangled <- df_grrf_scaled$df_grrf_scaled %>% 
  mutate(model = "Guided Regularization: rrf") %>% 
  bind_rows(df_rrf_scaled$df_rrf_scaled %>%
              mutate(model = "Constant Regularization: rrf") %>% 
              select(-type_depth)) %>% 
  group_by(model, mtry_df, lambdas_df) %>% 
  dplyr::summarise_if(is.numeric, mean) %>% 
  ungroup() %>% 
  mutate(mtry_df = as.integer(mtry_df), 
         model = as.factor(model)) 

max_rmse <- max(df_wrangled$rmse)
min_rmse <- min(df_wrangled$rmse)
seq_rmse <- round(seq(min_rmse, max_rmse, length.out = 12), 2)

# Tile plot for the RMSEs with the RRF and GRRF
p1 <- df_wrangled %>% 
  split(.$model) %>% 
  map({
    ~.x %>% 
      ggplot() + 
      facet_wrap(~model, nrow = 1) +
      geom_tile(aes(x = lambdas_df, 
                    y = factor(mtry_df), 
                    fill = rmse))  +
      scale_fill_gradient(high = "#FFA07A", low  = "#026440", 
                          limits = c(min_rmse, max_rmse), 
                          breaks = seq_rmse) +
      theme_minimal(22) +
      theme(panel.grid=element_blank(), 
            legend.position="top", 
            legend.text = element_text(angle = 90, hjust = 1), 
            axis.text.x = element_text(size = 14), 
            strip.text  = element_text(colour = "black", 
                                       size = 20, face = "bold"),
            text = element_text(family = "Times")) +
      labs(x = 
             expression("Constant "~lambda[i]~"(RRF)"), 
           y = expression(italic(mtry)), 
           fill = expression("RMSE"[test])) +
      coord_cartesian(expand=FALSE) +
      guides(fill = guide_colourbar(barwidth = 18, label.position = "bottom", 
                                    title.position = "top", title.hjust = 0.5))
  })

max_pars <- max(df_wrangled$pars)
min_pars <- min(df_wrangled$pars)

# Tile plot for the the number of parameters 
p2 <- df_wrangled %>% 
  split(.$model) %>% 
  map({
    ~.x %>% 
      ggplot() + 
      facet_wrap(~model, nrow = 1) +
      geom_tile(aes(x = lambdas_df, 
                    y = factor(mtry_df), 
                    fill = pars))  +
      scale_fill_gradient(low = "#f8eed1", high  = "#3B9AB2", 
                          limits  = c(min_pars, max_pars), 
                          breaks = c(0, 1, 10, 25, 50, 75, 100, 
                                     125, 150, 175, 200, 
                                     225, 250, 256)) +  
      theme_minimal(22) +
      theme(panel.grid = element_blank(), 
            legend.position="top", 
            legend.text = element_text(angle = 90, hjust = 1), 
            axis.text.x = element_text(size = 14), 
            strip.text  = element_text(colour = "black", 
                                       size = 20, face = "bold"),
            text = element_text(family="Times")
      ) +
      labs(x = 
             expression("Constant "~lambda[i]), 
           y = expression(italic(mtry)), 
           fill = "# Features") +
      coord_cartesian(expand=FALSE) +
      guides(fill = guide_colourbar(barwidth = 18, label.position="bottom", 
                                    title.position = "top", title.hjust = 0.5))
  })

# Saving results in one pdf file -------
pdf("plots/regression/tile_rmse_scaled.pdf", height = 12.8,
    width = 9.5, paper='special')   
(p1[[1]]) + (p2[[1]] + labs(y = "")) + 
  (p1[[2]] + guides(fill = FALSE) + 
     labs(x = expression("Varying "~gamma~" (GRRF)"))) +
  (p2[[2]] + guides(fill = FALSE) +  
     labs(y = "", 
          x = expression("Varying "~gamma~" (GRRF)"))) + plot_layout(nrow = 2)
dev.off()

# Generalized regularization plots --------------------------
levels(df_mix$model) <- c("Boosted[RF],", "Boosted[SVM],", "Correlation,")

max_rmse <- round(max(df_mix$rmse), 3)
min_rmse <- round(min(df_mix$rmse), 3)

p <- df_mix %>% 
  group_by(model, mtry_df, lambdas_df, gamma_df) %>% 
  summarise(mean_pars = mean(pars), 
            max_pars = max(pars), 
            min_pars = min(pars), 
            mean_rmse = mean(rmse), 
            max_rmse_df = max(rmse), 
            min_rmse_df = min(rmse)) %>% 
  ungroup() %>% 
  mutate(lambdas_df_name =
           paste0("lambda_0: ", lambdas_df)) %>% 
  ggplot(aes(mean_pars, mean_rmse)) +
  geom_hline(aes(yintercept = 0.5, linetype = " "),
             color = "red", show.legend = TRUE) +
  geom_errorbar(aes(ymin = min_rmse_df, ymax = max_rmse_df),
                position = position_dodge(.9), 
                width = 0.1) +
  geom_point(aes(fill = factor(mtry_df)), 
             size = 3, shape = 23) +
  facet_grid(gamma_df ~  model + lambdas_df, 
             #scales = "free",
             switch="y",
             labeller = 
               label_bquote(gamma==.(round(gamma_df, 3)),
                            atop(.(as.character(model)), 
                                 lambda[0]==.(lambdas_df)))) +
  scale_fill_manual(values = cols)+
  scale_x_continuous(
    trans = "log", 
    breaks = c(0, 1, 3, 5, 10, 25, 50, 125, 250),
    expand = c(0.1, 0.1)
  ) +
  ylim(min_rmse, max_rmse) + 
  labs(y = expression("RMSE"[test]),
       x = expression("Number of selected features"),  
       fill = expression(italic(mtry)), 
       shape = expression(gamma)) +
  theme_minimal(23) +
  theme( 
    text = element_text(family="Times"), 
    legend.title = element_text(color = "black", size = 18),
    legend.position="bottom", 
    strip.placement = "outside", 
    axis.text.x = element_text(
      size = 12, 
      angle = 90, hjust = 1, 
      color = "black"),
    axis.title.x = element_text(margin = unit(c(10, 0, 0, 0), "mm")),
    axis.text.x.bottom = element_text(size = 16)) +
  scale_linetype_discrete(name = expression("RMSE"[test]~"= 0.5")) +
  guides(fill = guide_legend(override.aes = list(linetype = 0)))  

pdf("plots/regression/mixture_all.pdf", 
    height = 10, width = 19, paper='special') 
p
dev.off()
#----------------------------------------------------------------
# End of code
#-------------------------------------------------------------------