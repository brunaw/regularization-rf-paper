# Package loading -------
library(tidyverse)
library(patchwork)
library(RColorBrewer)

# Setting a few colors 
colfunc <- colorRampPalette(c("#efe1c6", "#3B9AB2"))
cols <- colfunc(9)

# Loading previous results files
# Standard RF
std_rf <- read_rds("auxiliar_results/std_rf.rds")

# Constant lambda_0: for ranger and the rrf
df_ct <- read.table("auxiliar_results/df_ct.txt")
df_rrf_scaled <- readRDS("auxiliar_results/df_rrf_scaled.Rds")

# Mixture models 
df_mix <- read.table("auxiliar_results/df_mix.txt")

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
       y = "# Variables", 
       x = expression("RMSE"[test])) +
  theme( 
    strip.placement = "outside",
    legend.title = element_text(color = "black", size = 16),
    legend.position = "bottom") 
dev.off()

# Figure 2. Using only lambda_0 with ranger and the rrf ------
# First: merging the ranger and rrf results and wrangling data
df_wrangled <- df_ct %>% 
  mutate(model = "Generalized Regularization") %>% 
  select(-type_depth)  %>% 
  bind_rows(df_rrf_scaled$df_rrf_scaled %>% mutate(model = "rrf") %>% 
              select(-type_depth)) %>% 
  group_by(model, mtry_df, lambdas_df) %>% 
  dplyr::summarise_if(is.numeric, mean) %>% 
  ungroup() %>% 
  mutate(mtry_df = as.integer(mtry_df), 
         model = as.factor(model)) 

max_rmse <- max(df_wrangled$rmse)
min_rmse <- min(df_wrangled$rmse)

# Tile plot for the RMSEs with the RRF
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
                          limits = c(min_rmse, max_rmse)) +
      theme_minimal(20) +
      
      theme(panel.grid=element_blank(), 
            legend.position="top", 
            legend.text = element_text(angle = 90, hjust = 1), 
            axis.text.x = element_text(size = 12), 
            strip.text  = element_text(colour = "black", 
                                       size = 20, face = "bold")) +
      labs(x = 
             expression("Constant "~lambda[i]), 
           y = expression(italic(mtry)), 
           fill = expression("RMSE"[test])) +
      coord_cartesian(expand=FALSE)
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
                          limits  = c(min_pars, max_pars))+
      theme_minimal(20) +
      theme(panel.grid=element_blank(), 
            legend.position="top", 
            legend.text = element_text(angle = 90, hjust = 1), 
            axis.text.x = element_text(size = 12), 
            strip.text  = element_text(colour = "black", 
                                       size = 20, face = "bold")
      ) +
      labs(x = 
             expression("Constant "~lambda[i]), 
           y = expression(italic(mtry)), 
           fill = "# Variables") +
      coord_cartesian(expand=FALSE)
  })

# Saving results in one pdf file -------
pdf("plots/regression/tile_rmse_scaled.pdf", height = 12.8,
    width = 9.5, paper='special')   
(p1[[1]]) + (p2[[1]] + labs(y = "") ) + 
  (p1[[2]] + guides(fill = FALSE)) +
  
  (p2[[2]] + guides(fill = FALSE) +  labs(y = "")) + plot_layout(nrow = 2)
dev.off()


# Mixed models plots --------------------------
max_rmse <- 0.615
min_rmse <- 0.43

# Figure 3. Results using g(x_i) = corr(y, x_i)
p <- df_mix %>% 
  filter(model == "correlation") %>% 
  group_by(mtry_df, lambdas_df, gamma_df) %>% 
  summarise(mean_pars = mean(pars), 
            max_pars = max(pars), 
            min_pars = min(pars), 
            mean_rmse = mean(rmse), 
            max_rmse_df = max(rmse), 
            min_rmse_df = min(rmse)) %>% 
  ungroup() %>% 
  mutate(lambdas_df_name =
           paste0("lambda_0: ", lambdas_df)) %>% 
  ggplot(aes(log(mean_pars), mean_rmse)) +
  geom_hline(aes(yintercept = 0.5, linetype = " "),
             color = "red", show.legend = TRUE) +
  geom_errorbar(aes(ymin = min_rmse_df, ymax = max_rmse_df),
                position = position_dodge(.9), 
                width = 0.1) +
  geom_point(aes(fill = factor(mtry_df)), 
             size = 3, shape = 23) +
  facet_grid(gamma_df~lambdas_df, 
             scales = "free",
             switch="y", 
             labeller = 
               label_bquote(gamma==.(round(gamma_df, 3)),
                            lambda[0]==.(lambdas_df))) +
  scale_fill_manual(values = cols)+
  scale_x_continuous(
    limits = c(0, 5.53), 
    breaks = seq(0, 5.5, by = 0.5), 
    expand = c(0.02, 0.02)) +
  ylim(min_rmse, max_rmse) + 
  labs(y = expression("RMSE"[test]),
       x = expression(italic(log)~" number of selected variables"), 
       fill = expression(italic(mtry)), 
       shape = expression(gamma)) +
  theme_minimal(22) +
  theme( 
    legend.title = element_text(color = "black", size = 16),
    legend.position="bottom", 
    strip.placement = "outside", 
    axis.text.x = element_text(
      angle = 90, hjust = 1, 
      color = "black"),
    axis.title.x = element_text(margin = unit(c(10, 0, 0, 0), "mm")),
    #legend.text = element_text(angle = 90, hjust = 1, 
    #                           color = "black", size = 12), 
    axis.text.x.bottom = element_text(size = 12)) +
  scale_linetype_discrete(name = expression("RMSE"[test]~"= 0.5")) +
  guides(fill = guide_legend(override.aes = list(linetype = 0)))

pdf("plots/regression/mixture_corr_pars_rmse.pdf", 
    height = 11.5, width = 12, paper='special') 
p
dev.off()
#----------------------------------------------------------------------
# Figure 4. Results using g(x_i) = Boosted_{RF}

p <- df_mix %>% 
  filter(model == "boosted_rf") %>% 
  group_by(mtry_df, lambdas_df, gamma_df) %>% 
  summarise(mean_pars = mean(pars), 
            max_pars = max(pars), 
            min_pars = min(pars), 
            mean_rmse = mean(rmse), 
            max_rmse_df = max(rmse), 
            min_rmse_df = min(rmse)) %>% 
  ungroup() %>% 
  mutate(lambdas_df_name =
           paste0("lambda_0: ", lambdas_df)) %>% 
  ggplot(aes(log(mean_pars), mean_rmse)) +
  geom_hline(aes(yintercept = 0.5, linetype = " "),
             color = "red", show.legend = TRUE) +
  geom_errorbar(aes(ymin = min_rmse_df, ymax = max_rmse_df),
                position = position_dodge(.9), 
                width = 0.1) +
  geom_point(aes(fill = factor(mtry_df)), 
             size = 3, shape = 23) +
  facet_grid(gamma_df~lambdas_df, 
             scales = "free",
             switch="y", 
             labeller = 
               label_bquote(gamma==.(round(gamma_df, 3)),
                            lambda[0]==.(lambdas_df))) +
  scale_fill_manual(values = cols)+
  scale_x_continuous(
    #breaks = scales::pretty_breaks(8),
    limits = c(0, 5.53), 
    breaks = seq(0, 5.5, by = 0.5), 
    expand = c(0.02, 0.02)) +
  ylim(min_rmse, max_rmse) + 
  labs(y = expression("RMSE"[test]),
       x = expression(italic(log)~" number of selected variables"), 
       fill = expression(italic(mtry)), 
       shape = expression(gamma)) +
  theme_minimal(22) +
  theme( 
    legend.title = element_text(color = "black", size = 16),
    legend.position="bottom", 
    strip.placement = "outside", 
    axis.text.x = element_text(
      angle = 90, hjust = 1, 
      color = "black"),
    axis.title.x = element_text(margin = unit(c(10, 0, 0, 0), "mm")),
    #legend.text = element_text(angle = 90, hjust = 1, 
    #                           color = "black", size = 12), 
    axis.text.x.bottom = element_text(size = 12)) +
  scale_linetype_discrete(name = expression("RMSE"[test]~"= 0.5")) +
  guides(fill = guide_legend(override.aes = list(linetype = 0)))

pdf("plots/regression/mixture_guided_pars_rmse.pdf", height = 11.5,
    width = 12, paper='special') 
p
dev.off()

#----------------------------------------------------------------------
# Figure 4. Results using g(x_i) = Boosted_{RF}

p <- df_mix %>% 
  filter(model == "boosted_svm") %>% 
  group_by(mtry_df, lambdas_df, gamma_df) %>% 
  summarise(mean_pars = mean(pars), 
            max_pars = max(pars), 
            min_pars = min(pars), 
            mean_rmse = mean(rmse), 
            max_rmse_df = max(rmse), 
            min_rmse_df = min(rmse)) %>% 
  ungroup() %>% 
  mutate(lambdas_df_name =
           paste0("lambda_0: ", lambdas_df)) %>% 
  ggplot(aes(log(mean_pars), mean_rmse)) +
  geom_hline(aes(yintercept = 0.5, linetype = " "),
             color = "red", show.legend = TRUE) +
  geom_errorbar(aes(ymin = min_rmse_df, ymax = max_rmse_df),
                position = position_dodge(.9), 
                width = 0.1) +
  geom_point(aes(fill = factor(mtry_df)), 
             size = 3, shape = 23) +
  facet_grid(gamma_df~lambdas_df, 
             scales = "free",
             switch="y", 
             labeller = 
               label_bquote(gamma==.(round(gamma_df, 3)),
                            lambda[0]==.(lambdas_df))) +
  scale_fill_manual(values = cols)+
  scale_x_continuous(
    limits = c(0, 5.53), 
    breaks = seq(0, 5.5, by = 0.5), 
    expand = c(0.02, 0.02)) +
  ylim(min_rmse, max_rmse) + 
  labs(y = expression("RMSE"[test]),
       x = expression(italic(log)~" number of selected variables"), 
       fill = expression(italic(mtry)), 
       shape = expression(gamma)) +
  theme_minimal(22) +
  theme( 
    legend.title = element_text(color = "black", size = 16),
    legend.position="bottom", 
    strip.placement = "outside", 
    axis.text.x = element_text(
      angle = 90, hjust = 1, 
      color = "black"),
    axis.title.x = element_text(margin = unit(c(10, 0, 0, 0), "mm")),
    #legend.text = element_text(angle = 90, hjust = 1, 
    #                           color = "black", size = 12), 
    axis.text.x.bottom = element_text(size = 12)) +
  scale_linetype_discrete(name = expression("RMSE"[test]~"= 0.5")) +
  guides(fill = guide_legend(override.aes = list(linetype = 0)))

pdf("plots/regression/mixture_svm_norm_pars_rmse.pdf", height = 11.5,
    width = 12, paper='special') 
p
dev.off()

#-------------------------------------------------------------------
# End of code
#-------------------------------------------------------------------