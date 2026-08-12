##### FIGURES ######################################



## --- clinical thresholds (Figure 2 data) -----------------------------------


thr_pairs <- list(c("Ozbilen_Standard","Ganzoni","Standard vs Ganzoni"),
                  c("Ozbilen_Advanced","Ganzoni","Advanced vs Ganzoni"),
                  c("Ozbilen_Standard","SDT","Standard vs SDT"),
                  c("Ozbilen_Advanced","SDT","Advanced vs SDT"))
threshold_tab <- do.call(rbind, lapply(thr_pairs, function(p)
  do.call(rbind, lapply(c("F","M"), function(g){
    sub<-df[df$Gender==g,]; d<-abs(sub[[p[1]]]-sub[[p[2]]])
    data.frame(Comparison=p[3], Sex=g,
               pct_ge_250=round(mean(d>=250)*100,1),
               pct_ge_500=round(mean(d>=500)*100,1),
               pct_ge_1000=round(mean(d>=1000)*100,1), row.names=NULL)
  }))))

write.csv(threshold_tab, "Table_thresholds.csv"), row.names=FALSE))

print(threshold_tab)





## --- ferritin target sensitivity -------------------------------------------


ferr_sens <- do.call(rbind, lapply(c(50,70,100), function(tgt){
  store <- (tgt-5) * 5.78 * df$BSA
  std <- df$Ozbilen_RBC + store; adv <- std + df$Loss
  data.frame(Target=tgt, Sex=c("F","M"),
             Std_median=c(median(std[df$Gender=="F"]), median(std[df$Gender=="M"])),
             Adv_median=c(median(adv[df$Gender=="F"]), median(adv[df$Gender=="M"])),
             Store_median=round(median(store)))
}))

write.csv(ferr_sens,"Table_ferritin_sensitivity.csv", row.names=FALSE))
print(ferr_sens)




##############################################################################
## --- FIGURE 1: clinical thresholds -----------------------------------------


fig1_data <- threshold_tab %>%
  pivot_longer(c(pct_ge_250,pct_ge_500), names_to="threshold", values_to="pct") %>%
  mutate(threshold=recode(threshold, pct_ge_250="Difference >= 250 mg",
                          pct_ge_500="Difference >= 500 mg"),
         Sex=recode(Sex, F="Female", M="Male"),
         Comparison=factor(Comparison,
                           levels=c("Standard vs Ganzoni","Advanced vs Ganzoni",
                                    "Standard vs SDT","Advanced vs SDT")))
fig1 <- ggplot(fig1_data, aes(pct, Comparison, fill=Sex)) +
  geom_col(position=position_dodge(.7), width=.6) +
  geom_text(aes(label=sprintf("%.0f",pct)), position=position_dodge(.7),
            hjust=-0.2, size=3) +
  facet_wrap(~threshold) +
  scale_fill_manual(values=c(Female=C_F, Male=C_M)) +
  scale_x_continuous(limits=c(0,110), name="% of body profiles") +
  labs(y=NULL, title="Proportion of profiles with clinically meaningful dose differences") +
  theme_pub

ggsave(file.path("Figure1_clinical_thresholds.png"), fig1, width=9.5, height=4.2, dpi=300)





##############################################################################
## --- FIGURE 2: deviation from four-formula mean ----------------------------

df$row_mean <- rowMeans(df[, formulas])
dev_long <- df %>% select(Gender, all_of(formulas), row_mean) %>%
  pivot_longer(all_of(formulas), names_to="Formula", values_to="Dose") %>%
  mutate(Deviation=Dose-row_mean,
         Formula=factor(pretty_names[Formula],
                        levels=c("Ganzoni","SDT","Standard Version","Advanced Version")),
         Sex=recode(Gender, F="Female", M="Male"))
fig2 <- ggplot(dev_long, aes(Formula, Deviation, fill=Formula)) +
  geom_hline(yintercept=0, colour="grey20") +
  geom_boxplot(outlier.shape=NA, width=.6, alpha=.85) +
  facet_wrap(~Sex) +
  scale_fill_manual(values=c("Ganzoni"=C_GANZ,"SDT"=C_SDT,
                             "Standard Version"=C_STD,"Advanced Version"=C_ADV)) +
  coord_cartesian(ylim=c(-700,650)) +
  labs(x=NULL, y="Deviation from group mean (mg)",
       title="Systematic dose deviation from the four-formula average",
       caption="Above the line: higher dose estimates | Below the line: lower estimates") +
  theme_pub + theme(legend.position="none",
                    axis.text.x=element_text(angle=20, hjust=1))

ggsave(file.path("Figure2_deviation.png"), fig2, width=9.5, height=4.6, dpi=300)




##############################################################################
## --- FIGURE 3: uncovered menstrual iron loss -------------------------------

fig3_data <- expand.grid(PBAC=c(200,300,500,800,1000), Hb=c(7,9,11)) %>%
  mutate(uncovered=((Hb*5.5)+42)*(PBAC/100),
         Hb_label=factor(paste0("Hb ",Hb," g/dL")))
fig3 <- ggplot(fig3_data, aes(PBAC, uncovered, colour=Hb_label)) +
  geom_line(linewidth=1.1) + geom_point(size=2.5) +
  scale_colour_manual(values=c("Hb 7 g/dL"="#8B0000","Hb 9 g/dL"="#C44E52",
                               "Hb 11 g/dL"="#E8A598"), name="Baseline haemoglobin") +
  labs(x="PBAC score", y="Menstrual iron not covered over 3 cycles (mg)",
       title="Iron deficit left uncovered by menstruation-blind formulas") +
  theme_pub

ggsave(file.path("Figure3_uncovered_menstrual.png"), fig3, width=8, height=5, dpi=300)




##############################################################################
## --- FIGURE 4: iron budget (single-patient scenario) -----------------------

h_m<-1.60; wt<-60
eBV_pt<-0.3561*h_m^3+0.03308*wt+0.1833
BSA_pt<-0.007184*(160^0.725)*(60^0.425)
cost_ferr<-5.78*BSA_pt; iron_per_dHb<-eBV_pt*3.46*10
DOSE<-1000; base_hb<-10; target_hb<-13; base_ferr<-10
iron_hb_needed<-(target_hb-base_hb)*iron_per_dHb
ferr_ceiling<-base_ferr+(DOSE-iron_hb_needed)/cost_ferr
mens_loss<-function(hb,pbac) ((hb*5.5)+42)*(pbac/100)
scen<-do.call(rbind, lapply(c(200,300,500,800,1000), function(pbac){
  ml<-mens_loss(base_hb,pbac); after<-DOSE-ml
  if(after>=iron_hb_needed){ hb_mg<-iron_hb_needed; hb_fin<-target_hb
  ferr_mg<-after-iron_hb_needed; ferr_fin<-base_ferr+ferr_mg/cost_ferr
  } else { hb_mg<-max(after,0); hb_fin<-base_hb+hb_mg/iron_per_dHb
  ferr_mg<-0; ferr_fin<-base_ferr }
  data.frame(PBAC=pbac, mens_mg=min(ml,DOSE), hb_mg=hb_mg, ferr_mg=ferr_mg,
             hb_ach=(hb_fin-base_hb)/(target_hb-base_hb)*100,
             ferr_ach=max((ferr_fin-base_ferr)/(ferr_ceiling-base_ferr)*100,0))
}))
write.csv(scen, file.path("Table_iron_budget_scenario.csv"), row.names=FALSE)

ach <- scen %>% select(PBAC, Haemoglobin=hb_ach, Ferritin=ferr_ach) %>%
  pivot_longer(c(Haemoglobin,Ferritin), names_to="metric", values_to="pct")
alloc <- scen %>% transmute(PBAC,
                            `Menstrual loss`=-mens_mg/DOSE*100,
                            `Haemoglobin`=-hb_mg/DOSE*100,
                            `Iron stores`=-ferr_mg/DOSE*100) %>%
  pivot_longer(c(`Menstrual loss`,`Haemoglobin`,`Iron stores`),
               names_to="component", values_to="pct") %>%
  mutate(component=factor(component, levels=c("Menstrual loss","Haemoglobin","Iron stores")))
fig4 <- ggplot() +
  geom_col(data=ach, aes(factor(PBAC), pct, fill=metric),
           position=position_dodge(.7), width=.6) +
  geom_col(data=alloc, aes(factor(PBAC), pct, fill=component), width=.72) +
  geom_hline(yintercept=0, colour="grey20", linewidth=.7) +
  geom_hline(yintercept=100, colour="#2a9d5c", linetype="dashed") +
  scale_fill_manual(values=c("Haemoglobin"="#4C72B0","Ferritin"="#C44E52",
                             "Menstrual loss"="#8C6D31","Iron stores"="#E0A0A2"), name=NULL) +
  scale_y_continuous(breaks=c(-100,-50,0,50,100),
                     labels=c("100%","50%","0","50%","100%"),
                     name="% dose allocated  |  % target achieved") +
  labs(x="PBAC score", title="The iron budget: dose allocation and target attainment",
       caption=paste0("Illustrative single-patient scenario (Hb 10, ferritin 10, 160 cm, 60 kg; ",
                      "attainable ferritin ", round(ferr_ceiling), " ug/L). Requires clinical validation.")) +
  theme_pub

ggsave(file.path("Figure4_iron_budget.png"), fig4, width=10, height=6.5, dpi=300)





###############################################################################
## SUPPLEMENTARY FIGURE S1 - HEATMAPS
###############################################################################


## ---------------------------------------------------------------------------
## Base heatmap function (0.5 g/dL Hb bins, 5 kg weight bins)
## ---------------------------------------------------------------------------


heatmap_plot <- function(data, gender, variable, title, save_path = NULL) {
  data$hb_bin     <- cut(data$Hb,     breaks = seq(5, 12, by = 0.5), include.lowest = TRUE)
  data$weight_bin <- cut(data$Weight, breaks = seq(50, 100, by = 5), include.lowest = TRUE)
  data$hb_bin     <- factor(data$hb_bin, levels = rev(levels(data$hb_bin)))
  
  filtered_data <- data %>% filter(Gender == gender)
  
  heatmap_data <- filtered_data %>%
    group_by(hb_bin, weight_bin) %>%
    summarise(avg_diff = mean(.data[[variable]], na.rm = TRUE), .groups = "drop") %>%
    na.omit()
  
  p <- ggplot(heatmap_data, aes(x = weight_bin, y = hb_bin, fill = avg_diff)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    theme_minimal() +
    labs(title = title,
         x = "Weight (kg)",
         y = "Haemoglobin (g/dL)",
         fill = "mg") +
    geom_text(aes(label = round(avg_diff, 1)), color = "black", size = 3)
  
  if (!is.null(save_path)) {
    if (!grepl("\\.png$", save_path)) save_path <- paste0(save_path, ".png")
    ggsave(save_path, plot = p, width = 6, height = 6, dpi = 600, bg = "white")
  }
  return(p)
}


################################################################################
## Combined heatmaps: 5 variables x 2 sexes (2 columns, 5 rows)

variables <- c("Diff_RBC", "Diff_OzbStan_Ganz", "Diff_OzbAdv_Ganz",
               "Diff_OzbStan_SDT", "Diff_OzbAdv_SDT")

plots <- list()

for (variable in variables) {
  for (gender in c("F", "M")) {
    p <- heatmap_plot(
      data = df,
      gender = gender,
      variable = variable,
      title = paste0(plot_labels[[variable]], " - ", ifelse(gender == "F", "Female", "Male")),
      save_path = NULL
    )
    plots[[paste(variable, gender, sep = "_")]] <- p
  }
}

final_plot <- (
  plots[["Diff_RBC_F"]]          + plots[["Diff_RBC_M"]] +
    plots[["Diff_OzbStan_Ganz_F"]] + plots[["Diff_OzbStan_Ganz_M"]] +
    plots[["Diff_OzbAdv_Ganz_F"]]  + plots[["Diff_OzbAdv_Ganz_M"]] +
    plots[["Diff_OzbStan_SDT_F"]]  + plots[["Diff_OzbStan_SDT_M"]] +
    plots[["Diff_OzbAdv_SDT_F"]]   + plots[["Diff_OzbAdv_SDT_M"]]
) + plot_layout(ncol = 2)

ggsave("Combined_Heatmaps_5Vars.png", plot = final_plot,
       width = 16, height = 30, dpi = 600, bg = "white")

###############################################################################