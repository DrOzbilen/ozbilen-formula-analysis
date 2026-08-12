###############################################################################
## DETERMINISTIC ANALYSIS + TABLES
###############################################################################



formulas <- c("Ganzoni","Ozbilen_Standard","SDT","Ozbilen_Advanced")
pretty_names <- c(Ganzoni="Ganzoni", SDT="SDT",
                  Ozbilen_Standard="Standard Version",
                  Ozbilen_Advanced="Advanced Version")



## --- TABLE 2: Descriptive statistics of dataset-----------------------------

table2 <- df %>%
  group_by(Gender) %>% 
  summarise(
    
    across(.cols = c(Height, Weight, Hb, BSA, BMI, eBV), 
           .fns = list(
             Mean = ~mean(., na.rm = TRUE),
             Median = ~median(., na.rm = TRUE),
             Min = ~min(., na.rm = TRUE),
             Max = ~max(., na.rm = TRUE),
             SD = ~sd(., na.rm = TRUE)
           ),
           .names = "{.col}_{.fn}" 
    ),
    N = n()
  )

print(table2)

write.csv(table2, "Table2.csv")



## --- TABLE 3: dose by formula and sex --------------------------------------

table3 <- df %>%
  pivot_longer(all_of(formulas), names_to="Formula", values_to="Dose") %>%
  group_by(Formula, Gender) %>%
  summarise(summary = med_iqr_range(Dose), .groups="drop") %>%
  pivot_wider(names_from=Gender, values_from=summary) %>%
  mutate(Formula = pretty_names[Formula]) %>%
  arrange(match(Formula, c("Ganzoni","SDT","Standard Version","Advanced Version")))

print(table3)

write.csv(table3, "Table3.csv")




## --- TABLE 4: agreement (mean diff + 95% LoA; descriptive) -----------------

agree_one <- function(a,b,sub){ d<-sub[[a]]-sub[[b]]; m<-mean(d); s<-sd(d)
c(mean_diff=m, lo=m-1.96*s, hi=m+1.96*s) }
pairs_ab <- list(c("Ozbilen_Standard","Ganzoni","Standard vs Ganzoni"),
                 c("Ozbilen_Advanced","Ganzoni","Advanced vs Ganzoni"),
                 c("Ozbilen_Standard","SDT","Standard vs SDT"),
                 c("Ozbilen_Advanced","SDT","Advanced vs SDT"))
table4 <- do.call(rbind, lapply(pairs_ab, function(p)
  do.call(rbind, lapply(c("F","M"), function(g){
    sub<-df[df$Gender==g,]; r<-agree_one(p[1],p[2],sub)
    data.frame(Comparison=p[3], Sex=g, Mean_diff=round(r["mean_diff"],1),
               LoA_low=round(r["lo"],0), LoA_high=round(r["hi"],0), row.names=NULL)
  }))))

print(table4)

write.csv(table4, "Table4.csv")



## --- TABLE 5: Ganzoni 2.4 coefficient decomposition ------------------------

df$Coeff <- df$eBV * 34.6 / df$Weight
table5 <- df %>% group_by(Gender) %>%
  summarise(Median=round(median(Coeff),2),
            IQR=sprintf("%.2f-%.2f", q1(Coeff), q3(Coeff)),
            Range=sprintf("%.2f-%.2f", min(Coeff), max(Coeff)),
            `2.4 higher by (%)`=round((2.4-median(Coeff))/median(Coeff)*100,1),
            .groups="drop")
table5 <- rbind(table5, data.frame(Gender="All",
                                   Median=round(median(df$Coeff),2),
                                   IQR=sprintf("%.2f-%.2f", q1(df$Coeff), q3(df$Coeff)),
                                   Range=sprintf("%.2f-%.2f", min(df$Coeff), max(df$Coeff)),
                                   `2.4 higher by (%)`=round((2.4-median(df$Coeff))/median(df$Coeff)*100,1),
                                   check.names=FALSE))

print(table5)

write.csv(table5, "Table5.csv")

###############################################################################
