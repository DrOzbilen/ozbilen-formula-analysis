
##############################################################
## DATASET GENERATION  (author's original code)
##############################################################



Height <- 150:200                                                 # 51 values
Weight <- 40:100                                                  # 61 values
hgb_female <- seq(5.0, 11.9, length.out = 70)                     # female Hb, 70 values
hgb_male <- seq(5.0, 12.9, length.out = 80)                       # male   Hb, 80 values

df_female <- expand.grid(Gender = "F", Height = Height, Weight = Weight, Hb = hgb_female)
df_male <- expand.grid(Gender = "M", Height = Height, Weight = Weight, Hb = hgb_male)

df <- rbind(df_female, df_male)

df$BMI <- df$Weight / ((df$Height / 100)^2)

df <- df[df$BMI >= 17, ]                                          # eligibility filter -> 367,650 rows
  
df$BSA <- 0.007184 * (df$Height^0.725) * (df$Weight^0.425)        # Du Bois & Du Bois equation

df$eBV <- ifelse(df$Gender == "F",                                #Nadler's equation

                   0.3561 * (df$Height/100)^3 + 0.03308 * df$Weight + 0.1833,
                   0.3669 * (df$Height/100)^3 + 0.03219 * df$Weight + 0.6041)


df$T_HGB <- ifelse(df$Gender == "F", 13, 14)                      # target Hb

df$Store <- 549.1 * df$BSA                                        # storage: dSF 95 * 5.78 * BSA


df <- df %>%
  mutate(Loss = ifelse(Gender == "F", ((Hb * 5.5 + 42) * 3), 0))  # PBAC 300, 3 cycles

df$BSA   <- round(df$BSA, 2)
df$eBV   <- round(df$eBV, 2)
df$BMI   <- round(df$BMI, 1)
df$Store <- round(df$Store)



# --- FORMULAS --------------------------------------------------------------

df$Ganzoni_RBC       <- df$Weight * (df$T_HGB - df$Hb) * 2.4
df$Ozbilen_RBC       <- (df$T_HGB - df$Hb) * df$eBV * 34.6   # 3.46 * 10
df$Ganzoni           <- df$Ganzoni_RBC + 500
df$Ozbilen_Standard  <- df$Ozbilen_RBC + df$Store

df$SDT <- ifelse(df$Hb < 10,
                   ifelse(df$Weight >= 70, 2000, 1500),
                   ifelse(df$Weight <  70, 1000, 1500))

df$Ozbilen_Advanced  <- df$Ozbilen_Standard + df$Loss



# --- differences -----------------------------------------------------------

df$Diff_OzbStan_Ganz <- df$Ozbilen_Standard - df$Ganzoni
df$Diff_OzbAdv_Ganz  <- df$Ozbilen_Advanced - df$Ganzoni
df$Diff_OzbStan_SDT  <- df$Ozbilen_Standard - df$SDT
df$Diff_OzbAdv_SDT   <- df$Ozbilen_Advanced - df$SDT




# round dose columns up (ceiling)

dose_cols <- c("Ganzoni_RBC","Ganzoni","Ozbilen_RBC","Ozbilen_Standard",
               "Ozbilen_Advanced","Store","Loss","SDT")

df <- df %>% mutate(across(all_of(dose_cols), ceiling))


#################################################################################
