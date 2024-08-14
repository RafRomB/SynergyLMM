
df <- read.csv("inst/extdata/Melmet5_long.csv")

head(df)
unique(df$Treatment)

df.lmm <- LMM_Model(data = df, Mouse = "Mouse", Day = "Day", Treatment = "Treatment",
                    TV = "TV", C = "Ctrl", A = "Rabusertib", B = "Cytarabine", AB = "Rabu_Cytar")

df.lmm.var <- LMM_Model(data = df, Mouse = "Mouse", Day = "Day", Treatment = "Treatment",
                    TV = "TV", C = "Ctrl", A = "Rabusertib", B = "Cytarabine", AB = "Rabu_Cytar",
                    weights = varIdent(form = ~1|Mouse))


ranef_diag(df.lmm)
ranef_diag(df.lmm.var)

resid_diag(df.lmm)
resid_diag(df.lmm.var)

ObsvsPred_plot(df.lmm, 4, 5)
ObsvsPred_plot(df.lmm.var, 4, 5)

logLik_cont(model = df.lmm, lLik_thrh = 0)
logLik_cont(model = df.lmm.var, lLik_thrh = 0.2, varName = "Mouse")

loo_logLik_disp(model = df.lmm, disp_thrh = -32.2)
loo_logLik_disp(model = df.lmm.var, disp_thrh = -32, varName = "Mouse")

Cooks_dist(df.lmm, cook_thr = 0.25)
Cooks_dist(df.lmm.var, cook_thr = 0.25)
