library(NicheMapR)

Ww_g <- 500
Usrhyt <- 0.05
alpha <- 0.85
T_F_min <- 33
T_F_max <- 43
T_B_min <- 18
CT_max <- 48
shape_b <- 1/5
shape_c <- 1/5
rho_body <- 1000
c_body <- 3073
q <- 0
k_flesh <- 0.5
geom <- 2

loc <- c(130, -25)
maxshade <- 90
micro <- micro_global(loc = loc, maxshade = maxshade, Usrhyt = Usrhyt)
metout <- as.data.frame(micro$metout)
soil <- as.data.frame(micro$soil)
shadmet <- as.data.frame(micro$shadmet)
shadsoil <- as.data.frame(micro$shadsoil)
elevation <- micro$elev
press <- 101325 * ((1 - (0.0065 * elevation / 288)) ^ (1 / 0.190284))

DOYs <- unique(metout$DOY)
simday <- DOYs[1]  # January (DOY 15), hottest month at this southern-hemisphere site

metout_in <- subset(metout, DOY == simday)
shadmet_in <- subset(shadmet, DOY == simday)
soil_in <- subset(soil, DOY == simday)
shadsoil_in <- subset(shadsoil, DOY == simday)

trans <- trans_behav(Ww_g = Ww_g, alpha = alpha, T_F_min = T_F_min, T_F_max = T_F_max,
                      CT_max = CT_max, T_B_min = T_B_min, geom = geom, shape_b = shape_b, shape_c = shape_c,
                      rho_body = rho_body, k_flesh = k_flesh, q = q, lump = 1,
                      metout = metout_in, shadmet = shadmet_in, soil = soil_in, shadsoil = shadsoil_in,
                      press = press, alpha_sub = 1 - micro$REFL, shade = micro$maxshade[1])

day_results <- as.data.frame(trans$day_results)
sum_stats <- as.data.frame(trans$sum_stats)
act_window <- as.data.frame(trans$act_window)

outdir <- "../data/trans_behav"
write.csv(metout_in, file.path(outdir, "trans_behav_metout.csv"), row.names = FALSE)
write.csv(shadmet_in, file.path(outdir, "trans_behav_shadmet.csv"), row.names = FALSE)
write.csv(soil_in, file.path(outdir, "trans_behav_soil.csv"), row.names = FALSE)
write.csv(shadsoil_in, file.path(outdir, "trans_behav_shadsoil.csv"), row.names = FALSE)
write.csv(day_results, file.path(outdir, "trans_behav_day_results.csv"), row.names = FALSE)
write.csv(sum_stats, file.path(outdir, "trans_behav_sum_stats.csv"), row.names = FALSE)
write.csv(act_window, file.path(outdir, "trans_behav_act_window.csv"), row.names = FALSE)

params <- data.frame(
  Ww_g=Ww_g, alpha=alpha, T_F_min=T_F_min, T_F_max=T_F_max, T_B_min=T_B_min, CT_max=CT_max,
  shape_b=shape_b, shape_c=shape_c, rho_body=rho_body, c_body=c_body, q=q, k_flesh=k_flesh, geom=geom,
  press=press, alpha_sub=1-micro$REFL, shade=micro$maxshade[1], DOY=simday, loc_lon=loc[1], loc_lat=loc[2],
  Usrhyt=Usrhyt
)
write.csv(params, file.path(outdir, "trans_behav_params.csv"), row.names = FALSE)
