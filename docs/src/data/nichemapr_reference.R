# Reference output of NicheMapR for the tutorials "A mammal across air temperatures" and
# "A human that thermoregulates".
#
# Run once with NicheMapR installed; the docs read the CSV file this writes.
#   Rscript nichemapr_reference.R

library(NicheMapR)

# endoR_devel: the default 65 kg animal with thermoregulation, limits as in
# `example_thermoregulation_limits()`.
# DELTAR = 100 makes the exhaled air leave at lung temperature, as HeatExchange.jl does.

TAs <- seq(0, 50, 2)
endo <- function(TA) {
  out <- endoR_devel(TA = TA, THERMOREG = 1, TREGMODE = 1, DELTAR = 100,
                     TC = 37, TC_MAX = 39, TC_INC = 0.1,
                     AK1 = 0.9, AK1_MAX = 2.8, AK1_INC = 0.1,
                     SHAPE_B = 1.1, SHAPE_B_MAX = 5, UNCURL = 0.1,
                     PANT = 1, PANT_MAX = 10, PANT_INC = 0.1, PANT_MULT = 1.05,
                     PCTWET = 0.5, PCTWET_MAX = 100, PCTWET_INC = 0.1)
  data.frame(air_temperature_C = TA,
             metabolic_W = out$enbal[, "QGEN"], evaporation_W = out$enbal[, "QEVAP"],
             core_C = out$treg[, "TC"], skin_C = out$treg[, "TSKIN_D"], fur_C = out$treg[, "TFA_D"],
             lung_C = out$treg[, "TLUNG"], shape_b = out$treg[, "SHAPE_B"],
             flesh_conductivity = out$treg[, "K_FLESH"], pant = out$treg[, "PANT"],
             skin_wetness_pct = out$treg[, "PCTWET"],
             respiratory_water_g_h = out$masbal[, "H2OResp_g"], cutaneous_water_g_h = out$masbal[, "H2OCut_g"])
}
write.csv(do.call(rbind, lapply(TAs, endo)), "endoR_thermoreg_reference.csv", row.names = FALSE)

# HomoTherm: the default human from the cold to the heat, with its thermoregulation.

parts <- c("head", "trunk", "arm", "leg")
homo <- function(TA) {
  out <- HomoTherm(TA = TA, VEL = 0.1, RH = 50)
  part_rows <- do.call(rbind, lapply(parts, function(p) {
    treg <- out[[paste0(p, ".treg")]]
    enbal <- out[[paste0(p, ".enbal")]]
    data.frame(air_temperature_C = TA, part = p,
               core_C = treg["T_CORE"], skin_dorsal_C = treg["TSKIN_D"], skin_ventral_C = treg["TSKIN_V"],
               flesh_conductivity = treg["K_FLESH"], skin_wetness_pct = treg["PCTWET"],
               heat_generated_W = enbal["QMETAB"], evaporation_W = enbal["QEVAP"])
  }))
  b <- out$balance
  whole <- data.frame(air_temperature_C = TA, metabolic_W = b["QMETAB"], core_C = b["T_CORE"], lung_C = b["T_LUNG"],
                      skin_C = b["T_SKIN"], surface_C = b["T_CLO"], flesh_conductivity = b["K_FLESH"],
                      skin_wetness_pct = b["PCTWET"], cutaneous_water_L_h = b["EVAP_CUT_L"],
                      respiratory_water_L_h = b["EVAP_RESP_L"], sweat_L_h = b["SWEAT_L"])
  list(parts = part_rows, whole = whole)
}
homo_out <- lapply(seq(-10, 46, 2), homo)
write.csv(do.call(rbind, lapply(homo_out, `[[`, "whole")), "homotherm_whole.csv", row.names = FALSE)
