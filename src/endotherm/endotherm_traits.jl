"""
    EndothermThermoregulationParameters <: AbstractBehaviouralParameters

Parameters controlling endotherm thermoregulation behaviour.

# Parameters
- Thermoregulation control and limits
- Insulation (piloerection)
- Body shape adjustments
- Tissue thermal properties
- Core temperature regulation
- Panting and evaporative cooling
- Skin wetness
"""
Base.@kwdef struct EndothermThermoregulationParameters{
    TM, TO, MI,
    QM, QMR,
    IDD, IDV, IDDM, IDVM, IDDR, IDVR, IDS,
    SB, SBS, SBM,
    KF, KFS, KFM,
    TC, TCS, TCM, TCR,
    PA, PAS, PAM, PC, PM,
    SW, SWS, SWM
} <: AbstractBehaviourParameters

    # --- Control ---
    thermoregulation_mode::TM
    tolerance::TO
    max_iterations::MI

    # --- Metabolic limits ---
    Q_minimum::QM
    Q_minimum_ref::QMR

    # --- Insulation / piloerection ---
    insulation_depth_dorsal::IDD
    insulation_depth_ventral::IDV
    insulation_depth_dorsal_max::IDDM
    insulation_depth_ventral_max::IDVM
    insulation_depth_dorsal_ref::IDDR
    insulation_depth_ventral_ref::IDVR
    insulation_step::IDS

    # --- Shape change ---
    shape_b::SB
    shape_b_step::SBS
    shape_b_max::SBM

    # --- Tissue conductivity ---
    k_flesh::KF
    k_flesh_step::KFS
    k_flesh_max::KFM

    # --- Core temperature ---
    T_core::TC
    T_core_step::TCS
    T_core_max::TCM
    T_core_ref::TCR

    # --- Panting / evaporative cooling ---
    pant::PA
    pant_step::PAS
    pant_max::PAM
    pant_cost::PC
    pant_multiplier::PM

    # --- Skin wetness ---
    skin_wetness::SW
    skin_wetness_step::SWS
    skin_wetness_max::SWM
end
