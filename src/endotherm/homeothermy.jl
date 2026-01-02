function endotherm_thermoregulation_original(
    Q_gen,
    T_skin,
    T_insulation,
    organism,
    thermoregulation_pars,
    environment,
    model_pars,
)
    endotherm_out = nothing # initialise

    (; thermoregulation_mode, tolerance, max_iterations,
        Q_minimum, Q_minimum_ref,
        insulation_depth_dorsal, insulation_depth_ventral,
        insulation_depth_dorsal_max, insulation_depth_ventral_max,
        insulation_depth_dorsal_ref, insulation_depth_ventral_ref, insulation_step,
        shape_b, shape_b_step, shape_b_max,
        k_flesh, k_flesh_step, k_flesh_max,
        T_core, T_core_step, T_core_max, T_core_ref,
        pant, pant_step, pant_max, pant_cost, pant_multiplier,
        skin_wetness, skin_wetness_step, skin_wetness_max) = thermoregulation_pars

    #check if starting in piloerect state
    if insulation_step > 0.0 && (insulation_depth_dorsal + insulation_depth_ventral) > 0u"mm"
        # start with erect insulation
        insulation_depth_dorsal, insulation_depth_ventral, organism =
            piloerect!(
                insulation_depth_dorsal_max,
                insulation_depth_ventral_max,
                insulation_depth_dorsal_ref,
                insulation_depth_ventral_ref,
                0, # don't step it down
                organism)
    end

    endotherm_out = solve_metabolic_rate(T_skin, T_insulation, organism, environment, model_pars)
    T_skin = endotherm_out.thermoregulation.T_skin
    T_insulation = endotherm_out.thermoregulation.T_insulation
    Q_gen = endotherm_out.energy_fluxes.Q_gen
        
    iteration = 0

    # start of thermoregulation loop
    while Q_gen < Q_minimum * (1 - tolerance)

        iteration =+ 1
        if iteration > max_iterations
            @warn "max_iterations exceeded"
            return
        end
        # -----------------------------------------------------------------------------
        # 1. Reduce insulation (piloerection)
        # -----------------------------------------------------------------------------        
        if (insulation_depth_dorsal > insulation_depth_dorsal_ref) &&
           (insulation_depth_ventral > insulation_depth_ventral_ref)

            insulation_depth_dorsal, insulation_depth_ventral, organism =
                piloerect!(
                    insulation_depth_dorsal,
                    insulation_depth_ventral,
                    insulation_depth_dorsal_ref,
                    insulation_depth_ventral_ref,
                    insulation_step,
                    organism)

        # -------------------------------------------------------------------------------
        # 2. Uncurl (increase surface area)
        # -------------------------------------------------------------------------------                   
        elseif shape_b < shape_b_max
            shape_b, organism = uncurl!(shape_b, shape_b_step, shape_b_max, organism)

        # -------------------------------------------------------------------------------
        # 3. Vasodilate (increase k_flesh)
        # -------------------------------------------------------------------------------            
        elseif k_flesh < k_flesh_max
            k_flesh, organism = vasodilate!(k_flesh, k_flesh_step, k_flesh_max, organism)

        # -------------------------------------------------------------------------------
        # 4. Allow core temperature to rise (and possibly pant and sweat in parallel)
        # -------------------------------------------------------------------------------               
        elseif T_core < T_core_max
            T_core, Q_minimum, organism = hyperthermia!(T_core, T_core_step, T_core_max,
                T_core_ref, Q_minimum_ref, pant_cost, organism)
            if thermoregulation_mode >= 2 && pant < pant_max
                # pant in parallel to allowing core temperature to rise
                pant, pant_cost, Q_minimum, organism = pant!(pant, pant_step, pant_max,
                    T_core_ref, Q_minimum_ref, pant_multiplier, organism)
            end
            if thermoregulation_mode == 3
                # sweat in parallel to allowing core temperature to rise and panting
                if (skin_wetness > skin_wetness_max) || (skin_wetness_step <= 0)
                    @warn "All thermoregulatory options exhausted"
                    return
                end
                skin_wetness, organism = sweat!(skin_wetness, skin_wetness_step,
                    skin_wetness_max, organism)
            end

        # -------------------------------------------------------------------------------
        # 5. Pant to dump heat evaporatively (and possibly sweat in parallel)
        # ------------------------------------------------------------------------------- 
        elseif pant < pant_max
            pant, pant_cost, Q_minimum, organism = pant!(pant, pant_step, pant_max,
                T_core_ref, Q_minimum_ref, pant_multiplier, organism)
            if thermoregulation_mode == 3
                if (skin_wetness > skin_wetness_max) || (skin_wetness_step <= 0)
                    @warn "All thermoregulatory options exhausted"
                    return
                end
                # sweat in parallel to panting
                skin_wetness, organism = sweat!(skin_wetness, skin_wetness_step,
                    skin_wetness_max, organism)
            end

        # -------------------------------------------------------------------------------
        # 6. Sweat to dump heat evaporatively
        # -------------------------------------------------------------------------------            
        else
            if (skin_wetness > skin_wetness_max) || (skin_wetness_step <= 0)
                return
            end
            skin_wetness, organism = sweat!(skin_wetness, skin_wetness_step,
                skin_wetness_max, organism)
        end

        endotherm_out = solve_metabolic_rate(T_skin, T_insulation, organism, environment, model_pars)
        T_skin = endotherm_out.thermoregulation.T_skin
        T_insulation = endotherm_out.thermoregulation.T_insulation
        Q_gen = endotherm_out.energy_fluxes.Q_gen

    end
    return endotherm_out
end