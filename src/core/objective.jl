################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

"expected cost of active power generation"
function objective_min_expected_generation_cost(pm::AbstractPowerModel; kwargs...)
    gen_cost = Dict()

    T2 = pm.data["T2"]

    for (g, gen) in _PM.ref(pm, :gen, nw=1)
        pg = Dict(nw => _PM.var(pm, nw, :pg, g) for nw in _PM.nw_ids(pm))

        if length(gen["cost"]) == 1
            gen_cost[g] = gen["cost"][1]
        elseif length(gen["cost"]) == 2
            gen_cost[g] = gen["cost"][1]*pg[1] + 
                          gen["cost"][2]
        elseif length(gen["cost"]) == 3
            gen_cost[g] = gen["cost"][1]*sum(T2.get([n-1,n-1]) * pg[n]^2 for n in _PM.nw_ids(pm)) + 
                          gen["cost"][2]*pg[1] + 
                          gen["cost"][3]
        else
            gen_cost[g] = 0.0
        end
    end

    return JuMP.@objective(pm.model, Min,
            sum(gen_cost[g] for g in _PM.ids(pm, :gen, nw=1))
    )
end

"expected max PV generation"
function objective_max_PV(pm::AbstractPowerModel; kwargs...)
    p_size = Dict()


    for (p, PV) in _PM.ref(pm, :PV, nw=1)
        p_size[p] = Dict(nw => _PM.var(pm, nw, :p_size, p) for nw in [1])
    end

    return JuMP.@objective(pm.model, Max,
            sum(p_size[p][1] for p in _PM.ids(pm, :PV, nw=1))
    )
end


"expected max PV generation"
function objective_max_PV_det(pm::AbstractPowerModel; kwargs...)
    p_size = Dict()


    for (p, PV) in _PM.ref(pm, :PV)
        p_size[p] = _PM.var(pm, :p_size,p)
    end

    return JuMP.@objective(pm.model, Max,
            sum(p_size[p] for p in _PM.ids(pm, :PV))
    )
end

"expected min PV generation"
function objective_min_PV_det(pm::AbstractPowerModel; kwargs...)
    p_size = Dict()


    for (p, PV) in _PM.ref(pm, :PV)
        p_size[p] = _PM.var(pm, :p_size,p)
    end

    return JuMP.@objective(pm.model, Min,
            sum(p_size[p] for p in _PM.ids(pm, :PV))
    )
end

"objective_max_PV_equal_for_all_consumer"
function objective_max_PV_equal_for_all_consumer(pm::AbstractPowerModel; kwargs...)
    p_size = Dict()


    for (p, PV) in _PM.ref(pm, :PV, nw=1)
        p_size[p] = Dict(nw => _PM.var(pm, nw, :p_size, p) for nw in [1])
    end

    return JuMP.@objective(pm.model, Max,
            p_size[1][1]
    )
end

"objective_max_PV_equal_for_all_consumer"
function objective_max_PV_det_equal_for_all_consumer(pm::AbstractPowerModel; kwargs...)
    p_size = Dict()


    for (p, PV) in _PM.ref(pm, :PV)
        p_size[p] = _PM.var(pm, :p_size,p)
    end

    # return JuMP.@objective(pm.model, Max,
    #         sum(p_size[p] for p in _PM.ids(pm, :PV))
    # )
    if length(p_size)==0
        return JuMP.@objective(pm.model, Max,
        0)
    else
        return JuMP.@objective(pm.model, Max,
        p_size[1])
    end
    

end