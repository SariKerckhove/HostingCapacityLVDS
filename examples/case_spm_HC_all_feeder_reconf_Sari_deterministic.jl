"""
################################################################################
#  A. Koirala, T. V. Aacker, M. U. Hashmi, R. D’hulst, and D. V. Hertem,       #
“Chance-constrained optimization based PV hosting capacity calculation using   #
general Polynomial Chaos,”                                                     #
Submitted to IEE Transactions of Power Systems.                                #
[Online]: https://tinyurl.com/yckbp58j                                         #
################################################################################
# Based on StochasticPowerModels.jl                                            #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################
# This example is for Case study on Section V of the paper      #
# The plot 1 shows the relationship between computation time and sources of uncertainties 
# in the feeder. Which is Fig. 11 in the paper above.
################################################################################

# variables
"""

using Pkg
Pkg.activate(".")
using JuMP
using Ipopt
using PowerModels
using StochasticPowerModels
using PowerModelsDistribution
using JSON
using DataFrames
using CSV
using Statistics
using Plots
# constants 
const PM = PowerModels
const SPM = StochasticPowerModels

# solvers
# ipopt_solver = Ipopt.Optimizer

"Converts JSON file of three phase DN to single phase equivalent"
function build_mathematical_model_reconfiguration(dir, config_file_name, load_dist_csv, pv_dist_csv; t_s=52, pd=0.0, qd=0.0, scale_factor=1.0, curt=0.0, cross_area_fact=1.0)
    #configuration = "star"

    # build_mathematical_model_single_phase(file, feeder,load_file, pv_file, t_s= 59)
    # file  = joinpath(BASE_DIR, "test/data/Spanish/") --> dir
    # feeder="All_feeder/"*b.conf --> config file name
    # load_file= "beta_lm_2016_8_6.csv" --> load_dist_csv
    # pv_file = "beta_pm_2016_8_6.csv" --> pv_dist_cv

    """
    Specify voltage and power base, voltage base should be the phase-to-ground voltage
    of the feeder studied (kV), the power base can be arbitrairly chosen (MW)
    """
    voltage_base = 0.230  # (kV)
    power_base = 0.5  # (MW)
    Z_base = voltage_base^2 / power_base # (Ohm)
    current_base = power_base / (voltage_base * 1e-3) # (A)

    mwpu = 1 / power_base
    kwpu = (1e-3) / power_base

    network_model = Dict{String,Any}()
    configuration_json_dict = Dict{Any,Any}()
    # device_df=CSV.read(dir*config_file_name[1:length(config_file_name)-19]*".csv", DataFrame) #csv file with categories.

    # dist_lv = CSV.read(dir * load_dist_csv, DataFrame) # load_dist_csv
    dist_pv = CSV.read(dir * pv_dist_csv, DataFrame)
    # dist_pv=CSV.read(dir*"beta_pm_2022_181"*".csv", DataFrame)
    dist_pv_ts = dist_pv[in([t_s]).(dist_pv.timeslot), :]
    # dist_lv_ts = dist_lv[in([t_s]).(dist_lv.timeslot), :]
    # # requires the csv file with categories.

    # dist_lv_ts_feeder = dist_lv_ts[in(unique(device_df.category)).(dist_lv_ts.cluster), :]

    s_dict = Dict()

    i = 1
    # for dist in eachrow(dist_lv_ts_feeder)
    #     s = Dict()
    #     s["dst"] = "Beta"
    #     s["dst_id"] = dist["cluster"]
    #     s["pa"] = dist["alpha"]
    #     s["pb"] = dist["beta"]
    #     s["pc"] = dist["lower"]
    #     s["pd"] = dist["lower"] + dist["upper"]
    #     s_dict[string(i)] = s
    #     i = i + 1
    # end


    ##add Irradiance if day time or there is some Irradiance
    if dist_pv_ts.upper[1] > 0
        s = Dict()
        s["dst"] = "Beta"
        s["dst_id"] = 55
        s["pa"] = dist_pv_ts[!, "alpha"][1]
        s["pb"] = dist_pv_ts[!, "beta"][1]
        s["pc"] = dist_pv_ts[!, "lower"][1]
        s["pd"] = dist_pv_ts[!, "lower"][1] + dist_pv_ts[!, "upper"][1]
        s_dict[string(i)] = s
    end


    #network_model["is_kron_reduced"] = true
    network_model["p_factor"] = 0.95
    network_model["q_factor"] = sqrt(1 - 0.95^2)
    network_model["dcline"] = Dict{String,Any}()
    network_model["switch"] = Dict{String,Any}()
    #network_model["is_projected"] = true
    network_model["per_unit"] = true
    #network_model["data_model"] = MATHEMATICAL
    network_model["shunt"] = Dict{String,Any}()
    network_model["transformer"] = Dict{String,Any}()
    network_model["bus"] = Dict{String,Any}()
    network_model["map"] = Dict{String,Any}()
    #network_model["conductors"] = 1
    network_model["baseMVA"] = power_base
    network_model["basekv"] = voltage_base
    network_model["bus_lookup"] = Dict{Any,Int64}()
    network_model["run_type"] = 1
    network_model["load"] = Dict{String,Any}()
    network_model["gen"] = Dict{String,Any}("1" => Dict{String,Any}(
        "pg" => 0.2,
        "model" => 2,
        #"connections"   => [1, 2, 3],
        "shutdown" => 0.0,
        "startup" => 0.0,
        #"configuration" => WYE,
        "name" => "virtual_generator",
        "qg" => 0.0,
        "gen_bus" => 1, # always starting from 1 for slackbus
        "vbase" => voltage_base,
        "source_id" => Any["gen", 1],
        "index" => 1,
        "cost" => [20000.0, 1400.0, 0.0], #?? (@Arpan, what is this cost?)
        "gen_status" => 1,
        "qmax" => 1.275,
        "qmin" => -1.275,
        "pmax" => 1.5,
        "pmin" => -1.5,
        "ncost" => 3,
        "λpmin" => 1.65, #1.03643 ,
        "λpmax" => 1.65, #1.03643 ,
        "λqmin" => 1.65, #1.03643 ,
        "λqmax" => 1.65 #1.03643
    ))
    network_model["settings"] = Dict{String,Any}(
        "sbase_default" => power_base,
        "vbases_default" => Dict{String,Any}(), #No default is specified for now, since default is never used
        "voltage_scale_factor" => 1E3, #Voltages are thus expressed in kV
        "sbase" => power_base,
        "power_scale_factor" => 1E6, #Power is expressed in MW
        "base_frequency" => 50.0 #Hertz
    )
    network_model["branch"] = Dict{String,Any}()
    network_model["storage"] = Dict{String,Any}()
    open(dir * config_file_name, "r") do io
        configuration_json_dict = JSON.parse(io)
    end
    #voltage_base = configuration_json_dict["gridConfig"]["basekV"]
    #power_base = configuration_json_dict["gridConfig"]["baseMVA"]
    sub_dir = splitpath(config_file_name)[1]
    configuration_json_dict
    configuration = configuration_json_dict["gridConfig"]["connection_configuration"]
    branches_file_name = sub_dir * "/" * splitpath(configuration_json_dict["gridConfig"]["branches_file"])[2]
    buses_file_name = sub_dir * "/" * splitpath(configuration_json_dict["gridConfig"]["buses_file"])[2]
    devices_file_name = sub_dir * "/" * splitpath(configuration_json_dict["gridConfig"]["devices_file"])[2]


    open(dir * buses_file_name, "r") do io
        buses_json_dict = JSON.parse(io)
        for bus in buses_json_dict
            id = bus["busId"] + 1 #Indexing starts at one in Julia (only works with integers)
            id_s = string(id)
            network_model["bus_lookup"][id_s] = id
            network_model["settings"]["vbases_default"][id_s] = voltage_base

            if id == 1 #Settings for slack bus
                network_model["bus"][id_s] = Dict{String,Any}(
                    "name" => "slack",
                    "bus_type" => 3,
                    ##"grounded"  => Bool[0, 0, 0],
                    #"terminals" => [1, 2, 3],
                    "vbase" => voltage_base,
                    "index" => id,
                    "bus_i" => id,
                    "λvmin" => 1.65, #1.03643,
                    "λvmax" => 1.65, #1.03643,
                    "vmin" => 1.0,
                    "vmax" => 1,
                    "va" => 0.0,
                    "vm" => 1,
                    #"LPp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    #"LPq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    #"LQp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    #"LQq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    #"run_type"  => 1
                )
            else
                network_model["bus"][id_s] = Dict{String,Any}(
                    "name" => id_s,
                    "bus_type" => 1,
                    #"grounded"  => Bool[0, 0, 0],
                    #"terminals" => [1, 2, 3],
                    "vbase" => voltage_base,
                    "index" => id,
                    "bus_i" => id,
                    "λvmin" => 1.65,#1.03643,
                    "λvmax" => 1.65, #1.03643,
                    "vmin" => 0.95, # free to change, hard boundary
                    "vmax" => 1.05,
                    #"LPp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    #"LPq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    #"LQp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    #"LQq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    "run_type" => 1
                )
            end
        end
    end

    #print(device_df)
    open(dir * devices_file_name, "r") do io
        devices_json_dict = JSON.parse(io)
        for device in devices_json_dict["LVcustomers"]
            id = device["deviceId"] + 1 #Indexing starts at one in Julia
            # d = device_df[in(id - 1).(device_df.dev_id), :]
            id_s = string(id)
            # μ = dist_lv_ts_feeder[in(d[!, "category"][1]).(dist_lv_ts_feeder.cluster), :][!, "lower"][1] ## @Arpan, are mu and sigma here simply lower and upper?
            # σ = dist_lv_ts_feeder[in(d[!, "category"][1]).(dist_lv_ts_feeder.cluster), :][!, "upper"][1]
            cons = convert(Float64, device["yearlyNetConsumption"])
            network_model["load"][id_s] = Dict{String,Any}(
                #"connections"   => vec(Int.(device["phases"])),
                "name" => id_s * "-" * device["coded_ean"],
                "status" => 1,
                "vbase" => voltage_base,
                "vnom_kv" => 1.0,
                "source_id" => device["coded_ean"],
                "load_bus" => device["busId"] + 1,
                "dispatchable" => 0,
                "index" => id,
                "yearlyNetConsumption" => cons,
                #"phases"        => device["phases"],
                "pd" => 0.1 / 1e3 / power_base / 3,
                "qd" => 0.01 / 1e3 / power_base / 3 / 10,
                # "pd" => max(0.1, μ) / 1e3 / power_base / 3,
                # "qd" => max(0.01, μ) / 1e3 / power_base / 3 / 10,
                "p_inj" => 0.0,
                "q_inj" => 0.0,
                "conn_cap_kW" => device["connectionCapacity"],
                "dst_id" => 0, # d[!, "category"][1],
                "cluster_id" => 0, #findall(x -> x == 1, [s_dict["$i"]["dst_id"] == d[!, "category"][1] for i = 1:length(s_dict)])[1],
                "μ" => 0, #μ,
                "σ" => 0) # σ)
        end
    end

    open(dir * branches_file_name, "r") do io
        branches_json_dict = JSON.parse(io)
        impedance_dict = Dict{String,Any}( # Units are ohm/km!!!
            "BT - Desconocido BT" => [0.21, 0.075],
            "BT - MANGUERA" => [0.3586, 0.089], #[1.23, 0.08],
            "BT - RV 0,6/1 KV 2*16 KAL" => [2.14, 0.09], #16 = 20 A
            "BT - RV 0,6/1 KV 2*25 KAL" => [1.34, 0.097],
            "BT - RV 0,6/1 KV 3(1*150 KAL) + 1*95 KAL" => [0.2309, 0.085],
            "BT - RV 0,6/1 KV 3(1*240 KAL) + 1*150 KAL" => [0.1602, 0.079],
            "BT - RV 0,6/1 KV 3(1*240 KAL) + 1*95 KAL" => [0.1602, 0.079],
            "BT - RV 0,6/1 KV 4*25 KAL" => [1.34, 0.097],
            "BT - RV 0,6/1 KV 4*50 KAL" => [0.71849, 0.093],
            "BT - RV 0,6/1 KV 4*95 KAL" => [0.3586, 0.089],
            "BT - RX 0,6/1 KV 2*16 Cu" => [1.23, 0.08],
            "BT - RX 0,6/1 KV 2*2 Cu" => [9.9, 0.075],
            "BT - RX 0,6/1 KV 2*4 Cu" => [4.95, 0.075],
            "BT - RX 0,6/1 KV 2*6 Cu" => [3.3, 0.075],
            "BT - RZ 0,6/1 KV 2*16 AL" => [2.14, 0.09],
            "BT - RZ 0,6/1 KV 3*150 AL/80 ALM" => [0.2309, 0.85],
            "BT - RZ 0,6/1 KV 3*150 AL/95 ALM" => [0.2309, 0.85],
            "BT - RZ 0,6/1 KV 3*25 AL/54,6 ALM" => [1.34, 0.097],
            "BT - RZ 0,6/1 KV 3*35 AL/54,6 ALM" => [0.9073, 0.095],
            "BT - RZ 0,6/1 KV 3*50 AL/54,6 ALM" => [0.718497, 0.093],
            "BT - RZ 0,6/1 KV 3*70 ALM/54,6 AL" => [0.4539, 0.091],
            "BT - RZ 0,6/1 KV 3*95 AL/54,6 ALM" => [0.3586, 0.089],
            "BT - RZ 0,6/1 KV 4*16 AL" => [2.14, 0.09],
            "aansluitkabel" => [1.15, 0.150],
            "unknown" => [0.01, 0.001]
        )

        open(dir * branches_file_name, "r") do io
            branches_json_dict = JSON.parse(io)
            currentmax_dict = Dict{String,Any}( # what are the units? amperes i guess
                "BT - Desconocido BT" => 200,
                "BT - MANGUERA" => 150, #200 certain  40.18#150
                "BT - RV 0,6/1 KV 2*16 KAL" => 75,
                "BT - RV 0,6/1 KV 2*25 KAL" => 100,
                "BT - RV 0,6/1 KV 3(1*150 KAL) + 1*95 KAL" => 305,
                "BT - RV 0,6/1 KV 3(1*240 KAL) + 1*150 KAL" => 344,
                "BT - RV 0,6/1 KV 3(1*240 KAL) + 1*95 KAL" => 344,
                "BT - RV 0,6/1 KV 4*25 KAL" => 100,
                "BT - RV 0,6/1 KV 4*50 KAL" => 150,
                "BT - RV 0,6/1 KV 4*95 KAL" => 230,
                "BT - RX 0,6/1 KV 2*16 Cu" => 75,#95,
                "BT - RX 0,6/1 KV 2*2 Cu" => 40, #30,
                "BT - RX 0,6/1 KV 2*4 Cu" => 60,#40,
                "BT - RX 0,6/1 KV 2*6 Cu" => 80, #50,
                "BT - RZ 0,6/1 KV 2*16 AL" => 20, #75,
                "BT - RZ 0,6/1 KV 3*150 AL/80 ALM" => 305, #264,
                "BT - RZ 0,6/1 KV 3*150 AL/95 ALM" => 305,#264,
                "BT - RZ 0,6/1 KV 3*25 AL/54,6 ALM" => 100, #78.98,
                "BT - RZ 0,6/1 KV 3*35 AL/54,6 ALM" => 120,
                "BT - RZ 0,6/1 KV 3*50 AL/54,6 ALM" => 150, #118.47,
                "BT - RZ 0,6/1 KV 3*70 ALM/54,6 AL" => 160,
                "BT - RZ 0,6/1 KV 3*95 AL/54,6 ALM" => 230, # 182.21,
                "BT - RZ 0,6/1 KV 4*16 AL" => 75,
                "aansluitkabel" => 120, #200, 
                "unknown" => 400 # i put something high so it doesn't restrict the solution
            )


            for branch in branches_json_dict
                id = branch["branchId"] + 1
                id_s = string(id)
                network_model["branch"][id_s] = Dict{String,Any}(
                    "shift" => 0.0,
                    #"f_connections" => [1, 2, 3],
                    "name" => id_s,
                    "switch" => false,
                    "g_to" => 0.0,
                    "c_rating_a" => currentmax_dict[branch["cableType"]],
                    "vbase" => voltage_base,
                    "g_fr" => 0.0,
                    #"t_connections" => [1, 2, 3],
                    "f_bus" => branch["upBusId"] + 1,
                    "b_fr" => 0.0,
                    "c_rating_b" => currentmax_dict[branch["cableType"]],
                    "br_status" => 1,
                    "t_bus" => branch["downBusId"] + 1,
                    "b_to" => 0.0,
                    "index" => id,
                    "angmin" => -1.0472,
                    "angmax" => 1.0472,
                    "transformer" => false,
                    "tap" => 1.0,
                    "c_rating_c" => currentmax_dict[branch["cableType"]],
                    "λcmax" => 1.65 #1.65 #2.5 #1.03643    
                )

                if haskey(impedance_dict, branch["cableType"])
                    # network_model["branch"][id_s]["br_x"] = (cross_area_fact^(1/4)).* impedance_dict[branch["cableType"]][2] .* (branch["cableLength"]/1000+1E-6)./  Z_base 
                    # network_model["branch"][id_s]["br_r"] = (1/cross_area_fact).* impedance_dict[branch["cableType"]][1] .* (branch["cableLength"]/1000+1E-6)./  Z_base
                    network_model["branch"][id_s]["br_x"] = (1 / cross_area_fact) .* impedance_dict[branch["cableType"]][2] .* (branch["cableLength"] / 1000 + 1E-6) ./ Z_base
                    network_model["branch"][id_s]["br_r"] = (1 / cross_area_fact) .* impedance_dict[branch["cableType"]][1] .* (branch["cableLength"] / 1000 + 1E-6) ./ Z_base
                end


                if haskey(currentmax_dict, branch["cableType"])
                    network_model["branch"][id_s]["rate_a"] = cross_area_fact .* ((currentmax_dict[branch["cableType"]] * voltage_base) / 1e3) / power_base

                    #network_model["branch"][id_s]["I_rating"] = currentmax_dict[branch["cableType"]]/current_base

                end
            end
        end
    end

    network_model["sdata"] = s_dict
    network_model["curt"] = curt

    network_model["PV"] = deepcopy(network_model["load"])
    [network_model["PV"][d]["μ"] = s_dict[string(length(s_dict))]["pc"] for d in keys(network_model["PV"])]
    [network_model["PV"][d]["σ"] = s_dict[string(length(s_dict))]["pd"] for d in keys(network_model["PV"])]
    [network_model["PV"][d]["pd"] = s_dict[string(length(s_dict))]["pd"] / 1e6 / power_base / 3 for d in keys(network_model["PV"])]
    return network_model
end;

"Converts JSON file of three phase DN to single phase equivalent"
function build_mathematical_model_single_phase(dir, config_file_name, load_dist_csv, pv_dist_csv; t_s=52, pd=0.0, qd=0.0, scale_factor=1.0, curt=0.0, cross_area_fact=1.0)
    #configuration = "star"

    # build_mathematical_model_single_phase(file, feeder,load_file, pv_file, t_s= 59)
    # file  = joinpath(BASE_DIR, "test/data/Spanish/") --> dir
    # feeder="All_feeder/"*b.conf --> config file name
    # load_file= "beta_lm_2016_8_6.csv" --> load_dist_csv
    # pv_file = "beta_pm_2016_8_6.csv" --> pv_dist_cv

    """
    Specify voltage and power base, voltage base should be the phase-to-ground voltage
    of the feeder studied (kV), the power base can be arbitrairly chosen (MW)
    """
    voltage_base = 0.230  # (kV)
    power_base = 0.5  # (MW)
    Z_base = voltage_base^2 / power_base # (Ohm)
    current_base = power_base / (voltage_base * 1e-3) # (A)

    mwpu = 1 / power_base
    kwpu = (1e-3) / power_base

    network_model = Dict{String,Any}()
    configuration_json_dict = Dict{Any,Any}()
    device_df=CSV.read(dir*config_file_name[1:length(config_file_name)-19]*".csv", DataFrame) #csv file with categories.

    dist_lv = CSV.read(dir * load_dist_csv, DataFrame) # load_dist_csv
    dist_pv = CSV.read(dir * pv_dist_csv, DataFrame)
    # dist_pv=CSV.read(dir*"beta_pm_2022_181"*".csv", DataFrame)
    dist_pv_ts = dist_pv[in([t_s]).(dist_pv.timeslot), :]
    dist_lv_ts = dist_lv[in([t_s]).(dist_lv.timeslot), :]
    # requires the csv file with categories.

    global dist_lv_ts_feeder = dist_lv_ts[in(unique(device_df.category)).(dist_lv_ts.cluster), :]
    global s_dict = Dict()
    i = 1
    for dist in eachrow(dist_lv_ts_feeder)
        s = Dict()
        s["dst"] = "Beta"
        s["dst_id"] = dist["cluster"]
        s["pa"] = dist["alpha"]
        s["pb"] = dist["beta"]
        s["pc"] = dist["lower"]
        s["pd"] = dist["lower"] + dist["upper"]
        s_dict[string(i)] = s
        i = i + 1
    end


    ##add Irradiance if day time or there is some Irradiance
    if dist_pv_ts.upper[1] > 0
        s = Dict()
        s["dst"] = "Beta"
        s["dst_id"] = 55
        s["pa"] = dist_pv_ts[!, "alpha"][1]
        s["pb"] = dist_pv_ts[!, "beta"][1]
        s["pc"] = dist_pv_ts[!, "lower"][1]
        s["pd"] = dist_pv_ts[!, "lower"][1] + dist_pv_ts[!, "upper"][1]
        s_dict[string(i)] = s
    end


    #network_model["is_kron_reduced"] = true
    network_model["p_factor"] = 0.95
    network_model["q_factor"] = sqrt(1 - 0.95^2)
    network_model["dcline"] = Dict{String,Any}()
    network_model["switch"] = Dict{String,Any}()
    #network_model["is_projected"] = true
    network_model["per_unit"] = true
    #network_model["data_model"] = MATHEMATICAL
    network_model["shunt"] = Dict{String,Any}()
    network_model["transformer"] = Dict{String,Any}()
    network_model["bus"] = Dict{String,Any}()
    network_model["map"] = Dict{String,Any}()
    #network_model["conductors"] = 1
    network_model["baseMVA"] = power_base
    network_model["basekv"] = voltage_base
    network_model["bus_lookup"] = Dict{Any,Int64}()
    network_model["run_type"] = 1
    network_model["load"] = Dict{String,Any}()
    network_model["gen"] = Dict{String,Any}("1" => Dict{String,Any}(
        "pg" => 0.2,
        "model" => 2,
        #"connections"   => [1, 2, 3],
        "shutdown" => 0.0,
        "startup" => 0.0,
        #"configuration" => WYE,
        "name" => "virtual_generator",
        "qg" => 0.0,
        "gen_bus" => 1, # always starting from 1 for slackbus
        "vbase" => voltage_base,
        "source_id" => Any["gen", 1],
        "index" => 1,
        "cost" => [20000.0, 1400.0, 0.0], #?? (@Arpan, what is this cost?)
        "gen_status" => 1,
        "qmax" => 1.275,
        "qmin" => -1.275,
        "pmax" => 1.5,
        "pmin" => -1.5,
        "ncost" => 3,
        "λpmin" => 1.65, #1.03643 ,
        "λpmax" => 1.65, #1.03643 ,
        "λqmin" => 1.65, #1.03643 ,
        "λqmax" => 1.65 #1.03643
    ))
    network_model["settings"] = Dict{String,Any}(
        "sbase_default" => power_base,
        "vbases_default" => Dict{String,Any}(), #No default is specified for now, since default is never used
        "voltage_scale_factor" => 1E3, #Voltages are thus expressed in kV
        "sbase" => power_base,
        "power_scale_factor" => 1E6, #Power is expressed in MW
        "base_frequency" => 50.0 #Hertz
    )
    network_model["branch"] = Dict{String,Any}()
    network_model["storage"] = Dict{String,Any}()
    open(dir * config_file_name, "r") do io
        configuration_json_dict = JSON.parse(io)
    end
    #voltage_base = configuration_json_dict["gridConfig"]["basekV"]
    #power_base = configuration_json_dict["gridConfig"]["baseMVA"]
    sub_dir = splitpath(config_file_name)[1]
    configuration_json_dict
    configuration = configuration_json_dict["gridConfig"]["connection_configuration"]
    branches_file_name = sub_dir * "/" * splitpath(configuration_json_dict["gridConfig"]["branches_file"])[2]
    buses_file_name = sub_dir * "/" * splitpath(configuration_json_dict["gridConfig"]["buses_file"])[2]
    devices_file_name = sub_dir * "/" * splitpath(configuration_json_dict["gridConfig"]["devices_file"])[2]


    open(dir * buses_file_name, "r") do io
        buses_json_dict = JSON.parse(io)
        for bus in buses_json_dict
            id = bus["busId"] + 1 #Indexing starts at one in Julia (only works with integers)
            id_s = string(id)
            network_model["bus_lookup"][id_s] = id
            network_model["settings"]["vbases_default"][id_s] = voltage_base

            if id == 1 #Settings for slack bus
                network_model["bus"][id_s] = Dict{String,Any}(
                    "name" => "slack",
                    "bus_type" => 3,
                    ##"grounded"  => Bool[0, 0, 0],
                    #"terminals" => [1, 2, 3],
                    "vbase" => voltage_base,
                    "index" => id,
                    "bus_i" => id,
                    "λvmin" => 1.65, #1.03643,
                    "λvmax" => 1.65, #1.03643,
                    "vmin" => 1.0,
                    "vmax" => 1,
                    "va" => 0.0,
                    "vm" => 1,
                    #"LPp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    #"LPq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    #"LQp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    #"LQq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    #"run_type"  => 1
                )
            else
                network_model["bus"][id_s] = Dict{String,Any}(
                    "name" => id_s,
                    "bus_type" => 1,
                    #"grounded"  => Bool[0, 0, 0],
                    #"terminals" => [1, 2, 3],
                    "vbase" => voltage_base,
                    "index" => id,
                    "bus_i" => id,
                    "λvmin" => 1.65,#1.03643,
                    "λvmax" => 1.65, #1.03643,
                    "vmin" => 0.95, # free to change, hard boundary
                    "vmax" => 1.05,
                    #"LPp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    #"LPq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    #"LQp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    #"LQq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    "run_type" => 1
                )
            end
        end
    end

    #print(device_df)
    open(dir * devices_file_name, "r") do io
        devices_json_dict = JSON.parse(io)
        for device in devices_json_dict["LVcustomers"]
            id = device["deviceId"] + 1 #Indexing starts at one in Julia
            global d = device_df[in(id - 1).(device_df.dev_id), :]
            id_s = string(id)
            μ = dist_lv_ts_feeder[in(d[!, "category"][1]).(dist_lv_ts_feeder.cluster), :][!, "lower"][1] ## @Arpan, are mu and sigma here simply lower and upper?
            σ = dist_lv_ts_feeder[in(d[!, "category"][1]).(dist_lv_ts_feeder.cluster), :][!, "upper"][1]
            cons = convert(Float64, device["yearlyNetConsumption"])
            network_model["load"][id_s] = Dict{String,Any}(
                #"connections"   => vec(Int.(device["phases"])),
                "name" => id_s * "-" * device["coded_ean"],
                "status" => 1,
                "vbase" => voltage_base,
                "vnom_kv" => 1.0,
                "source_id" => device["coded_ean"],
                "load_bus" => device["busId"] + 1,
                "dispatchable" => 0,
                "index" => id,
                "yearlyNetConsumption" => cons,
                #"phases"        => device["phases"],
                "pd" => max(0.1, μ) / 1e3 / power_base / 3,
                "qd" => max(0.01, μ) / 1e3 / power_base / 3 / 10,
                "p_inj" => 0.0,
                "q_inj" => 0.0,
                "conn_cap_kW" => device["connectionCapacity"],
                "dst_id" => d[!, "category"][1],
                "cluster_id" => findall(x -> x == 1, [s_dict["$i"]["dst_id"] == d[!, "category"][1] for i = 1:length(s_dict)])[1],
                "μ" => μ,
                "σ" => σ)
        end
    end

    open(dir * branches_file_name, "r") do io
        branches_json_dict = JSON.parse(io)
        impedance_dict = Dict{String,Any}(
            "BT - Desconocido BT" => [0.21, 0.075],
            "BT - MANGUERA" => [0.3586, 0.089], #[1.23, 0.08],
            "BT - RV 0,6/1 KV 2*16 KAL" => [2.14, 0.09], #16 = 20 A
            "BT - RV 0,6/1 KV 2*25 KAL" => [1.34, 0.097],
            "BT - RV 0,6/1 KV 3(1*150 KAL) + 1*95 KAL" => [0.2309, 0.085],
            "BT - RV 0,6/1 KV 3(1*240 KAL) + 1*150 KAL" => [0.1602, 0.079],
            "BT - RV 0,6/1 KV 3(1*240 KAL) + 1*95 KAL" => [0.1602, 0.079],
            "BT - RV 0,6/1 KV 4*25 KAL" => [1.34, 0.097],
            "BT - RV 0,6/1 KV 4*50 KAL" => [0.71849, 0.093],
            "BT - RV 0,6/1 KV 4*95 KAL" => [0.3586, 0.089],
            "BT - RX 0,6/1 KV 2*16 Cu" => [1.23, 0.08],
            "BT - RX 0,6/1 KV 2*2 Cu" => [9.9, 0.075],
            "BT - RX 0,6/1 KV 2*4 Cu" => [4.95, 0.075],
            "BT - RX 0,6/1 KV 2*6 Cu" => [3.3, 0.075],
            "BT - RZ 0,6/1 KV 2*16 AL" => [2.14, 0.09],
            "BT - RZ 0,6/1 KV 3*150 AL/80 ALM" => [0.2309, 0.85],
            "BT - RZ 0,6/1 KV 3*150 AL/95 ALM" => [0.2309, 0.85],
            "BT - RZ 0,6/1 KV 3*25 AL/54,6 ALM" => [1.34, 0.097],
            "BT - RZ 0,6/1 KV 3*35 AL/54,6 ALM" => [0.9073, 0.095],
            "BT - RZ 0,6/1 KV 3*50 AL/54,6 ALM" => [0.718497, 0.093],
            "BT - RZ 0,6/1 KV 3*70 ALM/54,6 AL" => [0.4539, 0.091],
            "BT - RZ 0,6/1 KV 3*95 AL/54,6 ALM" => [0.3586, 0.089],
            "BT - RZ 0,6/1 KV 4*16 AL" => [2.14, 0.09],
            "aansluitkabel" => [1.15, 0.150]
        )

        open(dir * branches_file_name, "r") do io
            branches_json_dict = JSON.parse(io)
            currentmax_dict = Dict{String,Any}(
                "BT - Desconocido BT" => 200,
                "BT - MANGUERA" => 150, #200 certain  40.18#150
                "BT - RV 0,6/1 KV 2*16 KAL" => 75,
                "BT - RV 0,6/1 KV 2*25 KAL" => 100,
                "BT - RV 0,6/1 KV 3(1*150 KAL) + 1*95 KAL" => 305,
                "BT - RV 0,6/1 KV 3(1*240 KAL) + 1*150 KAL" => 344,
                "BT - RV 0,6/1 KV 3(1*240 KAL) + 1*95 KAL" => 344,
                "BT - RV 0,6/1 KV 4*25 KAL" => 100,
                "BT - RV 0,6/1 KV 4*50 KAL" => 150,
                "BT - RV 0,6/1 KV 4*95 KAL" => 230,
                "BT - RX 0,6/1 KV 2*16 Cu" => 75,#95,
                "BT - RX 0,6/1 KV 2*2 Cu" => 40, #30,
                "BT - RX 0,6/1 KV 2*4 Cu" => 60,#40,
                "BT - RX 0,6/1 KV 2*6 Cu" => 80, #50,
                "BT - RZ 0,6/1 KV 2*16 AL" => 20, #75,
                "BT - RZ 0,6/1 KV 3*150 AL/80 ALM" => 305, #264,
                "BT - RZ 0,6/1 KV 3*150 AL/95 ALM" => 305,#264,
                "BT - RZ 0,6/1 KV 3*25 AL/54,6 ALM" => 100, #78.98,
                "BT - RZ 0,6/1 KV 3*35 AL/54,6 ALM" => 120,
                "BT - RZ 0,6/1 KV 3*50 AL/54,6 ALM" => 150, #118.47,
                "BT - RZ 0,6/1 KV 3*70 ALM/54,6 AL" => 160,
                "BT - RZ 0,6/1 KV 3*95 AL/54,6 ALM" => 230, # 182.21,
                "BT - RZ 0,6/1 KV 4*16 AL" => 75,
                "aansluitkabel" => 120 #200
            )


            for branch in branches_json_dict
                id = branch["branchId"] + 1
                id_s = string(id)
                network_model["branch"][id_s] = Dict{String,Any}(
                    "shift" => 0.0,
                    #"f_connections" => [1, 2, 3],
                    "name" => id_s,
                    "switch" => false,
                    "g_to" => 0.0,
                    "c_rating_a" => currentmax_dict[branch["cableType"]],
                    "vbase" => voltage_base,
                    "g_fr" => 0.0,
                    #"t_connections" => [1, 2, 3],
                    "f_bus" => branch["upBusId"] + 1,
                    "b_fr" => 0.0,
                    "c_rating_b" => currentmax_dict[branch["cableType"]],
                    "br_status" => 1,
                    "t_bus" => branch["downBusId"] + 1,
                    "b_to" => 0.0,
                    "index" => id,
                    "angmin" => -1.0472,
                    "angmax" => 1.0472,
                    "transformer" => false,
                    "tap" => 1.0,
                    "c_rating_c" => currentmax_dict[branch["cableType"]],
                    "λcmax" => 1.65 #1.65 #2.5 #1.03643    
                )

                if haskey(impedance_dict, branch["cableType"])
                    # network_model["branch"][id_s]["br_x"] = (cross_area_fact^(1/4)).* impedance_dict[branch["cableType"]][2] .* (branch["cableLength"]/1000+1E-6)./  Z_base 
                    # network_model["branch"][id_s]["br_r"] = (1/cross_area_fact).* impedance_dict[branch["cableType"]][1] .* (branch["cableLength"]/1000+1E-6)./  Z_base
                    network_model["branch"][id_s]["br_x"] = (1 / cross_area_fact) .* impedance_dict[branch["cableType"]][2] .* (branch["cableLength"] / 1000 + 1E-6) ./ Z_base
                    network_model["branch"][id_s]["br_r"] = (1 / cross_area_fact) .* impedance_dict[branch["cableType"]][1] .* (branch["cableLength"] / 1000 + 1E-6) ./ Z_base
                end


                if haskey(currentmax_dict, branch["cableType"])
                    network_model["branch"][id_s]["rate_a"] = cross_area_fact .* ((currentmax_dict[branch["cableType"]] * voltage_base) / 1e3) / power_base

                    #network_model["branch"][id_s]["I_rating"] = currentmax_dict[branch["cableType"]]/current_base

                end
            end
        end
    end


    network_model["sdata"] = s_dict
    network_model["curt"] = curt

    network_model["PV"] = deepcopy(network_model["load"])
    [network_model["PV"][d]["μ"] = s_dict[string(length(s_dict))]["pc"] for d in keys(network_model["PV"])]
    [network_model["PV"][d]["σ"] = s_dict[string(length(s_dict))]["pd"] for d in keys(network_model["PV"])]
    [network_model["PV"][d]["pd"] = s_dict[string(length(s_dict))]["pd"] / 1e6 / power_base / 3 for d in keys(network_model["PV"])]
    return network_model
end;



# input
deg  = 2
aux  = true
r = false
dir=BASE_DIR*"/test/data/Spanish/All_feeder" #data directory
all_feeder=CSV.read(dir*"/ts_all_feeder.csv",DataFrame,header=["conf", "ts","HC1","HC2"]) 
load_file= "beta_lm_2016_8_6.csv"
pv_file = "beta_pm_2016_8_6.csv"

all_feeder[!,"HC1"]=zeros(nrow(all_feeder))
all_feeder[!,"HC2"]=zeros(nrow(all_feeder))
all_feeder[!,"HC3"]=zeros(nrow(all_feeder))
hc1=[]
hc2=[]
hc3=[]
t_cc=[]
t_opf=[]
nodes=[]
consumers=[]
unc=[]

device_ean_list = []
p_size_list = []

global i=0
for b in eachrow(all_feeder)
    feeder="All_feeder/"*b.conf
    
    ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>5, "max_cpu_time"=>500.0)
    
    s2 = Dict("output" => Dict("duals" => true))


    global i=i+1
    print("Feeder no: $i \n")
    #feeder="All_feeder/"*all_feeder[1,"conf"]
    file  = joinpath(BASE_DIR, "test/data/Spanish/")
    data  = build_mathematical_model_reconfiguration(file, feeder,load_file, pv_file, t_s= 59) # @Arpan: is only timestep 59 used regardless of which timestep was listed in the all_feeder.csv file?

    push!(unc,length(data["sdata"])) # nr of lv+pv distributions considered --> unc
    push!(nodes,length(data["bus"])) # nr of buses --> nodes
    push!(consumers,length(data["load"])) # nr of devices --> consumers
    
    [data["bus"]["$i"]["vmin"]=0.8 for i=1:length(data["bus"])] # @arpan why do we set all vmin to 0.8?

    result_hc= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=r; setting=s2, stochastic=false)
   
    #start up
    [data["PV"]["$i"]["p_size_start"]=result_hc["solution"]["PV"]["$i"]["p_size"] for i=1:length(data["PV"])]

    # # ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>5, "max_iter"=>3000)
    # result_hc_2= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=r; setting=s2)
    # e=1;
    # m=1

    # # if result_hc_2["termination_status"]== PM.LOCALLY_SOLVED
    # # for i=1:length(data["bus"])
    # #     if -result_hc_2["solution"]["nw"]["1"]["bus"]["$i"]["dual_voltage_max"]>500
    # #         l_old=data["bus"]["$i"]["λvmax"]
    # #         m=sample(result_hc_2, "bus", i, "vs"; sample_size=100000)
    # #         if quantile(m,[0.95])[1]>1.001*data["bus"]["$i"]["vmax"]^2
    # #             data["bus"]["$i"]["λvmax"]=2*l_old-(data["bus"]["$i"]["vmax"]^2-mean(m))/std(m)+0.3
    # #         elseif quantile(m,[0.95])[1]>1.0001*data["bus"]["$i"]["vmax"]^2
    # #                 data["bus"]["$i"]["λvmax"]=2*l_old-(data["bus"]["$i"]["vmax"]^2-mean(m))/std(m)+0.1
    # #         elseif quantile(m,[0.95])[1]< 0.99*data["bus"]["$i"]["vmax"]^2
    # #             data["bus"]["$i"]["λvmax"]= l_old-0.8
    # #         elseif quantile(m,[0.95])[1]< 1*data["bus"]["$i"]["vmax"]^2
    # #             data["bus"]["$i"]["λvmax"]= l_old-0.2
    # #         end
    # #     end
    # # end
    # if result_hc_2["termination_status"]== PM.LOCALLY_SOLVED # @Arpan, what is  data["bus"]["$i"]["λvmax"] ?
    #     for i=1:length(data["bus"])
    #         if -result_hc_2["solution"]["nw"]["1"]["bus"]["$i"]["dual_voltage_max"]>500
    #             l_old=data["bus"]["$i"]["λvmax"]
    #             m=sample(result_hc_2, "bus", i, "vs"; sample_size=100000)
    #             if quantile(m,[0.95])[1]>1.001*data["bus"]["$i"]["vmax"]^2
    #                 data["bus"]["$i"]["λvmax"]=l_old+0.4
    #             elseif quantile(m,[0.95])[1]>1.0001*data["bus"]["$i"]["vmax"]^2
    #                     data["bus"]["$i"]["λvmax"]=l_old+0.2
    #             elseif quantile(m,[0.95])[1]< 0.99*data["bus"]["$i"]["vmax"]^2
    #                 data["bus"]["$i"]["λvmax"]= l_old-0.6
    #             elseif quantile(m,[0.95])[1]< 1*data["bus"]["$i"]["vmax"]^2
    #                 data["bus"]["$i"]["λvmax"]= l_old-0.2
    #             end
    #         end
    #     end
    # result_hc_1= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=r; setting=s2) # @Arpan, is this the homogeneous solution?

    
    #     push!(hc1,result_hc_2["objective"]) #stochastic
    #     push!(hc2,result_hc_1["objective"]) #homogeneous?
    #     push!(t_cc, result_hc_1["solve_time"])
    # else
    #     push!(hc1,-1)
    #     push!(hc2,-1)
    #     push!(t_cc, result_hc_2["solve_time"])
    # end
     
    if result_hc["termination_status"]== PM.LOCALLY_SOLVED
        push!(hc3,result_hc["objective"]) # deterministic
        push!(t_opf, result_hc["solve_time"])
    else
        push!(hc3,-1)
        push!(t_opf, result_hc["solve_time"])
    end

    for k in keys(data["PV"])
        push!(device_ean_list, data["PV"][k]["source_id"])
        push!(p_size_list, result_hc["solution"]["PV"][k]["p_size"])
    end
end

device_psize_df = DataFrame(device_ean = device_ean_list, p_size = p_size_list)
CSV.write("PV_HC_per_device_deterministic.csv",string.(device_psize_df))


# all_feeder[!,"HC0"]=hc1 # stochastic
# all_feeder[!,"HC_CC"]=hc2 # homogeneous?
all_feeder[!,"HC_OPF"]=hc3 # deterministic
all_feeder[!,"t_opf"]=t_opf
# all_feeder[!,"t_cc"]=t_cc
all_feeder[!,"consumers"]=consumers
all_feeder[!,"nodes"]= nodes
all_feeder[!, "unc"] = unc

CSV.write("PV_HC_feeders_with_start_deterministic.csv",string.(all_feeder))

# scatter(all_feeder[all_feeder[!, :t_cc] .<500, :].unc, log10.(all_feeder[all_feeder[!, :t_cc] .<500, :].t_cc),label="gPC-CC-OPF", figsize=(28,8))
# # scatter!([1,2,3,4,5,6,7],[result_hc["solution"]["PV"]["$i"]["p_size"] for i=1:length(data["load"])],label="OPF HC", figsize=(28,8))
# plot!(xlabel="Uncertainties [-]")
# plot!(ylabel="log10 of computation time [sec]")
# plot!(title="Fig. 11: Boxplot of log_{10} of computational time for \text{gPC-\gls{cc}-\gls{opf}} in real LV feeders with respect to the number of uncertainties considered.")

"""
#deterministic

for b in eachrow(all_feeder)
    feeder="All_feeder/"*b.conf
    file  = joinpath(BASE_DIR, "test/data/Spanish/")
    data  = SPM.build_mathematical_model_single_phase(file, feeder, t_s= b.ts)

    s2 = Dict("output" => Dict("duals" => true))
    result_hc= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red; setting=s2, stochastic=false)
    e=1;
    if result_hc["termination_status"]== PM.LOCALLY_SOLVED
        push!(hc3,result_hc["objective"])
        #push!(hc2,result_hc_1["objective"])
    else
        push!(hc3,-1)
        #push!(hc2,-1)
    end
end

all_feeder[!,"HC3"]=hc3

"""