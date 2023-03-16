using DataFrames
using CSV
using Plots



device_psize_default_df=CSV.read("PV_HC_per_device_deterministic_default.csv",DataFrame,header=1) 
device_psize_sc2_df=CSV.read("PV_HC_per_device_deterministic_sc2.csv", DataFrame,header=1) 
# I discovered that some device names occur twice per dataframe, i don't understand why.
# this is the reason that the innerjoin has more lines than the separate dataframes.
# I'm going to ignore this, it's normally not going to hugely effect the outcome. 

# maybe take the inner set? DT[in([1,4]).(DT.ID), :]
function dropduplicates(df, cols; keep = :first)
    keep in [:first, :last] || error("keep parameter should be :first or :last")
    combine(groupby(df, cols)) do sdf
        if nrow(sdf) == 1 
            DataFrame()
        else
            DataFrame(
              filter(
                r->rownumber(r)==(keep == :first ? 1 : nrow(sdf)), 
                eachrow(sdf)
              )
            )
        end
    end
end


#See what are the duplicates
#filter(x->x.nrow!=1,combine(groupby(device_psize_default_df,:device_ean), nrow))
#To see the range of values in between in the device id_s
#filter(x->x.nrow!=1,combine(groupby(device_psize_default_df,:device_ean), :p_size=>maximum, :p_size=>minimum, nrow))

function keep_unique(p_df)
   device_ids_that_appear_once = filter(x->x.nrow==1,combine(groupby(p_df,:device_ean), nrow))[!,:device_ean]
   return filter(x->x.device_ean in device_ids_that_appear_once,p_df)
end

device_psize_compare_df = innerjoin(keep_unique(device_psize_default_df), keep_unique(device_psize_sc2_df), on = :device_ean, renamecols = "_default" => "_sc2")






HC_improvement_fractions = device_psize_compare_df.p_size_sc2 .- device_psize_compare_df.p_size_default


# plot cdf
n = length(HC_improvement_fractions)

p = plot(sort(HC_improvement_fractions), (1:n)./n, 
    xlabel = "PV HC improvement of individual device", ylabel = "fraction of devices", 
    title = "Cumulative Distribution", label = "")