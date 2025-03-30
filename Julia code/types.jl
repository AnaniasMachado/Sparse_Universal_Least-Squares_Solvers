# Ref. for obtaining peak memory usage: https://github.com/JuliaLang/julia/blob/master/test/netload/memtest.jl

using JuMP

struct DataInst
    A::Matrix{Float64}
    m::Int64
    n::Int64
    r::Int64
end

mutable struct GurobiInst
    model::JuMP.Model
    H::Matrix{VariableRef}
end

# struct DRSProjData
#     R::Matrix{Float64}
#     S::Matrix{Float64}
#     T::Matrix{Float64}
#     RMTSM::Matrix{Float64}
#     SSMP::Matrix{Float64}
#     RMPR::Matrix{Float64}
#     VG::Matrix{Float64}
# end

struct DRSProjDataSimple
    R::Matrix{Float64}
    S::Matrix{Float64}
    T::Matrix{Float64}
    RMP::Matrix{Float64}
    SMP::Matrix{Float64}
    T_factor::Matrix{Float64}
    U::Matrix{Float64}
end

struct DRSProjDataP123
    AMP::Matrix{Float64}
end

struct RUsage
    ru_utime_sec::Clong         #  user CPU time used
    ru_utime_usec::Clong        #  user CPU time used
    ru_stime_sec::Clong         #  system CPU time used
    ru_stime_usec::Clong        #  system CPU time used
    ru_maxrss::Clong            #  maximum resident set size
    ru_ixrss::Clong             #  integral shared memory size
    ru_idrss::Clong             #  integral unshared data size
    ru_isrss::Clong             #  integral unshared stack size
    ru_minflt::Clong            #  page reclaims (soft page faults)
    ru_majflt::Clong            #  page faults (hard page faults)
    ru_nswap::Clong             #  swaps
    ru_inblock::Clong           #  block input operations
    ru_oublock::Clong           #  block output operations
    ru_msgsnd::Clong            #  IPC messages sent
    ru_msgrcv::Clong            #  IPC messages received
    ru_nsignals::Clong          #  signals received
    ru_nvcsw::Clong             #  voluntary context switches
    ru_nivcsw::Clong            #  involuntary context switches
end