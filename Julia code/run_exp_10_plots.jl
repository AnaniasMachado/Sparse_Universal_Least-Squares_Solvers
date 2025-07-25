using CSV
using DataFrames
using Plots

grb_exp = 10
drs_exp = 1

gurobi_filepath = "Results/Experiment_$(grb_exp)/results_$(grb_exp)_Gurobi_Cal_P123.csv"
drs_filepath = "Results/DRS_Experiment_$(drs_exp)/results_$(drs_exp)_DRS_P123_table.csv"
drs_boyd_filepath = "Results/DRS_Experiment_$(drs_exp)/results_$(drs_exp)_DRS_Boyd_P123_table.csv"
drs_fp_filepath = "Results/DRS_Experiment_$(drs_exp)/results_$(drs_exp)_DRS_FP_P123_table.csv"
admme_filepath = "Results/DRS_Experiment_$(drs_exp)/results_$(drs_exp)_ADMMe_P123_table.csv"
admm_fp_filepath = "Results/DRS_Experiment_$(drs_exp)/results_$(drs_exp)_ADMM_P123_table.csv"

df_grb = CSV.read(gurobi_filepath, DataFrame)
df_drs = CSV.read(drs_filepath, DataFrame)
df_drs_boyd = CSV.read(drs_boyd_filepath, DataFrame)
df_drs_fp = CSV.read(drs_fp_filepath, DataFrame)
df_admme = CSV.read(admme_filepath, DataFrame)
df_admm= CSV.read(admm_fp_filepath, DataFrame)

df_grb = filter(row -> row[:time_mean] != -1, df_grb)
sort!(df_grb, :m)

plot(df_grb.m, df_grb.time_mean, label="Gurobi", xlabel="m", ylabel="time(s)", title="m x time")

plot!(df_drs.m, df_drs.time_mean, label="DRS")

plot!(df_drs_boyd.m, df_drs_boyd.time_mean, label="DRS_Boyd")

plot!(df_drs_fp.m, df_drs_fp.time_mean, label="DRS_FP")

plot!(df_admme.m, df_admme.time_mean, label="ADMM^epsilon")

plot!(df_admm.m, df_admm.time_mean, label="ADMM")

# display(current())

savefig("Plots/Experiment_$(grb_exp)/exp_$(grb_exp)_time_plot.png")

plot(df_grb.m, df_grb.H_div_AMP_norm_0_mean, label="Gurobi", xlabel="m", ylabel="||H||_0/||AMP||_0", title="m x ||H||_0/||AMP||_0")

plot!(df_drs.m, df_drs.H_div_AMP_norm_0_mean, label="DRS")

plot!(df_drs_boyd.m, df_drs_boyd.H_div_AMP_norm_0_mean, label="DRS_Boyd")

plot!(df_drs_fp.m, df_drs_fp.H_div_AMP_norm_0_mean, label="DRS_FP")

plot!(df_admme.m, df_admme.H_div_AMP_norm_0_mean, label="ADMM^epsilon")

plot!(df_admm.m, df_admm.H_div_AMP_norm_0_mean, label="ADMM")

# display(current())

savefig("Plots/Experiment_$(grb_exp)/exp_$(grb_exp)_norm_0_ratio_plot.png")

plot(df_grb.m, df_grb.H_div_AMP_norm_1_mean, label="Gurobi", xlabel="m", ylabel="||H||_1/||AMP||_1", title="m x ||H||_1/||AMP||_1")

plot!(df_drs.m, df_drs.H_div_AMP_norm_1_mean, label="DRS")

plot!(df_drs_boyd.m, df_drs_boyd.H_div_AMP_norm_1_mean, label="DRS_Boyd")

plot!(df_drs_fp.m, df_drs_fp.H_div_AMP_norm_1_mean, label="DRS_FP")

plot!(df_admme.m, df_admme.H_div_AMP_norm_1_mean, label="ADMM^epsilon")

plot!(df_admm.m, df_admm.H_div_AMP_norm_1_mean, label="ADMM")

# display(current())

savefig("Plots/Experiment_$(grb_exp)/exp_$(grb_exp)_norm_1_ratio_plot.png")