using CSV
using DataFrames
using Plots
using LaTeXStrings

grb_exp = 9
drs_exp = 2

gurobi_filepath = "Results/Experiment_$(grb_exp)/results_$(grb_exp)_Gurobi_PLS.csv"
drs_filepath = "Results/DRS_Experiment_$(drs_exp)/results_$(drs_exp)_DRS_P13_table.csv"
drs_boyd_filepath = "Results/DRS_Experiment_$(drs_exp)/results_$(drs_exp)_DRS_Boyd_P13_table.csv"
drs_fp_filepath = "Results/DRS_Experiment_$(drs_exp)/results_$(drs_exp)_DRS_FP_P13_table.csv"

df_grb = CSV.read(gurobi_filepath, DataFrame)
df_drs = CSV.read(drs_filepath, DataFrame)
df_drs_boyd = CSV.read(drs_boyd_filepath, DataFrame)
df_drs_fp = CSV.read(drs_fp_filepath, DataFrame)

df_grb = filter(row -> row[:time_mean] != -1, df_grb)
sort!(df_grb, :m)

plot(df_grb.m, df_grb.time_mean, label=L"\textrm{Gurobi}", xlabel=L"m", ylabel=L"\textrm{time(s)}")

plot!(df_drs.m, df_drs.time_mean, label=L"\textrm{DRS}^{\epsilon}")

plot!(df_drs_boyd.m, df_drs_boyd.time_mean, label=L"\textrm{DRS}^{\epsilon(r^0)}")

plot!(df_drs_fp.m, df_drs_fp.time_mean, label=L"\textrm{DRS}^{\textrm{FP}}")

# display(current())

savefig("Plots/Experiment_$(grb_exp)/exp_$(grb_exp)_time_plot.png")

plot(df_grb.m, df_grb.H_div_AMP_norm_0_mean, label=L"\textrm{Gurobi}", xlabel=L"m", ylabel=L"||\textrm{H}||_0/||\textrm{AMP}||_0")

plot!(df_drs.m, df_drs.H_div_AMP_norm_0_mean, label=L"\textrm{DRS}^{\epsilon}")

plot!(df_drs_boyd.m, df_drs_boyd.H_div_AMP_norm_0_mean, label=L"\textrm{DRS}^{\epsilon(r^0)}")

plot!(df_drs_fp.m, df_drs_fp.H_div_AMP_norm_0_mean, label=L"\textrm{DRS}^{\textrm{FP}}")

# display(current())

savefig("Plots/Experiment_$(grb_exp)/exp_$(grb_exp)_norm_0_ratio_plot.png")

plot(df_grb.m, df_grb.H_div_AMP_norm_1_mean, label=L"\textrm{Gurobi}", xlabel=L"m", ylabel=L"||\textrm{H}||_1/||\textrm{AMP}||_1")

plot!(df_drs.m, df_drs.H_div_AMP_norm_1_mean, label=L"\textrm{DRS}^{\epsilon}")

plot!(df_drs_boyd.m, df_drs_boyd.H_div_AMP_norm_1_mean, label=L"\textrm{DRS}^{\epsilon(r^0)}")

plot!(df_drs_fp.m, df_drs_fp.H_div_AMP_norm_1_mean, label=L"\textrm{DRS}^{\textrm{FP}}")

# display(current())

savefig("Plots/Experiment_$(grb_exp)/exp_$(grb_exp)_norm_1_ratio_plot.png")