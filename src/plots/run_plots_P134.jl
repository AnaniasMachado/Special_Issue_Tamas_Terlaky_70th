using CSV
using DataFrames
using Plots
using LaTeXStrings

problem = "P134"
start = 250
step = 50

admm_filepath = "results/problem_$(problem)/run_data_$(problem)_ADMM.csv"
drs_opt_filepath = "results/problem_$(problem)/run_data_$(problem)_DRS_Opt.csv"
drs_fp_filepath  = "results/problem_$(problem)/run_data_$(problem)_DRS_FP.csv"

df_admm = CSV.read(admm_filepath, DataFrame)
df_drs_opt = CSV.read(drs_opt_filepath, DataFrame)
df_drs_fp  = CSV.read(drs_fp_filepath, DataFrame)

df_admm = df_admm[start:step:end, :]
df_drs_opt = df_drs_opt[start:step:end, :]
df_drs_fp  = df_drs_fp[start:step:end, :]

default(markersize=3)

scatter(df_drs_fp.time, df_drs_fp.H_norm_0, label=L"\textrm{DRS}_{FP}", xlabel=L"\textrm{time}", ylabel=L"||\textrm{H}||_0")
scatter!(df_admm.time, df_admm.H_norm_0, label=L"\textrm{ADMM}")
savefig("Plots/problem_$(problem)/problem_$(problem)_H_norm_0.png")

scatter(df_drs_fp.time, df_drs_fp.H_norm_1, label=L"\textrm{DRS}_{FP}", xlabel=L"\textrm{time}", ylabel=L"||\textrm{H}||_1")
scatter!(df_admm.time, df_admm.H_norm_1, label=L"\textrm{ADMM}")
savefig("Plots/problem_$(problem)/problem_$(problem)_H_norm_1.png")

scatter(df_admm.time, df_admm.pri_res, label=L"\textrm{ADMM}", xlabel=L"\textrm{time}", ylabel=L"\textrm{r}_p")
savefig("Plots/problem_$(problem)/problem_$(problem)_ADMM_pri_res.png")

scatter(df_admm.time, df_admm.dual_res, label=L"\textrm{ADMM}", xlabel=L"\textrm{time}", ylabel=L"\textrm{r}_d")
savefig("Plots/problem_$(problem)/problem_$(problem)_ADMM_dual_res.png")

scatter(df_drs_fp.time, df_drs_fp.res, label=L"\textrm{DRS}_{FP}", xlabel=L"\textrm{time}", ylabel=L"\textrm{fp}_{res}")
savefig("Plots/problem_$(problem)/problem_$(problem)_DRS_fp_res.png")

scatter(df_drs_fp.time, df_drs_opt.pri_res, label=L"\textrm{DRS}_{FP}", xlabel=L"\textrm{time}", ylabel=L"\textrm{r}_p")
savefig("Plots/problem_$(problem)/problem_$(problem)_DRS_pri_res.png")

scatter(df_drs_fp.time, df_drs_opt.dual_res, label=L"\textrm{DRS}_{FP}", xlabel=L"\textrm{time}", ylabel=L"\textrm{r}_d")
savefig("Plots/problem_$(problem)/problem_$(problem)_DRS_dual_res.png")