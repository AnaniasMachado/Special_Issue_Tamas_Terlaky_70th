using CSV
using DataFrames
using Plots
using LaTeXStrings

problem = "P123"
time_begin = 51
time_end = 2001
start = 1
step = 50

admm_filepath = "results/problem_$(problem)/run_data_$(problem)_ADMM.csv"
drs_opt_filepath = "results/problem_$(problem)/run_data_$(problem)_DRS_Opt.csv"
drs_fp_filepath  = "results/problem_$(problem)/run_data_$(problem)_DRS_FP.csv"

df_admm = CSV.read(admm_filepath, DataFrame)
df_drs_opt = CSV.read(drs_opt_filepath, DataFrame)
df_drs_fp  = CSV.read(drs_fp_filepath, DataFrame)

df_admm = df_admm[(df_admm.time .>= time_begin) .& (df_admm.time .<= time_end), :]
df_drs_opt = df_drs_opt[(df_drs_opt.time .>= time_begin) .& (df_drs_opt.time .<= time_end), :]
df_drs_fp  = df_drs_fp[(df_drs_fp.time .>= time_begin) .& (df_drs_fp.time .<= time_end), :]

df_admm = df_admm[start:step:end, :]
df_drs_opt = df_drs_opt[start:step:end, :]
df_drs_fp  = df_drs_fp[start:step:end, :]

default(markersize=3)

scatter(df_drs_fp.time, df_drs_fp.H_norm_0, label=L"\textrm{DRS_{fp}}", xlabel=L"\textrm{time}", ylabel=L"\Vert\!\!\textrm{H}\Vert_0")
scatter!(df_admm.time, df_admm.H_norm_0, label=L"\textrm{ADMM}")
savefig("Plots/problem_$(problem)/problem_$(problem)_H_norm_0_sliced.png")

scatter(df_drs_fp.time, df_drs_fp.H_norm_1, label=L"\textrm{DRS_{fp}}", xlabel=L"\textrm{time}", ylabel=L"\Vert\!\!\textrm{H}\Vert_1")
scatter!(df_admm.time, df_admm.H_norm_1, label=L"\textrm{ADMM}")
savefig("Plots/problem_$(problem)/problem_$(problem)_H_norm_1_sliced.png")