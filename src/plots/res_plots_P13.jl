using CSV
using DataFrames
using Plots
using LaTeXStrings

problem = "P13"

drs_opt_eps_filepath = "results/problem_$(problem)/res_data_$(problem)_DRS_Opt_Eps.csv"
drs_opt_r0_filepath = "results/problem_$(problem)/res_data_$(problem)_DRS_Opt_r0.csv"
drs_fp_eps_filepath = "results/problem_$(problem)/res_data_$(problem)_DRS_FP_Eps.csv"
drs_fp_r0_filepath = "results/problem_$(problem)/res_data_$(problem)_DRS_FP_r0.csv"

df_drs_opt_eps = CSV.read(drs_opt_eps_filepath, DataFrame)
df_drs_opt_r0 = CSV.read(drs_opt_r0_filepath, DataFrame)
df_drs_fp_eps = CSV.read(drs_fp_eps_filepath, DataFrame)
df_drs_fp_r0 = CSV.read(drs_fp_r0_filepath, DataFrame)

x = 2500

plot(x:nrow(df_drs_opt_eps), df_drs_opt_eps.pri_res[x:end], label=L"\textrm{r}_p^k", xlabel=L"k", ylabel=L"\textrm{r}_p^k")

savefig("Plots/problem_$(problem)/problem_$(problem)_DRS_Opt_Eps_pri_res_from_k=$(x).png")

plot(x:nrow(df_drs_opt_eps), df_drs_opt_eps.dual_res[x:end], label=L"\textrm{r}_d^k", xlabel=L"k", ylabel=L"\textrm{r}_d^k")

savefig("Plots/problem_$(problem)/problem_$(problem)_DRS_Opt_Eps_dual_res_from_k=$(x).png")

plot(x:nrow(df_drs_opt_eps), df_drs_opt_eps.Hh_norm[x:end], label=L"\textrm{H}^{k+1/2}", xlabel=L"k", ylabel=L"normF")

plot!(x:nrow(df_drs_opt_eps), df_drs_opt_eps.H_norm[x:end], label=L"\\textrm{H}^{k}")

savefig("Plots/problem_$(problem)/problem_$(problem)_DRS_Opt_Eps_H_norm_from_k=$(x).png")

plot(x:nrow(df_drs_opt_r0), df_drs_opt_r0.res[x:end], label=L"\textrm{r}^k", xlabel=L"k", ylabel=L"\textrm{r}^k")

savefig("Plots/problem_$(problem)/problem_$(problem)_DRS_Opt_r0_res_from_k=$(x).png")

plot(x:nrow(df_drs_fp_eps), df_drs_fp_eps.res[x:end], label=L"\tau^k", xlabel=L"k", ylabel=L"\tau^k")

savefig("Plots/problem_$(problem)/problem_$(problem)_DRS_FP_Eps_res_from_k=$(x).png")

plot(x:nrow(df_drs_fp_r0), df_drs_fp_r0.res[x:end], label=L"\tau^k", xlabel=L"k", ylabel=L"\tau^k")

savefig("Plots/problem_$(problem)/problem_$(problem)_DRS_FP_r0_res_from_k=$(x).png")