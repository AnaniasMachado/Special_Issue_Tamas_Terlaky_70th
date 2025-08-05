using CSV
using DataFrames
using MAT

include("../types.jl")
include("../utility.jl")
include("../methods/drs.jl")

matrix_folder = "./instances/rectangular_dense_01"
mat_files = readdir(matrix_folder)
mat_file = mat_files[1]
mat_path = joinpath(matrix_folder, mat_file)
mat_data = matread(mat_path)
A = mat_data["matrix"]
A = Matrix(A)
AMP = pinv(A)
m, n = size(A)
r = rank(A)

println("m = $m, n = $n, r = $r")

data = DataInst(A, m, n, r, AMP=AMP)
constraints = ["P13"]
problem = "P13"
rho = 3.0
lambda = 10^(-2)

epsilon = 10^(-5)
eps_opt = epsilon
eps_abs = epsilon
eps_rel = 10^(-4)
fixed_tol = false

time_limit = 1200
stop_crits = ["Opt", "Fixed_Point"]
stop_crit = stop_crits[2]

results_folder = "results/problem_$problem"

DRS_time = @elapsed begin
    DRS_H, DRS_k, DRS_res_data, DRS_sol_data = drs_res_data(A, lambda, eps_abs, eps_rel, problem, fixed_tol, eps_opt, stop_crit, time_limit)
end
DRS_H_norm_0 = matrix_norm_0(DRS_H)
DRS_H_norm_1 = norm(DRS_H, 1)

println("Stop crit: $stop_crit")
println("Fixedtol: $fixed_tol")
println("DRS time: $DRS_time")
println("DRS k: $DRS_k")
println("DRS norm 1: $DRS_H_norm_1")
println("DRS norm 0: $DRS_H_norm_0")
println("DRS bound ratio: $(DRS_H_norm_0 / (m * r))")
println("Feasibility: $(norm(A' * A * DRS_H - A'))")

df = DataFrame()

for i in 1:length(DRS_res_data)
    if fixed_tol && stop_crit == "Opt"
        result = DataFrame(
            pri_res = [DRS_res_data[i][1]],
            dual_res = [DRS_res_data[i][2]],
            Hh_norm = [norm(DRS_sol_data[i][1])],
            H_norm = [norm(DRS_sol_data[i][2])]
        )
        append!(df, result)
    else
        result = DataFrame(
            res = [DRS_res_data[i]],
            Hh_norm = [norm(DRS_sol_data[i][1])],
            H_norm = [norm(DRS_sol_data[i][2])]
        )
        append!(df, result)
    end
end

if fixed_tol && stop_crit == "Opt"
    results_filename = "res_data_$(problem)_DRS_Opt_Eps.csv"
    results_filepath = joinpath(results_folder, results_filename)
    CSV.write(results_filepath, df)
elseif !fixed_tol && stop_crit == "Opt"
    results_filename = "res_data_$(problem)_DRS_Opt_r0.csv"
    results_filepath = joinpath(results_folder, results_filename)
    CSV.write(results_filepath, df)
elseif fixed_tol && stop_crit == "Fixed_Point"
    results_filename = "res_data_$(problem)_DRS_FP_Eps.csv"
    results_filepath = joinpath(results_folder, results_filename)
    CSV.write(results_filepath, df)
else
    results_filename = "res_data_$(problem)_DRS_FP_r0.csv"
    results_filepath = joinpath(results_folder, results_filename)
    CSV.write(results_filepath, df)
end