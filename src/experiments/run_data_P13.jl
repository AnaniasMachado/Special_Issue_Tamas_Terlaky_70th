using CSV
using DataFrames
using MAT

include("../types.jl")
include("../utility.jl")
include("../methods/admm.jl")
include("../methods/drs.jl")

m = 5000
n = Int(m / 2)
r = Int(m / 4)

matrix_folder = "./instances/rectangular_dense"
mat_file = "A_$(m)_$(n)_$(r).mat"
mat_path = joinpath(matrix_folder, mat_file)
mat_data = matread(mat_path)
A = mat_data["A"]
A = Matrix(A)

println("m = $m, n = $n, r = $r")

methods = ["DRS"]
method = methods[2]

problem = "P13"
rho = 3.0
lambda = 10^(-2)

epsilon = 10^(-5)
eps_opt = epsilon
eps_abs = epsilon
eps_rel = 10^(-3)
fixed_tol = false

time_limit = 7200
stop_crits = ["Opt", "Fixed_Point"]
stop_crit = stop_crits[1]

results_folder = "results/problem_$problem"

df = DataFrame()

if method == "DRS"
    DRS_H, DRS_k, DRS_res_data, DRS_sol_data, DRS_time_data = drs_run_data(A, lambda, eps_abs, eps_rel, problem, fixed_tol, eps_opt, stop_crit, time_limit)

    for i in 1:length(DRS_res_data)
        if stop_crit == "Opt"
            result = DataFrame(
                pri_res = [DRS_res_data[i][1]],
                dual_res = [DRS_res_data[i][2]],
                H_norm_0 = [DRS_sol_data[i][1]],
                H_norm_1 = [DRS_sol_data[i][2]],
                time = [DRS_time_data[i]]
            )
            append!(df, result)
        elseif stop_crit == "Fixed_Point"
            result = DataFrame(
                res = [DRS_res_data[i]],
                H_norm_0 = [DRS_sol_data[i][1]],
                H_norm_1 = [DRS_sol_data[i][2]],
                time = [DRS_time_data[i]]
            )
            append!(df, result)
        end
    end
end

if method == "DRS"
    if stop_crit == "Opt"
        results_filename = "run_data_$(problem)_DRS_Opt.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    elseif stop_crit == "Fixed_Point"
        results_filename = "run_data_$(problem)_DRS_FP.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    end
end