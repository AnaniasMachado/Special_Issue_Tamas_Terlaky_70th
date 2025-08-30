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

methods = ["ADMM", "DRS"]
method = methods[1]

problem = "P123"
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

if method == "ADMM"
    H, k, res_data, sol_data, time_data = admm_p123_run_data(A, rho, eps_abs, eps_rel, fixed_tol, eps_opt, time_limit)

    for i in 1:length(res_data)
        result = DataFrame(
            pri_res = [res_data[i][1]],
            dual_res = [res_data[i][2]],
            H_norm_0 = [sol_data[i][1]],
            H_norm_1 = [sol_data[i][2]],
            time = [time_data[i]]
        )
        append!(df, result)
    end
elseif method == "DRS"
    H, k, res_data, sol_data, time_data = drs_run_data(A, lambda, eps_abs, eps_rel, problem, fixed_tol, eps_opt, stop_crit, time_limit)

    for i in 1:length(res_data)
        if stop_crit == "Opt"
            result = DataFrame(
                pri_res = [res_data[i][1]],
                dual_res = [res_data[i][2]],
                H_norm_0 = [sol_data[i][1]],
                H_norm_1 = [sol_data[i][2]],
                time = [time_data[i]]
            )
            append!(df, result)
        elseif stop_crit == "Fixed_Point"
            result = DataFrame(
                res = [res_data[i]],
                H_norm_0 = [sol_data[i][1]],
                H_norm_1 = [sol_data[i][2]],
                time = [time_data[i]]
            )
            append!(df, result)
        end
    end
end

if method == "ADMM"
    results_filename = "run_data_$(problem)_ADMM.csv"
    results_filepath = joinpath(results_folder, results_filename)
    CSV.write(results_filepath, df)
elseif method == "DRS"
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