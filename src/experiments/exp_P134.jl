using CSV
using DataFrames
using MAT
using Base.GC
using Statistics

include("../types.jl")
include("../utility.jl")
include("../methods/solvers.jl")
include("../methods/solvers_cal.jl")
include("../methods/misc.jl")
include("../methods/admm.jl")
include("../methods/drs.jl")

methods = ["Gurobi", "Gurobi_Cal", "ADMM", "DRS"]
method = methods[3]

# Mixed parameters
problems = ["P134"]
problem = problems[1]
epsilon = 10^(-5)
eps_abs = epsilon
eps_rel = 10^(-3)
fixed_tol = false
eps_opt = epsilon
time_limit = 7200

# Gurobi parameters
constraints_set = [["P1", "P3", "P4"], ["PMN", "P3"], ["PLS", "PMN"], ["P13R", "P14R"], ["PMX"]]
constraints = constraints_set[4]

# ADMM parameters
rho = 3.0

# DRS parameters
lambda = 10^(-2)

stop_crits = ["Opt", "Fixed_Point"]
stop_crit = stop_crits[1]

matrices_folder = "./instances/square_dense"
m_values = [100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000]

results_folder = "results/problem_$(problem)"

solutions_folder = "./solutions/problem_$(problem)"

df = DataFrame()

min_unsolvable_m = Inf

for m in m_values
    n = m
    r = Int(m / 4)
    mat_file = "A_m$(m)_n$(n)_r$(r)_d100_idx1.mat"

    mat_path = joinpath(matrices_folder, mat_file)
    mat_data = matread(mat_path)

    println("Solving for matrix: $mat_path")
    
    A = mat_data["matrix"]
    A = Matrix(A)
    AMP = pinv(A)

    data = DataInst(A, m, n, r, AMP=AMP)

    norm_0 = -1.0
    norm_1 = -1.0
    rank = -1.0
    time = -1.0
    factor = -1.0
    if method == "Gurobi"
        try
            time = @elapsed begin
                H = gurobi_solver(data, constraints, eps_opt, time_limit)
            end
            norm_0 = matrix_norm_0(H)
            norm_1 = norm(H, 1)
            rank = calculate_rank(H)

            problem_label = join(constraints, "_")

            solution_filename = "Gurobi/problem_$(problem_label)_m_$(m)_n_$(n)"
            solution_filepath = joinpath(solutions_folder, solution_filename)
            matwrite(solution_filepath, Dict("H" => H, "time" => time))
        catch e
            if isa(e, ErrorException)
                global min_unsolvable_m = min(m, min_unsolvable_m)
            else
                throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
            end
        end
    elseif method == "Gurobi_Cal"
        try
            time = @elapsed begin
                H = gurobi_solver_cal(data, problem, eps_opt, time_limit)
            end
            norm_0 = matrix_norm_0(H)
            norm_1 = norm(H, 1)
            rank = calculate_rank(H)

            solution_filename = "Gurobi_Cal/problem_$(problem)_m_$(m)_n_$(n)"
            solution_filepath = joinpath(solutions_folder, solution_filename)
            matwrite(solution_filepath, Dict("H" => H, "time" => time))
        catch e
            if isa(e, ErrorException)
                global min_unsolvable_m = min(m, min_unsolvable_m)
            else
                throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
            end
        end
    elseif method == "ADMM"
        H, k, time, factor = compare_drs_fp(method, A, rho, lambda, eps_abs, eps_rel, problem, fixed_tol, eps_opt, stop_crit, time_limit)
        norm_0 = matrix_norm_0(H)
        norm_1 = norm(H, 1)
        rank = calculate_rank(H)

        if fixed_tol
            solution_filename = "ADMMe/problem_$(problem)_m_$(m)_n_$(n)"
            solution_filepath = joinpath(solutions_folder, solution_filename)
            matwrite(solution_filepath, Dict("H" => H, "time" => time))
        else
            solution_filename = "ADMM/problem_$(problem)_m_$(m)_n_$(n)"
            solution_filepath = joinpath(solutions_folder, solution_filename)
            matwrite(solution_filepath, Dict("H" => H, "time" => time))
        end
    elseif method == "DRS"
        H, k, time, factor = compare_drs_fp(method, A, rho, lambda, eps_abs, eps_rel, problem, fixed_tol, eps_opt, stop_crit, time_limit)
        norm_0 = matrix_norm_0(H)
        norm_1 = norm(H, 1)
        rank = calculate_rank(H)

        if fixed_tol && stop_crit == "Opt"
            solution_filename = "DRS_Opt_Eps/problem_$(problem)_m_$(m)_n_$(n)"
            solution_filepath = joinpath(solutions_folder, solution_filename)
            matwrite(solution_filepath, Dict("H" => H, "time" => time, "k" => k))
        elseif !fixed_tol && stop_crit == "Opt"
            solution_filename = "DRS_Opt_r0/problem_$(problem)_m_$(m)_n_$(n)"
            solution_filepath = joinpath(solutions_folder, solution_filename)
            matwrite(solution_filepath, Dict("H" => H, "time" => time, "k" => k))
        elseif fixed_tol && stop_crit == "Fixed_Point"
            solution_filename = "DRS_FP_Eps/problem_$(problem)_m_$(m)_n_$(n)"
            solution_filepath = joinpath(solutions_folder, solution_filename)
            matwrite(solution_filepath, Dict("H" => H, "time" => time, "k" => k))
        elseif !fixed_tol && stop_crit == "Fixed_Point"
            solution_filename = "DRS_FP_r0/problem_$(problem)_m_$(m)_n_$(n)"
            solution_filepath = joinpath(solutions_folder, solution_filename)
            matwrite(solution_filepath, Dict("H" => H, "time" => time, "k" => k))
        end
    else
        throw(ErrorException("Invalid method chose."))
    end

    GC.gc()

    result = DataFrame(
        m = [m],
        n = [n],
        r = [r],
        AMP_norm_0 = [matrix_norm_0(AMP)],
        AMP_norm_1 = [norm(AMP, 1)],
        norm_0 = [norm_0],
        norm_1 = [norm_1],
        rank = [rank],
        time = [time],
        factor = [factor]
    )

    append!(df, result)

    GC.gc()
end

if method == "Gurobi"
    problem_label = join(constraints, "_")
    results_filename = "results_$(problem_label)_Gurobi.csv"
    results_filepath = joinpath(results_folder, results_filename)
    CSV.write(results_filepath, df)
elseif method == "Gurobi_Cal"
    results_filename = "results_$(problem)_Gurobi_Cal.csv"
    results_filepath = joinpath(results_folder, results_filename)
    CSV.write(results_filepath, df)
elseif method == "ADMM"
    if fixed_tol
        results_filename = "results_$(problem)_ADMMe.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    else
        results_filename = "results_$(problem)_ADMM.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    end
elseif method == "DRS"
    if fixed_tol && stop_crit == "Opt"
        results_filename = "results_$(problem)_DRS_Opt_Eps.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    elseif !fixed_tol && stop_crit == "Opt"
        results_filename = "results_$(problem)_DRS_Opt_r0.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    elseif fixed_tol && stop_crit == "Fixed_Point"
        results_filename = "results_$(problem)_DRS_FP_Eps.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    else
        results_filename = "results_$(problem)_DRS_FP_r0.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    end
else
    throw(ErrorException("Invalid method chose."))
end