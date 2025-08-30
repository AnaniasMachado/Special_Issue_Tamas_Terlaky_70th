using MAT

include("utility.jl")
include("types.jl")
include("./methods/misc.jl")
include("./methods/admm.jl")
include("./methods/drs.jl")

m = 100
n = 50
r = 25

A = gen_random_rank_r_matrix(m, n, r)

matrix_folder = "./instances/rectangular_dense_02"
mat_files = readdir(matrix_folder)
mat_file = mat_files[2]
mat_path = joinpath(matrix_folder, mat_file)
mat_data = matread(mat_path)
A = mat_data["A"]
A = Matrix(A)
m, n = size(A)
r = rank(A)

println("m = $m, n = $n, r = $r")

problem = "P123"
rho = 3.0
lambda = 10^(-2)

epsilon = 10^(-5)
eps_opt = epsilon
eps_abs = epsilon
eps_rel = 10^(-3)
fixed_tol = false

time_limit = 1200
stop_crits = ["Opt", "Fixed_Point"]
stop_crit = stop_crits[2]

DRS_time = @elapsed begin
    DRS_H, DRS_k = drs(A, lambda, eps_abs, eps_rel, problem, fixed_tol, eps_opt, stop_crit, time_limit)
end
DRS_H_norm_0 = matrix_norm_0(DRS_H)
DRS_H_norm_1 = norm(DRS_H, 1)

println("Stop crit: $stop_crit")
println("Fixed tol: $fixed_tol")
println("DRS time: $DRS_time")
println("DRS k: $DRS_k")
println("DRS norm 1: $DRS_H_norm_1")
println("DRS norm 0: $DRS_H_norm_0")
println("DRS bound ratio: $(DRS_H_norm_0 / (m * r))")

stop_crits = ["Opt", "Fixed_Point"]
stop_crit = stop_crits[1]
method = "DRS"

DRS_Opt_H, DRS_Opt_k, DRS_Opt_time = compare_drs_fp(method, A, rho, lambda, eps_abs, eps_rel, problem, fixed_tol, eps_opt, stop_crit, time_limit)
DRS_Opt_H_norm_0 = matrix_norm_0(DRS_Opt_H)
DRS_Opt_H_norm_1 = norm(DRS_Opt_H, 1)

println("----------------")

println("Stop crit: $stop_crit")
println("Fixed tol: $fixed_tol")
println("DRS Opt time: $DRS_Opt_time")
println("DRS Opt k: $DRS_Opt_k")
println("DRS Opt norm 1: $DRS_Opt_H_norm_1")
println("DRS Opt norm 0: $DRS_Opt_H_norm_0")
println("DRS Opt bound ratio: $(DRS_Opt_H_norm_0 / (m * r))")
factor = abs(norm(DRS_Opt_H, 1) - norm(DRS_H, 1)) / norm(DRS_H, 1)
println("DRS Opt factor: $(factor)")

method = "ADMM"

ADMM_H, ADMM_k, ADMM_time = compare_drs_fp(method, A, rho, lambda, eps_abs, eps_rel, problem, fixed_tol, eps_opt, stop_crit, time_limit)
ADMM_H_norm_0 = matrix_norm_0(ADMM_H)
ADMM_H_norm_1 = norm(ADMM_H, 1)

println("----------------")

println("Fixed tol: $fixed_tol")
println("ADMM time: $ADMM_time")
println("ADMM k: $ADMM_k")
println("ADMM norm 1: $ADMM_H_norm_1")
println("ADMM norm 0: $ADMM_H_norm_0")
println("ADMM bound ratio: $(ADMM_H_norm_0 / (m * r))")
factor = abs(norm(ADMM_H, 1) - norm(DRS_H, 1)) / norm(DRS_H, 1)
println("ADMM factor: $(factor)")