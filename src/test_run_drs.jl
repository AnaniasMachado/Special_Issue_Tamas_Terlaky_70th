using MAT

include("utility.jl")
include("types.jl")
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