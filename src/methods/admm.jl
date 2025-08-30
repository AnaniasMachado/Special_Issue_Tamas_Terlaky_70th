using LinearAlgebra

epsilon = 10^(-5)

function soft_thresholding_matrix(X::Matrix{Float64}, lambda::Float64)
    return sign.(X) .* max.(abs.(X) .- lambda, 0)
end

function variables_initialization(V1::Matrix{Float64}, U1::Matrix{Float64}, D_inv::Matrix{Float64}, rho::Float64)
    Theta = (V1 * U1') / norm(V1 * U1', Inf)
    Lambda = Theta / rho
    E = V1 * D_inv * U1' + Lambda
    return Lambda, E
end

function admm_p123(A::Matrix{Float64}, rho::Float64, eps_abs::Float64, eps_rel::Float64, fixed_tol::Bool, eps_opt::Float64, time_limit::Int64)
    start_time = time()
    m, n = size(A)
    U, S, V = svd(A, full=true)
    S = Diagonal(S)
    r = count_singular_values(S)
    D = S[1:r, 1:r]
    D_inv = inv(D)
    U1 = U[:, 1:r]
    V1 = V[:, 1:r]
    V2 = V[:, r+1:end]

    V1DinvU1T = V1 * D_inv * U1'
    V1DinvU1T_F = norm(V1DinvU1T)

    V2V2T = V2 * V2'
    U1U1T = U1 * U1'

    Z = zeros(n - r, r)
    H = zeros(n, m)
    Lambda, Ekm = variables_initialization(V1, U1, D_inv, rho)
    k = 0
    while true
        k += 1
        V2ZU1T = V2V2T * (Ekm - Lambda) * U1U1T
        H = V1DinvU1T + V2ZU1T
        # Updates Ek
        Ek = soft_thresholding_matrix(H + Lambda, 1/rho)
        res_infeas = H - Ek
        Lambda += res_infeas
        # Calculates stop criterion variables
        rk_F = norm(res_infeas)
        sk_F = rho * norm(V2' * (Ek - Ekm) * U1)
        if !fixed_tol
            matrix_norms = [norm(Ek), norm(V2ZU1T), V1DinvU1T_F]
            primal_upper_bound = eps_abs * sqrt(m*n) + eps_rel * maximum(matrix_norms)
            dual_upper_bound = eps_abs * sqrt((n-r)*r) + eps_rel * rho * norm(V2' * Lambda * U1)
            if (rk_F <= primal_upper_bound) && (sk_F <= dual_upper_bound)
                break
            end
        else
            if (rk_F <= eps_opt) && (sk_F <= eps_opt)
                break
            end
        end
        # Makes Ek the new Ek-1
        Ekm = Ek
        # Checks time limit
        elapsed_time = time() - start_time
        if elapsed_time > time_limit
            println("TimeLimit: ADMM exceed time limit to solve the problem.")
            return "-"
        end
    end
    return H
end

function admm_p134(A::Matrix{Float64}, rho::Float64, eps_abs::Float64, eps_rel::Float64, fixed_tol::Bool, eps_opt::Float64, time_limit::Int64)
    start_time = time()
    m, n = size(A)
    U, S, V = svd(A, full=true)
    S = Diagonal(S)
    r = count_singular_values(S)
    D = S[1:r, 1:r]
    D_inv = inv(D)
    U1 = U[:, 1:r]
    U2 = U[:, r+1:end]
    V1 = V[:, 1:r]
    V2 = V[:, r+1:end]

    V1DinvU1T = V1 * D_inv * U1'
    V1DinvU1T_F = norm(V1DinvU1T)

    V2V2T = V2 * V2'
    U1U1T = U1 * U1'
    U2U2T = U2 * U2'

    Z = zeros(n - r, r)
    H = zeros(n, m)
    Lambda, Ekm = variables_initialization(V1, U1, D_inv, rho)
    k = 0
    while true
        k += 1
        V2WU2T = V2V2T * (Ekm - Lambda) * U2U2T
        H = V1DinvU1T + V2WU2T
        # Updates Ek
        Ek = soft_thresholding_matrix(H + Lambda, 1/rho)
        res_infeas = H - Ek
        Lambda += res_infeas
        # Calculates stop criterion variables
        rk_F = norm(res_infeas)
        sk_F = rho * norm(V2' * (Ek - Ekm) * U2)
        if !fixed_tol
            matrix_norms = [norm(Ek), norm(V2WU2T), V1DinvU1T_F]
            primal_upper_bound = eps_abs * sqrt(m*n) + eps_rel * maximum(matrix_norms)
            dual_upper_bound = eps_abs * sqrt((n-r)*r) + eps_rel * rho * norm(V2' * Lambda * U2)
            if (rk_F <= primal_upper_bound) && (sk_F <= dual_upper_bound)
                break
            end
        else
            if (rk_F <= eps_opt) && (sk_F <= eps_opt)
                break
            end
        end
        # Makes Ek the new Ek-1
        Ekm = Ek
        # Checks time limit
        elapsed_time = time() - start_time
        if elapsed_time > time_limit
            println("TimeLimit: ADMM exceed time limit to solve the problem.")
            return "-"
        end
    end
    return H
end

function admm_p123_run_data(A::Matrix{Float64}, rho::Float64, eps_abs::Float64, eps_rel::Float64, fixed_tol::Bool, eps_opt::Float64, time_limit::Int64)
    start_time = time()
    m, n = size(A)
    U, S, V = svd(A, full=true)
    S = Diagonal(S)
    r = count_singular_values(S)
    D = S[1:r, 1:r]
    D_inv = inv(D)
    U1 = U[:, 1:r]
    V1 = V[:, 1:r]
    V2 = V[:, r+1:end]

    V1DinvU1T = V1 * D_inv * U1'
    V1DinvU1T_F = norm(V1DinvU1T)

    V2V2T = V2 * V2'
    U1U1T = U1 * U1'

    Z = zeros(n - r, r)
    H = zeros(n, m)
    Lambda, Ekm = variables_initialization(V1, U1, D_inv, rho)
    res_data = []
    sol_data = []
    time_data = []
    k = 0
    while true
        k += 1
        V2ZU1T = V2V2T * (Ekm - Lambda) * U1U1T
        H = V1DinvU1T + V2ZU1T
        # Updates Ek
        Ek = soft_thresholding_matrix(H + Lambda, 1/rho)
        res_infeas = H - Ek
        Lambda += res_infeas
        # Calculates stop criterion variables
        rk_F = norm(res_infeas)
        sk_F = rho * norm(V2' * (Ek - Ekm) * U1)
        # Makes Ek the new Ek-1
        Ekm = Ek
        # Checks time limit
        elapsed_time = time() - start_time
        push!(res_data, [rk_F, sk_F])
        push!(sol_data, [matrix_norm_0(H), norm(H, 1)])
        push!(time_data, elapsed_time)
        if elapsed_time > time_limit
            println("TimeLimit: ADMM exceed time limit to solve the problem.")
            return H, k, res_data, sol_data, time_data
        end
    end
    return H, k, res_data, sol_data, time_data
end

function admm_p134_run_data(A::Matrix{Float64}, rho::Float64, eps_abs::Float64, eps_rel::Float64, fixed_tol::Bool, eps_opt::Float64, time_limit::Int64)
    start_time = time()
    m, n = size(A)
    U, S, V = svd(A, full=true)
    S = Diagonal(S)
    r = count_singular_values(S)
    D = S[1:r, 1:r]
    D_inv = inv(D)
    U1 = U[:, 1:r]
    U2 = U[:, r+1:end]
    V1 = V[:, 1:r]
    V2 = V[:, r+1:end]

    V1DinvU1T = V1 * D_inv * U1'
    V1DinvU1T_F = norm(V1DinvU1T)

    V2V2T = V2 * V2'
    U1U1T = U1 * U1'
    U2U2T = U2 * U2'

    Z = zeros(n - r, r)
    H = zeros(n, m)
    Lambda, Ekm = variables_initialization(V1, U1, D_inv, rho)
    res_data = []
    sol_data = []
    time_data = []
    k = 0
    while true
        k += 1
        V2WU2T = V2V2T * (Ekm - Lambda) * U2U2T
        H = V1DinvU1T + V2WU2T
        # Updates Ek
        Ek = soft_thresholding_matrix(H + Lambda, 1/rho)
        res_infeas = H - Ek
        Lambda += res_infeas
        # Calculates stop criterion variables
        rk_F = norm(res_infeas)
        sk_F = rho * norm(V2' * (Ek - Ekm) * U2)
        # Makes Ek the new Ek-1
        Ekm = Ek
        # Checks time limit
        elapsed_time = time() - start_time
        push!(res_data, [rk_F, sk_F])
        push!(sol_data, [matrix_norm_0(H), norm(H, 1)])
        push!(time_data, elapsed_time)
        if elapsed_time > time_limit
            println("TimeLimit: ADMM exceed time limit to solve the problem.")
            return H, k, res_data, sol_data, time_data
        end
    end
    return H, k, res_data, sol_data, time_data
end

function admm_p123_closeness(H_ref::Matrix{Float64}, A::Matrix{Float64}, rho::Float64, eps_abs::Float64, eps_rel::Float64, fixed_tol::Bool, eps_opt::Float64, time_limit::Int64)
    start_time = time()
    m, n = size(A)
    U, S, V = svd(A, full=true)
    S = Diagonal(S)
    r = count_singular_values(S)
    D = S[1:r, 1:r]
    D_inv = inv(D)
    U1 = U[:, 1:r]
    V1 = V[:, 1:r]
    V2 = V[:, r+1:end]

    V1DinvU1T = V1 * D_inv * U1'
    V1DinvU1T_F = norm(V1DinvU1T)

    V2V2T = V2 * V2'
    U1U1T = U1 * U1'

    Z = zeros(n - r, r)
    H = zeros(n, m)
    Lambda, Ekm = variables_initialization(V1, U1, D_inv, rho)
    k = 0
    while true
        k += 1
        V2ZU1T = V2V2T * (Ekm - Lambda) * U1U1T
        H = V1DinvU1T + V2ZU1T
        # Updates Ek
        Ek = soft_thresholding_matrix(H + Lambda, 1/rho)
        res_infeas = H - Ek
        Lambda += res_infeas
        factor = abs(norm(H, 1) - norm(H_ref, 1)) / norm(H_ref, 1)
        if (norm(H, 1) < norm(H_ref, 1)) || (factor <= 10^(-5))
            break
        end
        # Makes Ek the new Ek-1
        Ekm = Ek
        # Checks time limit
        elapsed_time = time() - start_time
        if elapsed_time > time_limit
            println("TimeLimit: ADMM exceed time limit to solve the problem.")
            return H, k
        end
    end
    return H, k
end

function admm_p123_max_iter(max_iter::Int64, A::Matrix{Float64}, rho::Float64, eps_abs::Float64, eps_rel::Float64, fixed_tol::Bool, eps_opt::Float64, time_limit::Int64)
    start_time = time()
    m, n = size(A)
    U, S, V = svd(A, full=true)
    S = Diagonal(S)
    r = count_singular_values(S)
    D = S[1:r, 1:r]
    D_inv = inv(D)
    U1 = U[:, 1:r]
    V1 = V[:, 1:r]
    V2 = V[:, r+1:end]

    V1DinvU1T = V1 * D_inv * U1'
    V1DinvU1T_F = norm(V1DinvU1T)

    V2V2T = V2 * V2'
    U1U1T = U1 * U1'

    Z = zeros(n - r, r)
    H = zeros(n, m)
    Lambda, Ekm = variables_initialization(V1, U1, D_inv, rho)
    k = 0
    while (k <= max_iter)
        k += 1
        V2ZU1T = V2V2T * (Ekm - Lambda) * U1U1T
        H = V1DinvU1T + V2ZU1T
        # Updates Ek
        Ek = soft_thresholding_matrix(H + Lambda, 1/rho)
        res_infeas = H - Ek
        Lambda += res_infeas
        # Calculates stop criterion variables
        rk_F = norm(res_infeas)
        sk_F = rho * norm(V2' * (Ek - Ekm) * U1)
        if !fixed_tol
            matrix_norms = [norm(Ek), norm(V2ZU1T), V1DinvU1T_F]
            primal_upper_bound = eps_abs * sqrt(m*n) + eps_rel * maximum(matrix_norms)
            dual_upper_bound = eps_abs * sqrt((n-r)*r) + eps_rel * rho * norm(V2' * Lambda * U1)
            if (rk_F <= primal_upper_bound) && (sk_F <= dual_upper_bound)
                continue
            end
        else
            if (rk_F <= eps_opt) && (sk_F <= eps_opt)
                continue
            end
        end
        # Makes Ek the new Ek-1
        Ekm = Ek
        # Checks time limit
        elapsed_time = time() - start_time
        if elapsed_time > time_limit
            println("TimeLimit: ADMM exceed time limit to solve the problem.")
            return H, k
        end
    end
    return H, k
end

function admm_p134_closeness(H_ref::Matrix{Float64}, A::Matrix{Float64}, rho::Float64, eps_abs::Float64, eps_rel::Float64, fixed_tol::Bool, eps_opt::Float64, time_limit::Int64)
    start_time = time()
    m, n = size(A)
    U, S, V = svd(A, full=true)
    S = Diagonal(S)
    r = count_singular_values(S)
    D = S[1:r, 1:r]
    D_inv = inv(D)
    U1 = U[:, 1:r]
    U2 = U[:, r+1:end]
    V1 = V[:, 1:r]
    V2 = V[:, r+1:end]

    V1DinvU1T = V1 * D_inv * U1'
    V1DinvU1T_F = norm(V1DinvU1T)

    V2V2T = V2 * V2'
    U1U1T = U1 * U1'
    U2U2T = U2 * U2'

    Z = zeros(n - r, r)
    H = zeros(n, m)
    Lambda, Ekm = variables_initialization(V1, U1, D_inv, rho)
    k = 0
    while true
        k += 1
        V2WU2T = V2V2T * (Ekm - Lambda) * U2U2T
        H = V1DinvU1T + V2WU2T
        # Updates Ek
        Ek = soft_thresholding_matrix(H + Lambda, 1/rho)
        res_infeas = H - Ek
        Lambda += res_infeas
        factor = abs(norm(H, 1) - norm(H_ref, 1)) / norm(H_ref, 1)
        if (norm(H, 1) < norm(H_ref, 1)) || (factor <= 10^(-5))
            break
        end
        # Makes Ek the new Ek-1
        Ekm = Ek
        # Checks time limit
        elapsed_time = time() - start_time
        if elapsed_time > time_limit
            println("TimeLimit: ADMM exceed time limit to solve the problem.")
            return H, k
        end
    end
    return H, k
end

function admm_p134_max_iter(max_iter::Int64, A::Matrix{Float64}, rho::Float64, eps_abs::Float64, eps_rel::Float64, fixed_tol::Bool, eps_opt::Float64, time_limit::Int64)
    start_time = time()
    m, n = size(A)
    U, S, V = svd(A, full=true)
    S = Diagonal(S)
    r = count_singular_values(S)
    D = S[1:r, 1:r]
    D_inv = inv(D)
    U1 = U[:, 1:r]
    U2 = U[:, r+1:end]
    V1 = V[:, 1:r]
    V2 = V[:, r+1:end]

    V1DinvU1T = V1 * D_inv * U1'
    V1DinvU1T_F = norm(V1DinvU1T)

    V2V2T = V2 * V2'
    U1U1T = U1 * U1'
    U2U2T = U2 * U2'

    Z = zeros(n - r, r)
    H = zeros(n, m)
    Lambda, Ekm = variables_initialization(V1, U1, D_inv, rho)
    k = 0
    while (k <= max_iter)
        k += 1
        V2WU2T = V2V2T * (Ekm - Lambda) * U2U2T
        H = V1DinvU1T + V2WU2T
        # Updates Ek
        Ek = soft_thresholding_matrix(H + Lambda, 1/rho)
        res_infeas = H - Ek
        Lambda += res_infeas
        # Calculates stop criterion variables
        rk_F = norm(res_infeas)
        sk_F = rho * norm(V2' * (Ek - Ekm) * U2)
        if !fixed_tol
            matrix_norms = [norm(Ek), norm(V2WU2T), V1DinvU1T_F]
            primal_upper_bound = eps_abs * sqrt(m*n) + eps_rel * maximum(matrix_norms)
            dual_upper_bound = eps_abs * sqrt((n-r)*r) + eps_rel * rho * norm(V2' * Lambda * U2)
            if (rk_F <= primal_upper_bound) && (sk_F <= dual_upper_bound)
                continue
            end
        else
            if (rk_F <= eps_opt) && (sk_F <= eps_opt)
                continue
            end
        end
        # Makes Ek the new Ek-1
        Ekm = Ek
        # Checks time limit
        elapsed_time = time() - start_time
        if elapsed_time > time_limit
            println("TimeLimit: ADMM exceed time limit to solve the problem.")
            return H, k
        end
    end
    return H, k
end