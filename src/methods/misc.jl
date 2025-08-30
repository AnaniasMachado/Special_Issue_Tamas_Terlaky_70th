
function compare_drs_fp(method::String, A::Matrix{Float64}, rho::Float64, lambda::Float64, eps_abs::Float64, eps_rel::Float64, problem::String, fixed_tol::Bool, eps_opt::Float64, stop_crit::String, time_limit::Int64)
    time_ref = @elapsed begin
        H_ref, k_ref = drs(A, lambda, eps_abs, eps_rel, problem, fixed_tol, eps_opt, "Fixed_Point", time_limit)
    end
    if (method == "DRS") && (stop_crit == "Fixed_Point")
        return H_ref, k_ref, time_ref, 0.0
    elseif method == "DRS"
        H, max_iter = drs_closeness(H_ref, A, lambda, eps_abs, eps_rel, problem, fixed_tol, eps_opt, stop_crit, time_limit)
        time = @elapsed begin
            H, k = drs_max_iter(max_iter, A, lambda, eps_abs, eps_rel, problem, fixed_tol, eps_opt, stop_crit, time_limit)
        end
        factor = abs(norm(H, 1) - norm(H_ref, 1)) / norm(H_ref, 1)
        return H, k, time, factor
    elseif (method == "ADMM") && (problem == "P123")
        H, max_iter = admm_p123_closeness(H_ref, A, rho, eps_abs, eps_rel, fixed_tol, eps_opt, time_limit)
        time = @elapsed begin
            H, k = admm_p123_max_iter(max_iter, A, rho, eps_abs, eps_rel, fixed_tol, eps_opt, time_limit)
        end
        factor = abs(norm(H, 1) - norm(H_ref, 1)) / norm(H_ref, 1)
        return H, k, time, factor
    elseif (method == "ADMM") && (problem == "P134")
        H, max_iter = admm_p134_closeness(H_ref, A, rho, eps_abs, eps_rel, fixed_tol, eps_opt, time_limit)
        time = @elapsed begin
            H, k = admm_p134_max_iter(max_iter, A, rho, eps_abs, eps_rel, fixed_tol, eps_opt, time_limit)
        end
        factor = abs(norm(H, 1) - norm(H_ref, 1)) / norm(H_ref, 1)
        return H, k, time, factor
    end
end