m_values = [100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000];
d = 1.0;
n_mtx = 1;
output_dir = '../instances/square_dense_02';

for m = m_values
    n = m;
    r = m / 4;

    generate_experiment_matrices(m, n, r, d, n_mtx, output_dir);
end