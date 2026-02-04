function [Z, stats] = BE_method(L, V, T, n_steps, tol, maxiter, method, animation_callback)
    % Metodo di Eulero all'indietro per equazione del calore
    % Input: L (Laplaciano), V (vertici iniziali N×3), T (tempo finale),
    %        n_steps, tol, maxiter, method ('jacobi'/'gauss-seidel')
    %        animation_callback (opzionale, function handle)
    % Output: Z (cell array soluzioni), stats (statistiche)
    
    % Parametro opzionale
    if nargin < 8
        animation_callback = [];
    end
    
    N = size(V, 1);
    dt = T / n_steps;
    
    % Prealloca output
    Z = cell(n_steps + 1, 1);
    Z{1} = V;
    
    % Statistiche
    stats.iterations = zeros(n_steps, 1);
    stats.residuals = zeros(n_steps, 1);
    stats.times = zeros(n_steps, 1);
    
    % Matrice sistema (costante per tutti i passi)
    M = speye(N) + dt * L;
    
    % Callback iniziale
    if ~isempty(animation_callback)
        animation_callback(V, 0, 0);
    end
    
    % Time stepping
    for j = 1:n_steps
        t_step = dt * j;
        
        tic;
        % Risolvi 3 sistemi (uno per coordinata x,y,z)
        V_new = zeros(N, 3);
        total_iter = 0;
        max_res = 0;
        
        for k = 1:3
            [V_new(:,k), iter_k, res_k] = MyMatrixSolver(M, Z{j}(:,k), tol, maxiter, method);
            total_iter = total_iter + iter_k;
            max_res = max(max_res, res_k(end));
        end
        
        stats.times(j) = toc;
        stats.iterations(j) = total_iter / 3;
        stats.residuals(j) = max_res;
        
        Z{j+1} = V_new;
        
        % Callback animazione
        if ~isempty(animation_callback)
            animation_callback(V_new, t_step, j);
        end
    end
    
    % Statistiche finali
    stats.final_residual = stats.residuals(end);
    stats.total_time = sum(stats.times);
    stats.avg_iterations = mean(stats.iterations);
end