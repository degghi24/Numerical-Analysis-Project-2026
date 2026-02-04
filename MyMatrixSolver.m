function [X, iter, res] = MyMatrixSolver(M, b, tol, maxiter, method)
    % Risolve sistema lineare M*X = b con metodo iterativo
    % Input: M (matrice), b (termine noto), tol, maxiter, method ('jacobi'/'gauss-seidel')
    % Output: X (soluzione), iter (num iterazioni), res (vettore residui)
    
    [n, k] = size(b);
    
    % Estrai componenti matrice (formato sparse)
    D = spdiags(diag(M), 0, n, n);
    d_vec = diag(M);
    
    % Verifica diagonale non nulla
    if any(abs(d_vec) < eps)
        error('Matrice ha zeri sulla diagonale');
    end
    
    D_inv = spdiags(1 ./ d_vec, 0, n, n);
    L = tril(M, -1);
    U = triu(M, 1);
    
    X = zeros(n, k);
    res = zeros(maxiter, 1);
    
    switch lower(method)
        case 'jacobi'
            for iter = 1:maxiter
                X_new = D_inv * (b - (L + U) * X);
                err = norm(X_new - X, 'fro') / (norm(X_new, 'fro') + eps);
                res(iter) = err;
                
                if err < tol
                    X = X_new;
                    res = res(1:iter);
                    return;
                end
                
                X = X_new;
            end
            
        case 'gauss-seidel'
            D_plus_L = D + L;
            
            for iter = 1:maxiter
                X_new = D_plus_L \ (b - U * X);
                err = norm(X_new - X, 'fro') / (norm(X_new, 'fro') + eps);
                res(iter) = err;
                
                if err < tol
                    X = X_new;
                    res = res(1:iter);
                    return;
                end
                
                X = X_new;
            end
            
        otherwise
            error('Metodo non riconosciuto: %s', method);
    end
    
    warning('Non convergenza dopo %d iterazioni (residuo: %.2e)', maxiter, res(end));
end