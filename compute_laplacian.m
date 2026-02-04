function L = compute_laplacian(A)
    % Calcola matrice Laplaciana da matrice adiacenza
    % Input: A (matrice sparse N×N di adiacenza)
    % Output: L (Laplaciana N×N sparse, L = D - A)
    
    N = size(A, 1);
    
    % Calcola gradi (numero vicini per ogni nodo)
    degrees = full(sum(A, 2));
    
    % Matrice diagonale gradi
    D = spdiags(degrees, 0, N, N);
    
    % Laplaciana
    L = D - A;
    
    fprintf('Laplaciana: %dx%d, nnz=%d (%.4f%% sparse)\n', ...
        N, N, nnz(L), 100*nnz(L)/(N*N));
    fprintf('Grado: min=%.1f, max=%.1f, medio=%.2f\n', ...
        min(degrees), max(degrees), mean(degrees));
end