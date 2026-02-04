function A = build_adjacency(F)
    % Costruisce matrice di adiacenza del grafo mesh
    % Input: F (matrice M×3 facce, ogni riga = triangolo [v1, v2, v3])
    % Output: A (matrice sparse N×N, A(i,j)=1 se i~j)
    
    N = max(F(:));
    num_faces = size(F, 1);
    num_edges_total = 3 * num_faces;
    
    % Pre-alloca liste edge
    edges_i = zeros(num_edges_total, 1);
    edges_j = zeros(num_edges_total, 1);
    
    % Estrai edge da tutte le facce
    for k = 1:num_faces
        base_idx = 3 * (k - 1);
        
        % Edge del triangolo: (v1,v2), (v2,v3), (v3,v1)
        edges_i(base_idx + 1) = F(k, 1);
        edges_j(base_idx + 1) = F(k, 2);
        
        edges_i(base_idx + 2) = F(k, 2);
        edges_j(base_idx + 2) = F(k, 3);
        
        edges_i(base_idx + 3) = F(k, 3);
        edges_j(base_idx + 3) = F(k, 1);
    end
    
    % Costruisci matrice sparse simmetrica
    all_i = [edges_i; edges_j];
    all_j = [edges_j; edges_i];
    
    A = sparse(all_i, all_j, ones(length(all_i), 1), N, N);
    A = double(A > 0);  % Binaria
    A = A - spdiags(diag(A), 0, N, N);  % Rimuovi diagonale
    
    % Info
    num_edges = nnz(A) / 2;
    sparsity = 100 * nnz(A) / (N * N);
    
    fprintf('Grafo: %d nodi, %d edge\n', N, num_edges);
    fprintf('Sparsita: %.4f%% (%.2e elementi non-nulli)\n', sparsity, nnz(A));
    fprintf('Grado medio: %.2f\n', 2 * num_edges / N);
end