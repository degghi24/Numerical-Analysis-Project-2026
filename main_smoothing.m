%% MAIN - Laplacian Smoothing 
% Progetto 8 - Calcolo Numerico
% 
% Autori: Filippo De Guio, Giulio Bottacin
% Data: Febbraio 2025

clear; close all; clc;

fprintf('====================================\n');
fprintf('LAPLACIAN SMOOTHING - Progetto 8\n');
fprintf('Autori: Filippo De Guio, Giulio Bottacin\n');
fprintf('====================================\n\n');

%% CONFIGURAZIONE

% Cartella risultati
results_folder = 'results';
if ~exist(results_folder, 'dir')
    mkdir(results_folder);
end

% Parametri output
gif_framerate = 15;
image_resolution = 300;
animation_skip = 10;  % Salva ogni N frame

%% PARAMETRI MESH

mesh_files = {
    'data/noisybunny.off';
    'data/maxplank.off';
    'data/noisyknot.off'
};

mesh_names = {
    'Noisy Bunny';
    'Max Planck Head';
    'Noisy Knot'
};

mesh_short_names = {
    'bunny';
    'planck';
    'knot'
};

% Viste per visualizzazione
mesh_views = {
    [0, 90];      % bunny
    [160, -45];   % planck
    [-10, 0]      % knot
};

% Parametri smoothing per mesh
T_values = [0.4; 1.0; 0.4];

%% PARAMETRI METODO

n_steps = 200;
tol = 1e-6;           
maxiter = 1000;       
method = 'gauss-seidel';

%% ELABORAZIONE

n_meshes = length(mesh_files);
results = cell(n_meshes, 1);

for idx = 1:n_meshes
    fprintf('\n[%d/%d] Elaborazione: %s\n', idx, n_meshes, mesh_names{idx});
    fprintf('------------------------------------\n');
    
    %% Caricamento mesh
    [V_original, F] = read_off(mesh_files{idx});
    N = size(V_original, 1);
    M = size(F, 1);
    
    % Calcola bounding box
    bbox = max(V_original) - min(V_original);
    bbox_size = norm(bbox);
    fprintf('Bounding box: [%.4f, %.4f, %.4f], size: %.4f\n', ...
        bbox(1), bbox(2), bbox(3), bbox_size);

    %% Costruzione grafo e Laplaciano
    A = build_adjacency(F);
    L = compute_laplacian(A);
    
    %% Verifica condizione convergenza
    max_degree = full(max(sum(A, 2)));
    dt = T_values(idx) / n_steps;
    dt_critico = 1 / max_degree;
    
    fprintf('\n>> VERIFICA CONVERGENZA:\n');
    fprintf('   dt         = %.6f\n', dt);
    fprintf('   1/d_max    = %.6f\n', dt_critico);
    
    if dt >= dt_critico
        warning('   ATTENZIONE: dt >= 1/d_max, convergenza non garantita!');
        fprintf('   Suggerimento: aumentare n_steps a %d o piu\n', ...
                ceil(T_values(idx) * max_degree) + 1);
    else
        fprintf('   OK: Condizione dt < 1/d_max soddisfatta\n');
    end
    fprintf('\n');
    
    %% Setup figura animazione
    gif_filename = fullfile(results_folder, ...
                            sprintf('%s_animation.gif', mesh_short_names{idx}));
    
    fig_anim = figure('Position', [100 100 900 700]);
    fig_anim.Name = sprintf('Animazione: %s', mesh_names{idx});
    
    h_patch = patch('Vertices', V_original, 'Faces', F, ...
                    'FaceColor', [1.0, 0.7, 0.7], ...
                    'EdgeColor', 'none', ...
                    'FaceLighting', 'gouraud');
    
    % Lighting: dal basso per Max Planck (volto verso il basso)
    if idx == 2  % Max Planck Head
        light('Position', [0 0 -1], 'Style', 'infinite');
        light('Position', [1 1 0.5], 'Style', 'infinite', 'Color', [0.3 0.3 0.3]);
    else
        camlight('headlight');
    end
    lighting gouraud;
    material dull;
    
    view(mesh_views{idx});
    axis equal; axis vis3d; axis off;
    
    h_title = title(sprintf('%s - t = 0.0000 / %.4f', ...
                           mesh_names{idx}, T_values(idx)), ...
                   'FontSize', 13, 'FontWeight', 'bold');
    
    %% Callback animazione (semplificato)
    % Usa global per evitare troppi parametri
    global ANIM_DATA;
    ANIM_DATA = struct();
    ANIM_DATA.h_patch = h_patch;
    ANIM_DATA.h_title = h_title;
    ANIM_DATA.mesh_name = mesh_names{idx};
    ANIM_DATA.T_final = T_values(idx);
    ANIM_DATA.gif_file = gif_filename;
    ANIM_DATA.gif_frame_count = 0;
    ANIM_DATA.framerate = gif_framerate;
    ANIM_DATA.skip = animation_skip;
    
    animation_callback = @update_animation;
    
    %% Esegui smoothing con BE_method
    tic_smooth = tic;
    fprintf('>> Esecuzione BE_method...\n');
    
    [Z, ~] = BE_method(L, V_original, T_values(idx), n_steps, ...
                      tol, maxiter, method, animation_callback);
    
    time_smoothing = toc(tic_smooth);
    V_smoothed = Z{end};
    
    fprintf('   Smoothing completato in %.2f secondi\n', time_smoothing);
    
    %% Calcola metriche
    % Spostamento medio
    displacements = sqrt(sum((V_smoothed - V_original).^2, 2));
    mean_displacement = mean(displacements) / bbox_size * 100;
    
    % Energia Laplaciana (riduzione rumore)
    E_init = norm(L * V_original, 'fro')^2 / N;
    E_final = norm(L * V_smoothed, 'fro')^2 / N;
    noise_reduction = 100 * (E_init - E_final) / E_init;
    
    fprintf('   Energia Laplaciana: %.2e -> %.2e\n', E_init, E_final);
    fprintf('   Riduzione rumore: %.1f%%\n', noise_reduction);
    
    % Salva risultati
    results{idx} = struct(...
        'name', mesh_names{idx}, ...
        'N', N, ...
        'M', M, ...
        'T', T_values(idx), ...
        'bbox_size', bbox_size, ...
        'mean_displacement', mean_displacement, ...
        'noise_reduction', noise_reduction, ...
        'time_smoothing', time_smoothing ...
    );
    
    %% Confronto visivo
    fig_comp = figure('Position', [50 50 1400 600]);
    fig_comp.Name = sprintf('Confronto: %s', mesh_names{idx});
    
    % Mesh iniziale
    subplot(1, 2, 1);
    patch('Vertices', V_original, 'Faces', F, ...
          'FaceColor', [1.0, 0.7, 0.7], ...
          'EdgeColor', 'none', ...
          'FaceLighting', 'gouraud');
    if idx == 2
        light('Position', [0 0 -1], 'Style', 'infinite');
        light('Position', [1 1 0.5], 'Style', 'infinite', 'Color', [0.3 0.3 0.3]);
    else
        camlight('headlight');
    end
    lighting gouraud;
    material dull;
    view(mesh_views{idx});
    axis equal; axis vis3d; axis off;
    title('Mesh Iniziale (Rumorosa)', 'FontSize', 13, 'FontWeight', 'bold');
    
    % Mesh smoothed
    subplot(1, 2, 2);
    patch('Vertices', V_smoothed, 'Faces', F, ...
          'FaceColor', [0.6, 0.85, 1.0], ...
          'EdgeColor', 'none', ...
          'FaceLighting', 'gouraud');
    if idx == 2
        light('Position', [0 0 -1], 'Style', 'infinite');
        light('Position', [1 1 0.5], 'Style', 'infinite', 'Color', [0.3 0.3 0.3]);
    else
        camlight('headlight');
    end
    lighting gouraud;
    material dull;
    view(mesh_views{idx});
    axis equal; axis vis3d; axis off;
    title(sprintf('Mesh Smoothed (T=%.2f)', T_values(idx)), ...
          'FontSize', 13, 'FontWeight', 'bold');
    
    % Salva confronto
    filename_comp = fullfile(results_folder, ...
                             sprintf('%s_comparison.png', mesh_short_names{idx}));
    print(fig_comp, filename_comp, '-dpng', sprintf('-r%d', image_resolution));
    
    fprintf('   Salvato: %s\n', filename_comp);
end

%% RESOCONTO FINALE

fprintf('\n');
fprintf('========================================\n');
fprintf('     RESOCONTO FINALE ELABORAZIONE     \n');
fprintf('========================================\n\n');

for idx = 1:n_meshes
    r = results{idx};
    
    fprintf('MESH %d: %s\n', idx, r.name);
    fprintf('----------------------------------------\n');
    fprintf('Dimensioni:\n');
    fprintf('  Vertici:             %d\n', r.N);
    fprintf('  Facce:               %d\n', r.M);
    fprintf('\n');
    fprintf('Parametri Smoothing:\n');
    fprintf('  Tempo finale (T):    %.4f\n', r.T);
    fprintf('  Passi temporali:     %d\n', n_steps);
    fprintf('  Metodo:              %s\n', method);
    fprintf('  Tolleranza:          %.1e\n', tol);
    fprintf('\n');
    fprintf('Risultati:\n');
    fprintf('  Spostamento medio:   %.2f%% (rispetto a bbox)\n', r.mean_displacement);
    fprintf('  Riduzione rumore:    %.1f%% (energia Laplaciana)\n', r.noise_reduction);
    fprintf('  Tempo elaborazione:  %.2f s\n', r.time_smoothing);
    fprintf('\n');
    fprintf('Output generati:\n');
    fprintf('  - %s_animation.gif\n', mesh_short_names{idx});
    fprintf('  - %s_comparison.png\n', mesh_short_names{idx});
    fprintf('========================================\n\n');
end

fprintf('Tutti i file salvati in: %s/\n', results_folder);
fprintf('Elaborazione completata!\n\n');

%% FUNCTION DI SUPPORTO

function update_animation(V_current, t_current, step)
    % Aggiorna visualizzazione mesh e salva frame GIF
    global ANIM_DATA;
    
    % Controlla se salvare questo frame
    if mod(step, ANIM_DATA.skip) ~= 0 && step ~= 1
        return;
    end
    
    % Aggiorna mesh
    set(ANIM_DATA.h_patch, 'Vertices', V_current);
    
    % Transizione colore rosso -> blu
    progress = min(max(t_current / ANIM_DATA.T_final, 0), 1);
    color_r = 1.0 - 0.4 * progress;
    color_g = 0.7 + 0.15 * progress;
    color_b = 0.7 + 0.3 * progress;
    set(ANIM_DATA.h_patch, 'FaceColor', [color_r, color_g, color_b]);
    
    % Aggiorna titolo
    set(ANIM_DATA.h_title, 'String', ...
        sprintf('%s - t = %.4f / %.4f (step %d)', ...
                ANIM_DATA.mesh_name, t_current, ANIM_DATA.T_final, step));
    drawnow;
    
    % Cattura e salva frame GIF
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    ANIM_DATA.gif_frame_count = ANIM_DATA.gif_frame_count + 1;
    
    if ANIM_DATA.gif_frame_count == 1
        % Primo frame
        imwrite(imind, cm, ANIM_DATA.gif_file, 'gif', ...
                'Loopcount', inf, 'DelayTime', 1/ANIM_DATA.framerate);
    else
        % Frame successivi
        imwrite(imind, cm, ANIM_DATA.gif_file, 'gif', ...
                'WriteMode', 'append', 'DelayTime', 1/ANIM_DATA.framerate);
    end
end