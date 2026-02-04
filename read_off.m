function [V, F] = read_off(filename)
    % Legge file formato .off (Object File Format)
    % Input: filename (percorso al file .off)
    % Output: V (matrice N×3 vertici), F (matrice M×3 facce)
    
    fid = fopen(filename, 'r');
    if fid == -1
        error('Impossibile aprire file: %s', filename);
    end
    
    % Verifica header OFF
    header = strtrim(fgetl(fid));
    if ~strcmp(header, 'OFF')
        fclose(fid);
        error('File non formato OFF valido (header: %s)', header);
    end
    
    % Leggi numero vertici e facce
    counts = fscanf(fid, '%d %d %d', 3);
    nV = counts(1);
    nF = counts(2);
    
    fprintf('Lettura %s: %d vertici, %d facce\n', filename, nV, nF);
    
    % Leggi vertici
    V = fscanf(fid, '%f %f %f', [3, nV])';
    
    % Leggi facce
    F = zeros(nF, 3);
    for i = 1:nF
        line = fscanf(fid, '%d', 4);
        
        if isempty(line) || length(line) < 4
            warning('Riga %d incompleta', i);
            continue;
        end
        
        n_verts = line(1);
        if n_verts ~= 3
            fclose(fid);
            error('Faccia %d ha %d vertici (solo triangoli supportati)', i, n_verts);
        end
        
        % Converti da 0-based a 1-based MATLAB
        F(i, :) = line(2:4)' + 1;
    end
    
    fclose(fid);
    
    % Verifica validità
    if max(F(:)) > nV
        error('Indici facce fuori range');
    end
    
    fprintf('Lettura completata\n');
end