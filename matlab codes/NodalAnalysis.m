% Parse and generate MNA matrices from the netlist with the given filename.
function [M, C, b] = NodalAnalysis(filename)
%function [Mupper, Mlower] = NodalAnalysis(filename)

G = zeros; C = zeros; L = zeros; Av = zeros; AL = zeros; Is = zeros; Vs = zeros; Avcvs = zeros; Bvcvs = zeros;
Nv = 0; Nvcvs = 0; NL = 0;

f = fopen(filename, 'r');

while 1
    % First char determines type of element:
    type = fscanf(f, '%c', 1);
    switch type
        case 'R'
            x = fscanf(f, '%*s %d %d %f ', 3);
            %disp('R'); disp(x'); % x(1) = node1, x(2) = node2, x(3) = value
            % Do something with R...
            n1 = x(1); n2 = x(2); R = x(3);
            
            m = max(n1, n2);
            if m > length(G), G(m, m) = 0; end % Expand if needed
            
            if (n1 > 0)&&(n2 > 0)
                G(n1, n1) = G(n1, n1) + 1/R; G(n2, n1) = G(n2, n1) - 1/R; 
                G(n1, n2) = G(n1, n2) - 1/R; G(n2, n2) = G(n2, n2) + 1/R; 
            elseif n1 > 0
                G(n1, n1) = G(n1, n1) + 1/R;
            else
                G(n2, n2) = G(n2, n2) + 1/R;
            end
        
        case 'C'
            x = fscanf(f, '%*s %d %d %f ', 3);
            n1 = x(1); n2 = x(2); Cap = x(3);                        
            m = max(n1, n2);
            if m > length(C), C(m, m) = 0; end % Expand if needed
            
            if (n1 > 0)&&(n2 > 0)
                C(n1, n1) = C(n1, n1) + Cap; C(n2, n1) = C(n2, n1) - Cap; 
                C(n1, n2) = C(n1, n2) - Cap; C(n2, n2) = C(n2, n2) + Cap; 
            elseif n1 > 0
                C(n1, n1) = C(n1, n1) + Cap;
            else
                C(n2, n2) = C(n2, n2) + Cap;
            end
            
        case 'I'
            x = fscanf(f, '%*s %d %d DC %f ', 3);
            %disp('I'); disp(x'); % x(1) = node1, x(2) = node2, x(3) = value
            % Do something with I...
            n1 = x(1); n2 = x(2); isb = x(3);
            
            m = max(n1, n2);
            if m > length(Is), Is(m) = 0; end % Expand if needed
            
            if n1 > 0, Is(n1) = Is(n1) - isb; end
            if n2 > 0, Is(n2) = Is(n2) + isb; end
            
        case 'L'
            x = fscanf(f, '%*s %d %d %f ', 3);
            n1 = x(1); n2 = x(2); ind = x(3);
            
            NL = NL + 1;
            if n1 > 0, AL(n1, NL) = 1; end
            if n2 > 0, AL(n2, NL) = -1; end
            L(NL) = ind;
            
        case 'V'
            x = fscanf(f, '%*s %d %d DC %f ', 3);
            %disp('V'); disp(x'); % x(1) = node+, x(2) = node-, x(3) = value
            % Do something with V...
            Nv = Nv + 1;
            n1 = x(1); n2 = x(2); vsb = x(3);
            if n1 > 0, Av(n1, Nv) = -1; end
            if n2 > 0, Av(n2, Nv) = 1; end
            Vs(Nv) = vsb;
            
        case 'E'
            x = fscanf(f, '%*s %d %d %d %d %f ', 5);
            %disp('E'); disp(x'); % x(1) = node+, x(2) = node-, x(3) = nodectrl+, x(4) = nodectrl-, x(5) = gain
            % Do something with E...
            Nvcvs = Nvcvs + 1;
            n1 = x(1); n2 = x(2); nc1 = x(3); nc2 = x(4); gain = x(5);

            m = max(nc1, nc2);
            [hBvcvs, wBvcvs] = size(Bvcvs); %#ok<ASGLU>
            if m > wBvcvs, Bvcvs(Nvcvs, m) = 0; end % Expand if needed
            
            if n1 > 0, Avcvs(n1, Nvcvs) = -1; Bvcvs(Nvcvs, n1) = 1; end
            if n2 > 0, Avcvs(n2, Nvcvs) = 1; Bvcvs(Nvcvs, n2) = -1; end
            if nc1 > 0, Bvcvs(Nvcvs, nc1) = Bvcvs(Nvcvs, nc1) - gain; end
            if nc2 > 0, Bvcvs(Nvcvs, nc2) = Bvcvs(Nvcvs, nc2) + gain; end
            
        otherwise
            break % Assume end of file
    end
end
fclose(f);

N = length(G);

b = Is';
if length(b) < N, b = padarray(b, N-length(b), 0, 'pos'); end

if Nv > 0
    b = [b; Vs'];
    [hAv, wAv] = size(Av);  %#ok < NASGU > 
    if hAv < N, Av = padarray(Av, N-hAv, 0, 'pos'); end
end

if Nvcvs > 0
    b = padarray(b, Nvcvs, 0, 'pos');
    
    [hAvcvs, wAvcvs] = size(Avcvs);   %#ok < NASGU > 
    if hAvcvs < N
        Avcvs = padarray(Avcvs, N-hAvcvs, 0, 'pos'); 
    end
    
    [hBvcvs, wBvcvs] = size(Bvcvs);   %#ok < ASGLU > 
    if wBvcvs < N
        Bvcvs = padarray(Bvcvs', N-wBvcvs, 0, 'pos'); 
        Bvcvs = Bvcvs'; 
    end

    %%%%%
    if Nv > 0
        Mupper = [G Av Avcvs];
        Mlower = [[-Av'; Bvcvs] zeros(Nv+Nvcvs)];
    else
        Mupper = [G Avcvs];
        Mlower = [Bvcvs zeros(Nvcvs)];
    end
    
    M = [Mupper; Mlower];
    %%%%%
else
    if Nv > 0
        M = [[G Av]; [-Av' zeros(Nv)]];
    else
        M = G;
    end
end

if NL > 0
    b = padarray(b, NL, 0, 'pos');

    [hAL, wAL] = size(AL);  %#ok < NASGU > 

    if hAL < length(M), AL = padarray(AL, length(M)-hAL, 0, 'pos'); end
    M = [[M AL]; [AL' zeros(NL)]];

    m = length(M);
    C(m, m) = 0;
    C = C + diag( [zeros(1, m-NL) -L] );
end
