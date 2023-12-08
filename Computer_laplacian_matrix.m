function L= Computer_laplacian_matrix(W)

if W == W'    
    fprintf('Run the first case undirected graph')
    % out-degree matrix
    D=diag(sum(W,2));                
    Wsymm=D^(-1/2)*W*D^(-1/2);       
    Wnew=Wsymm;
    n_ROI = size(W,1);
    % compute normalized undirected graph Laplacian: L^=D^(-1/2)*(D-W)*D^(-1/2)
    L=eye(n_ROI)-Wnew;      
       
else
    fprintf('Run the second case directed graph')
    % verify strong-connected graph
    P = verify_stron_con(W);
    
    % % =========================================================================
    % %  (1) Method 1: https://epfl-lts2.github.io/gspbox-html/doc/utils/gsp_create_laplacian.html  
    % % =========================================================================
%     G = gsp_graph(W);
%     G = gsp_create_laplacian(G, 'chung');
%     L_dir = G.L;
%     L = L_dir;
%     phi = G.phi;
    
    % % =========================================================================
    % %  (2) Method 2
    % % =========================================================================
    % out-degree matrix
    n_ROI = size(W,1);
    D=diag(sum(W,2)); 
    % transition matrix of random walk
    P = inv(D)* W;
    % Create Markov chain
    mc = dtmc(P);
    numstates = mc.NumStates;
    % Check Markov chain for ergodicity
    markov_ergodic = isergodic(mc);
    
    % Visually confirm that the Markov chain is not ergodic by plotting its eigenvalues on the complex plane.
    % h = figure;
    % eigplot(mc);
     
    % Determine Markov chain asymptotics
    xFix = asymptotics(mc);
    XFix = diag(xFix);
    Wnew = (XFix^(1/2)*P*XFix^(-1/2) + XFix^(-1/2)*P'*XFix^(1/2))./2;
    % Chung normalized directed graph Laplacian
    L = eye(n_ROI)-Wnew;
    
end
end