function [sim,distance] = my_siminet(matrice1, matrice2)

% Inputs:  matrice1: Adjacency matrix of the first graph
%          matrice2: Adjacency matrix of the second graph
% Outputs:  - sim: similarity score between 0 et 1  
%           - distance: distance between graph1 and graph2

    matrice1t=triu(matrice1);
    matrice2t=triu(matrice2);
    Etotal=or(matrice1t,matrice2t);
    DE=sum(sum(abs(matrice1t-matrice2t)));
    DE=DE/sum(sum(Etotal));
    distance=DE;
    sim=1/(1+distance);
    
end

