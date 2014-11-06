%{
Solves the optimization problem proposed in [Gress and Davidson 2014].

For questions or comments, please contact Aubrey Gress at
agress@ucdavis.edu.

The optimization problem is:

min         p'ZLZ'p
such that   a_1' X_1 D_11 X_1' a_1 = 1

This is solved by first solvering for a_1:

a_1 + lambda V_11 a_1 = 0

where:

    lambda: smallest eigenvalue of V
    V_11: Is the square upper left block of V
    V = (ZLZ')^-1 Z B Z'
    B = [D_11 0 ... 0]
        [0    0 ... 0]
        [............]
        [0 ........ 0]

(i.e. B is D_11 padded with zeros so it's the proper size)

After solving for a_1, the remaining projection vectors are solved for by:

[a_2]   
[...] = V_21 a_1
[a_k]   

where:
    v_21: Is the remaining lower left block of V

PARAMETERS:
X: Cell array of data sets.  The i'th data set is a matrix of size [n_i x p_i]
where n_i is the number of instances and p_i is the number of features

W: Inter-dataset and intra-dataset relationship matrix.  Is size [n x n]
where n = n_1 + n_2 + ... + n_k where k is the number of data sets in X.
Solves assumes it is symmetric.

options: Struct of options.  Can have the following fields:
    reg: Regularization parameters.  Defaults to 0. In [Gress and Davidson 2014] this was
    selected using cross validation.

    numVecs: number of projection vectors to use.  Defaults to using all 
    projection vectorsIn [Gress and Davidson 2014] this was selected using cross validation.

RETURN VALUES:
Xproj: Cell array of projected datasets.

projections: Cell array of projections vectors
%}


function [Xproj,projections] = HP(X, W, options) 
    if ~exist('options','var')
        options = struct();
    end
    if ~isfield(options,'reg')
        options.reg = 0;
    end
    if ~isfield(options,'numVecs')
        options.numVecs = inf;
    end
    assert(isequal(W,W'));
    assert(size(W,1) == size(W,2));
    
    XBlock = Helpers.makeBlockMatrix(X);
    assert(size(XBlock,1) == size(W,1));
    L = makeLaplacian(W);    
    instanceIDs = Helpers.getDataSetIDs(X,1);
    D11 = diag(sum(Helpers.getSubW(W,instanceIDs,1,1),2));
    
    X1 = X{1};      
    Q = XBlock'*L*XBlock;    
    B = X1'*D11*X1;
    B = B + options.reg*eye(size(B));
    [numRows,numCols] = size(XBlock);
    Bpadded = Helpers.padMatrix(B,numCols,numCols);
    Qinv = inv(Q);
    V = Qinv*Bpadded;                  
    
    X1dim = size(X1,2);
    V11 = V(1:X1dim,1:X1dim);    
    V21 = V(X1dim+1:end,1:X1dim);
    [A1,eigVals] = eig(V11);
    [sortedEigVals,I] = sort(diag(eigVals),'descend');
    A1 = A1(:,I);
        
    for i=1:size(A1,2)
        ai = A1(:,i);
        v = ai'*B*ai;
        A1(:,i) = ai/sqrt(v);
    end
    
    A2_m = V21*A1;
    Aall = [A1 ; A2_m];
    Xproj = cell(size(X));
    projections = cell(size(X));
    featureIDs = Helpers.getDataSetIDs(X,2);  
    
    for i=1:length(X)
        Ai = Aall(featureIDs==i,:);
        if ~isinf(options.numVecs);
            Ai = Ai(:,1:options.numVecs);
        end
        if i > 1
            numVecs = size(Ai,2);
            lambdas = inv(diag(sortedEigVals(1:numVecs)));
            for j=1:numVecs
                ai_j = Ai(:,j);
                Ai(:,j) = ai_j*lambdas(j,j);
            end
        end
        Xproj{i} = X{i}*Ai;
        projections{i} = Ai;
    end
    
end

function [L] = makeLaplacian(W)
    L = diag(sum(W,2)) - W;
end
