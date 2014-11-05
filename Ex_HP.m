%{
Example script for HP.m.  

For questions or comments, please contact Aubrey Gress at
agress@ucdavis.edu.
%}

function [] = Ex_HP()
    [X,W] = CreateRandomDataSet();
    options = struct();
    options.reg = .1;
    options.numVecs = 10;
        
    [Xproj,projections] = HP(X,W,options);
    D = DistanceMatrix(Xproj{1},Xproj{2});
    
    instanceIDs = Helpers.getDataSetIDs(Xproj,1);
    W12 = Helpers.getSubW(W,instanceIDs,1,2);
    meanD = mean(D(:));
    WD = W12.*D;
    meanWD = sum(WD(:))/nnz(WD);
    display(['Mean embedded inter-dataset distance: ' num2str(meanD)]);
    display(['Mean embedded inter-dataset distance with similar-to relationship: ' num2str(meanWD)]);    
end
%{
Creates two random datasets and a matrix of relationships between them.
Each dataset is a mixture of points sampled from two standard Normal distributions
with means [2 2 ... 2] and [-2 -2 ... -2]
Each normal distribution is associated with a label from which the inter-dataset
relationships are generated.  

RETURN VALUES:
X: Cell array of data sets.  Each data set is a matrix of size [n_i x p_i]
where n_i is the number of instances and p_i is the number of features

W: Inter-dataset and Intra-dataset relationship matrix.  Is size [n x n]
where n = n_1 + n_2 + ... + n_k where k is the number of data sets in X
%}
function [X,W] = CreateRandomDataSet()
    numInstances = [100 50];
    dim = [15 12];
    X = {};
    clusterIDs = {};
    numClusters = 2;
    for i=1:length(numInstances)        
        clusterCentroids = ones(numClusters,dim(i));
        clusterCentroids(1,:) = 2*clusterCentroids(1,:);
        clusterCentroids(2,:) = -2*clusterCentroids(2,:);
        IDs = randi(numClusters,1,numInstances(i));
        clusterIDs{i} = IDs;
        Xi = zeros(numInstances(i),dim(i));
        for j=1:size(Xi,1)
            Xi(j,:) = randn(1,size(Xi,2)) + clusterCentroids(IDs(j),:);
        end
        X{i} = Xi;
    end    
    W = zeros(sum(numInstances));
    clusterIDs1 = clusterIDs{1};
    clusterIDs2 = clusterIDs{2};
    for i=1:numInstances(1)
        Wi = clusterIDs1(i) == clusterIDs2;
        W(i,numInstances(1)+1:end) = Wi;
    end
    W = (W + W');
end

function [D] = DistanceMatrix(X1, X2)
    D = pdist2(X1,X2);
end