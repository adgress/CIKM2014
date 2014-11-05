classdef Helpers
    %HELPERS Helper functions for HP.m and Ex_HP.m [Gress and Davidson 2014].
    
    properties
    end
    
    methods(Static)
        function [Xpadded] = padMatrix(X,rows,cols)
            Xpadded = zeros(rows,cols);
            Xpadded(1:size(X,1),1:size(X,2)) = X;
        end
        
        function [Wij] = getSubW(W,dataSetIDs,i,j)
            Wij = W(dataSetIDs==i,dataSetIDs==j);
        end
        
        function [Xblock] = makeBlockMatrix(X)
            [numRows,numCols] = Helpers.getSizeOfBlockDiagonalMatrix(X);
            Xblock = zeros(numRows,numCols);
            rowIdx=0;
            colIdx=0;
            for i=1:length(X)
                Xi = X{i};
                XiRows = size(Xi,1);
                XiCols = size(Xi,2);
                Xblock(rowIdx+1:rowIdx+XiRows,colIdx+1:colIdx+XiCols) = Xi;
                rowIdx = rowIdx + XiRows;
                colIdx = colIdx + XiCols;
            end
        end
        
        function [dataSetIDs] = getDataSetIDs(X,dim)
            sizes = [0 0];
            [sizes(1),sizes(2)] = Helpers.getSizeOfBlockDiagonalMatrix(X);
            dataSetIDs = zeros(sizes(dim),1);
            numInstances = [0 cumsum(cellfun(@(x) size(x,dim),X))];
            for i=1:length(numInstances)-1
                dataSetIDs(numInstances(i)+1:numInstances(i+1)) = i;
            end
        end
        
        function [numRows,numCols] = getSizeOfBlockDiagonalMatrix(X)
            numRows = sum(cellfun(@(x) size(x,1),X));
            numCols = sum(cellfun(@(x) size(x,2),X));
        end
    end
    
end

