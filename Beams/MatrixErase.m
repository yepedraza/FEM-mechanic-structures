function [MatrixErased] = MatrixErase(Matrix,row,column)

[rows,columns]= size(Matrix);
min_rows = min(rows);
max_rows = max(rows);

if max_rows > rows || min_rows < 0
    
    disp("Ingrese una matriz cuadrada");
end

MatrixErased = Matrix;
MatrixErased(row,:) = [];
MatrixErased(:,column) = [];