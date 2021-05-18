function [U_Vector, F_Vector, InternalForces] = SpringSolver(k_vector, n_nodes, f_nodes, n_elements, ext_Force, I_matrix, ext_Node)

for count1 = 1:n_elements
    K_matrix(count1) = SpringElementStiffness(k_vector(count1));   
end

K = zeros(n_nodes);

for count2 = 1:n_elements
    g=1;
    K = SpringAssemble(K,K_matrix(count2),I_matrix(count2,g),I_matrix(count2,g+1));
end
Sf_nodes = sort(f_nodes, descend);
aux = size(f_nodes);
k = K;
for count3 = 1:aux
    k = MatrixErase(k,Sf_nodes(count3),Sf_nodes(count3));
end
vc = 1;
for count4 = 1:n_nodes
    for count5 = 1:aux
        if count4 ~= f_nodes(count5)
            if count4 == ext_Node
                f(vc)= ext_Force;
                U(count4) = ext_Force;
                vc = vc+1;
            else
                f(vc) = 0;
                U(count4) = 0;
                vc = vc+1;
            end
        end
    end
end
f_t = f';
u = k\f_t;

for count6 = l:n_nodes
    for count7 = 1:aux
        if count6 == f_nodes(count7)
            U(count6) = u(count7);
        end
    end
end
U_Vector = U';
F_Vector = K*U_Vector;

for count8 = 1:n_elements
    gg=1
    InternalForces(count8) = K_matrix(count8)*[U_Vector(I_matrix(count8,gg)),U_Vector(I_matrix(count8,gg+1))];
end

