function [U_Vector, F_Vector, InternalForces] = SpringSolver(k_vector, n_nodes, f_nodes, n_elements, ext_Force, I_matrix, ext_Node)
%Creation for each element stiffness matrix
for count1 = 1:n_elements
    K_matrix(:,:,count1) = SpringElementStiffness(k_vector(count1));   
end
K = zeros(n_nodes);
%Assemble
for count2 = 1:n_elements
    indx1=1;
    K = SpringAssemble(K,K_matrix(:,:,count2),I_matrix(count2,indx1),I_matrix(count2,indx1+1));
end
%Boundary and load conditions
Sf_nodes = sort(f_nodes, 'descend'); %Order the fixed nodes 
[p1, reduce] = size(f_nodes); %Get the number of nodes to eliminate
k = K;
%Reducing the matrix
for count3 = 1:reduce
    k = MatrixErase(k,Sf_nodes(count3),Sf_nodes(count3));
end
%Creating part of global displacement vector and auxiliary force vector
vc = 1;
for count4 = 1:n_nodes
    verify = false;
    for count5 = 1:reduce
        if count4 == f_nodes(count5) && verify == false
            verify=true;
            U(count4) = 0;
        end
    end
    if verify ~= true
        if count4 == ext_Node
                f(vc)= ext_Force;
                vc = vc+1;
            else
                f(vc) = 0;
                vc = vc+1;
        end
    end         
end
%Solving the auxiliary displacement vector 
f_t = f';
u = k\f_t;
vcc=1;
%Finishing global displacement vector
for count6 = 1:n_nodes
   verify = false;
   for count7 = 1:reduce
        if count6 == f_nodes(count7) && verify == false
            verify=true;
            U(count6) = 0;
        end
   end
   if verify ~= true
       U(count6) = u(vcc);
       vcc = vcc+1;
   end
end
%Solving the auxiliary displacement vector
U_Vector = U';
F_Vector = K*U_Vector;
%Finding the forces in each spring
for count8 = 1:n_elements
    gg=1;
    InternalForces(:,:,count8) = K_matrix(:,:,count8)*[U_Vector(I_matrix(count8,gg));U_Vector(I_matrix(count8,gg+1))];
end