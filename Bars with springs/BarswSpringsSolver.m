function [U_Vector, F_Vector, S_InternalForces, B_InternalForces, B_InternalStresses] = BarswSpringsSolver(k_vector, A_vector, E_vector, T_vector, bar_coords, n_nodes, f_nodes, n_springs, n_bars, ext_Force, I_matrix, ext_Node)
%Finding the length of each bar 
for count1 = 1:n_bars
    L_vector(count1)= PlaneTrussElementLength(bar_coords(count1, 1),bar_coords(count1, 2),bar_coords(count1, 3),bar_coords(count1, 4));
end
%Creation for each element stiffness matrix
t_count = 0;
for count2 = 1:n_bars
    BK_matrix(:,:,count2) = PlaneTrussElementStiffness(E_vector(count2), A_vector(count2), L_vector(count2), T_vector(count2));   
    t_count = t_count+1;
end
%Creation for each element stiffness matrix
for count3 = 1:n_springs
    SK_matrix(:,:,count3) = PlaneTrussElementStiffness(k_vector(count3),1,1,T_vector(count3+t_count));   
end
K = zeros(n_nodes*2);
%Assemble
indx1=0;
for count4 = 1:n_bars
    K = PlaneTrussAssemble(K,BK_matrix(:,:,count4),I_matrix(count4,1),I_matrix(count4,2));
    indx1=indx1+1;
end
for count5 = 1:n_springs
    K = PlaneTrussAssemble(K,SK_matrix(:,:,count5),I_matrix((indx1+count5),1),I_matrix((indx1+count5),2));
end
%Boundary and load conditions
Sf_nodes = sort(f_nodes, 'descend'); %Order the fixed nodes 
[p1, reduce] = size(f_nodes); %Get the number of nodes to eliminate
k = K;
%Reducing the matrix
for count6 = 1:reduce
    k = MatrixErase(k,Sf_nodes(count6),Sf_nodes(count6));
end
%Creating part of global displacement vector and auxiliary force vector
vc = 1;
for count7 = 1:(n_nodes*2)
    verify = false;
    for count8 = 1:reduce
        if count7 == f_nodes(count8) && verify == false
            verify=true;
            U(count7) = 0;
        end
    end
    if verify ~= true
        if count7 == ext_Node*2
            f(vc)= ext_Force(2);
            vc = vc+1;
        elseif count7 == (ext_Node*2-1)
            f(vc)= ext_Force(1);
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
for count9 = 1:(n_nodes*2)
   verify = false;
   for count10 = 1:reduce
        if count9 == f_nodes(count10) && verify == false
            verify=true;
            U(count9) = 0;
        end
   end
   if verify ~= true
       U(count9) = u(vcc);
       vcc = vcc+1;
   end
end
%Solving the auxiliary displacement vector
U_Vector = U';
F_Vector = K*U_Vector;
%Finding the forces in each bar
for count11 = 1:n_bars
    gg=1;
    B_InternalForces(count11) = PlaneTrussElementForce(E_vector(count11), A_vector(count11), L_vector(count11), T_vector(count11),[U_Vector((I_matrix(count11,gg)*2-1));U_Vector((I_matrix(count11,gg)*2));U_Vector((I_matrix(count11,gg+1)*2-1));U_Vector((I_matrix(count11,gg+1)*2))]);
end
%Finding the stresses in each bar
for count12 = 1:n_bars
    gg=1;
    B_InternalStresses(count12) = PlaneTrussElementStress(E_vector(count12), L_vector(count12), T_vector(count12),[U_Vector((I_matrix(count12,gg)*2-1));U_Vector((I_matrix(count12,gg)*2));U_Vector((I_matrix(count12,gg+1)*2-1));U_Vector((I_matrix(count12,gg+1)*2))]);
end
for count13 = 1:n_springs
    gg=1;
    S_InternalForces(:,:,count13) = SpringElementForces(k_vector(count13), [U_Vector(I_matrix(indx1+count13,gg));U_Vector(I_matrix(indx1+count13,gg+1))]);
end