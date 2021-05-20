function [U_Vector, F_Vector, InternalForces, InternalStresses] = Bar3DSolver(A_vector, E_vector, bar_coords, n_nodes, f_nodes, n_elements, ext_Force, I_matrix, ext_Node)
format shortG
%Finding the length and angles in x,y,z of each bar 
for count1 = 1:n_elements
    L_vector(count1)= SpaceTrussElementLength(bar_coords(count1, 1),bar_coords(count1, 2),bar_coords(count1, 3),bar_coords(count1, 4), bar_coords(count1, 5), bar_coords(count1, 6));
    Thetax(count1) = acosd( (bar_coords(count1, 1) - bar_coords(count1, 4))/ L_vector(count1) ); 
    Thetay(count1) = acosd( (bar_coords(count1, 2) - bar_coords(count1, 5))/ L_vector(count1) ); 
    Thetaz(count1) = acosd( (bar_coords(count1, 3) - bar_coords(count1, 6))/ L_vector(count1) ); 
end
%Creation for each element stiffness matrix
for count2 = 1:n_elements
    K_matrix(:,:,count2) = SpaceTrussElementStiffness(E_vector(count2), A_vector(count2), L_vector(count2), Thetax(count2), Thetay(count2), Thetaz(count2));   
end
K = zeros(n_nodes*3);
%Assemble
for count3 = 1:n_elements
    indx1=1;
    K = SpaceTrussAssemble(K,K_matrix(:,:,count3),I_matrix(count3,indx1),I_matrix(count3,indx1+1));
end
%Boundary and load conditions
Sf_nodes = sort(f_nodes, 'descend'); %Order the fixed nodes 
[p1, reduce] = size(f_nodes); %Get the number of nodes to eliminate
k = K;
%Reducing the matrix
for count4 = 1:reduce
    k = MatrixErase(k,Sf_nodes(count4),Sf_nodes(count4));
end
%Creating part of global displacement vector and auxiliary force vector
vc = 1;
for count5 = 1:(n_nodes*3)
    verify = false;
    for count6 = 1:reduce
        if count5 == f_nodes(count6) && verify == false
            verify=true;
            U(count5) = 0;
        end
    end
    if verify ~= true
        if count5 == ext_Node*3
            f(vc)= ext_Force(3);
            vc = vc+1;
        elseif count5 == (ext_Node*3-1)
            f(vc)= ext_Force(2);
            vc = vc+1;
        elseif count5 == (ext_Node*3-2)
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
for count7 = 1:(n_nodes*3)
   verify = false;
   for count8 = 1:reduce
        if count7 == f_nodes(count8) && verify == false
            verify=true;
            U(count7) = 0;
        end
   end
   if verify ~= true
       U(count7) = u(vcc);
       vcc = vcc+1;
   end
end
%Solving the auxiliary displacement vector
U_Vector = U';
F_Vector = K*U_Vector;
%Finding the forces in each bar
for count9 = 1:n_elements
    gg=1;
    InternalForces(count9) = SpaceTrussElementForce(E_vector(count9), A_vector(count9), L_vector(count9), Thetax(count9), Thetay(count9), Thetaz(count9), [U_Vector((I_matrix(count9,gg)*3-2));U_Vector((I_matrix(count9,gg)*3-1));U_Vector((I_matrix(count9,gg)*3));U_Vector((I_matrix(count9,gg+1)*3-2));U_Vector((I_matrix(count9,gg+1)*3-1));U_Vector((I_matrix(count9,gg+1)*3))]);
end
%Finding the stresses in each bar
for count10 = 1:n_elements
    gg=1;
    InternalStresses(count10) = SpaceTrussElementStress(E_vector(count10), L_vector(count10), Thetax(count10), Thetay(count10), Thetaz(count10),[U_Vector((I_matrix(count9,gg)*3-2));U_Vector((I_matrix(count9,gg)*3-1));U_Vector((I_matrix(count9,gg)*3));U_Vector((I_matrix(count9,gg+1)*3-2));U_Vector((I_matrix(count9,gg+1)*3-1));U_Vector((I_matrix(count9,gg+1)*3))]);
end