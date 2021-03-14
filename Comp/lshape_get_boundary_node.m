%Function which tells us when a point is on the boundary of the domain.
%If not on boundary it returns 0. If the boundary is Dirichlet
%returns 1, if the boundary is Neumann returns 2.

function boundary_node = lshape_get_boundary_node(coords,limits)

%Tolerance
eps = 1.0E-7;

if (abs(coords(1)-limits(1,1))<eps || abs(coords(2)-limits(2,1)) <eps || abs(coords(2)-limits(2,2)) <eps ...
         || (abs(coords(1)-(limits(1,1)+.5*(limits(1,2)-limits(1,1))))< eps && coords(2)>=limits(2,1)+.5*(limits(2,2)-limits(2,1)) ||...
        abs(coords(2)-(limits(2,1)+.5*(limits(2,2)-limits(2,1))))< eps && coords(1)<=limits(1,1)+.5*(limits(1,2)-limits(1,1))) || abs(coords(1)-1.0) < eps)
    %Dirichlet Boundary
    boundary_node = 1;
%elseif (abs(coords(1)-1.0) < eps)
%    %Neumann Boundary
%    boundary_node = 2;
else
    %Not on boundary
    boundary_node = 0;
end

end