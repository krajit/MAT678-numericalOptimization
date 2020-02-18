
function out=laplaceEqn(U)
nx = size(U,1);
ny = size(U,2);
%Number of steps in space(y)       
dx=1/(nx-1);                     %Width of space step(x)
dy=1/(ny-1);                     %Width of space step(y)
x=0:dx:1;                        %Range of x(0,2) and specifying the grid points
y=0:dy:1;                        %Range of y(0,2) and specifying the grid points

u=zeros(nx,ny);
e=zeros(nx,ny);
%U=zeros(nx,ny);   



%Calculating the coefficient matrix for the implicit scheme
Ex=sparse(2:nx,1:nx-1,1,nx,nx);
Cx=sparse(1:1,2:2,1,nx,nx);
Dx=sparse(nx:nx,nx-1:nx-1,1,nx,nx);
Ax=Ex+Ex'-2*speye(nx)+Cx+Dx;        
Ey=sparse(2:ny,1:ny-1,1,ny,ny);
Ay=Ey+Ey'-2*speye(ny)+Cx+Dx; 
Bx=-1*speye(nx);
%B(1,nx+1)=2/dx^2;B((nx)*(nx),(nx)*(nx)-(nx))=2/dy^2;
A=-kron(Ay/dy^2,speye(nx))-kron(speye(ny),Ax/dx^2)-kron(Bx,speye(nx));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=@(x,y)(1-min(1,max(0,(12.*((x-0.5).*(x-0.5)+(y-0.5).*(y-0.5)))-1.0/3)) ); 
%f=@(x,y)(20.*cos(3.*pi.*x).*sin(2.*pi.*y));
S = f(x',y)+U;

%Source term
S=reshape(S,[],1);
S=A\S;
S=reshape(S,nx,ny);
u(1:nx,1:ny)=S;
%Boundary conditions

M =@(x,y)(-142.0/3+12.*((x-0.5).*(x-0.5)+(y-0.5).*(y-0.5))) ;
yd=M(x',y);

e(1,:)=-12;
e(nx,:)=-12;
e(:,1)=-12;
e(:,ny)=-12;

 out= sum(sum(0.5*(u - yd).^2*dx*dy)+sum((e.*u)*dx)+0.5*sum((U).^2)*dx*dy);
%surf(x,y,u)


end