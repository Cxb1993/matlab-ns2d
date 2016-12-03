function s = init(s)
%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%SIMULATOR simulator class constructor.
%   s = Simulator(m) creates a simulator object from the mesh object

%Name: init
%Location: <path>/@Simulator
%Purpose: initialize the global matrix from elementary matrix

% modificado em 13/01/2006
% revisado   em 18/04/2006

[uc vc pc]=getBC(s.m);

IEN = getIEN(s.m);
X= getX(s.m);
Y=getY(s.m);
Z=getZ(s.m);

nele=size(IEN,1);
nnodes=size(X,1);
nvert=nnodes-nele;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% alocacao de memoria para matrizes e vetores                   %
% K, M, G e D montados para u e v                               %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
K11 = sparse(nnodes,nnodes);
K22 = sparse(nnodes,nnodes);
M = sparse(nnodes,nnodes);
Mz = sparse(nnodes,nnodes);
G1 = sparse(nnodes,nvert);
G2 = sparse(nnodes,nvert);
D1 = sparse(nvert,nnodes);
D2 = sparse(nvert,nnodes);

s.uant=sparse(nnodes*2+nvert,1);

%%% determinacao do tipo de elemento na classe TElement
element=FEMMiniElement2d();

for mele = 1:nele

    v1=IEN(mele,1);
    v2=IEN(mele,2);
    v3=IEN(mele,3);
    mele/nele

	radius = Y[v1]+Y[v2]+Y[v3]/3.0;

    [massele,kxx,kyy,kxy,kyx,kx,ky,gxele,gyele,dxele,dyele,dmass,ngleu,nglep,v]=getmgqAxi(element,mele,IEN,X,Y,Z);

    vp = v(1:3);
    %ngleu: numero de graus de liberdade do elemento associados a velocidade
    %nglep: idem associados a p

    % K= [11 0]
    % M= [0 22]

    K11(v,v) = K11(v,v)+(kxx+kyy-(1.0/radius)*ky);
    K22(v,v) = K22(v,v)+(kxx+kyy-(1.0/radius)*ky+(1.0/(radius*radius))*massele);
    M(v,v) = M(v,v)+massele;

    G1(v,vp)=G1(v,vp)+gxele;
    G2(v,vp)=G2(v,vp)+gyele;

    D1(vp,v)=G1(vp,v)+dxele;
    D2(vp,v)=G2(vp,v)+dyele+(1.0/radius)*dmass;

end;

s.nnodes=nnodes;
s.nvert=nvert;
s.M=[ M Mz;Mz M ];
s.K=[ K11 Mz; Mz K22 ];
s.G=[ G1;G2 ];
s.D=[ D1 D2 ];

s.us=uc;
s.vs=vc;
s.ps=pc;

s.convlin=sparse(nvert,nvert);

