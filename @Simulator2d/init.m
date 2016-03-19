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
K = sparse(nnodes,nnodes);
M = sparse(nnodes,nnodes);
Mz = sparse(nnodes,nnodes);
G1 = sparse(nnodes,nvert);
G2 = sparse(nnodes,nvert);

s.uant=sparse(nnodes*2+nvert,1);

%%% determinacao do tipo de elemento na classe TElement
element=FEMMiniElement2d();

for mele = 1:nele

    v1=IEN(mele,1);
    v2=IEN(mele,2);
    v3=IEN(mele,3);
    mele/nele

    [massele,kxx,kyy,kxy,kyx,gxele,gyele,ngleu,nglep,v]=getmgq(element,mele,IEN,X,Y,Z);

    vp = v(1:3);
    %ngleu: numero de graus de liberdade do elemento associados a velocidade
    %nglep: idem associados a p

    % K= [11 12]
    % M= [21 22]

    K(v,v) = K(v,v)+(kxx+kyy);
    M(v,v) = M(v,v)+massele;

    G1(v,vp)=G1(v,vp)+gxele;
    G2(v,vp)=G2(v,vp)+gyele;

end;

s.nnodes=nnodes;
s.nvert=nvert;
s.M=[ M Mz;Mz M ];
s.K=[ K Mz; Mz K ];
s.G=[ G1;G2 ];
s.D=s.G';

s.us=uc;
s.vs=vc;
s.ps=pc;

s.convlin=sparse(nvert,nvert);

