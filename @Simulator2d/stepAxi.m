function s = stepAxi(s,dt,comp,steptype,Re)

%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%SIMULATOR simulator class constructor.
%   s = Simulator(m) creates a simulator object from the mesh object

%Name: stepAxi
%Location: <path>/@Simulator
%Purpose: this is the main program,

% modificado em 13/03/2006
% revisado   em 09/04/2007

IEN = getIEN(s.m);
X= getX(s.m);
Y=getY(s.m);
Z=getZ(s.m);

nelem=size(IEN,1);
nvert=size(X,1)-nelem;
nnodes=size(X,1);

velu=s.us;
velv=s.vs;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%  alpha = 0   -> explicito                                     %
%  alpha = 0.5 -> crank-nicholson                               %
%  alpha = 1   -> implicito                                     %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

alpha=1;

%% coeficiente do laplaciano - viscosidade e concentracao
k=1/Re;

%%% calculo do termo convectivo para o met semi-lagrangeano
%[up,vp] = convectLin(s,dt);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% semi-lagrangiano  --  montagem dos vetores e matriz           %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

%va=((1/dt)*s.M-(1-alpha)*k*s.K)*[up;vp]; 
va=((1/dt)*s.M-(1-alpha)*k*s.K)*[up;vp]-s.G*s.ps; % com correcao

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Lagrangiano  --  montagem dos vetores e matriz           %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%va=((1/dt)*s.M-(1-alpha)*k*s.K)*[velu;velv]-s.G*s.ps;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% metodo acoplado                                               %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if(strcmp(steptype,'coupled'))

    mat=1/dt*s.M+alpha*k*s.K;
    % A=[K G]=[Ku  0 Gx][ ]=[us]  [us]=[ ]
    %   [D 0] [ 0 Kv Gy][u]=[vs]  [vs]=[b]
    %         [Dx Dy  0][ ]=[ps]  [ps]=[ ]

    A=sparse([mat s.G; s.D 0*(s.D*s.G)]);
    b=[va;zeros(nvert,1)];
    [A b] = setCoupledBC(s,A,b);
    u=A\b;
    us=u(1:nnodes,1);
    vs=u(1+nnodes:nnodes*2,1);
    ps=u(1+2*nnodes:2*nnodes+nvert,1);
    s.uant=u;

end;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Metodo da Projecao discreto baseado em decomposicao LU        %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if(strcmp(steptype,'uncoupled'))

    if(comp)
        mat=1/dt*s.M+alpha*k*s.K;
        [At G D E b1 b2 ip] = setUncoupledBC(s,mat,va*0);

        inva=diag(sparse(1./sum(mat,2)));
        % uzawa
        %inva=inv(At);

        E=E-D*inva*G;

        s.At=At;
        s.Gt=G;
        s.Dt=D;
        s.Et=E;
        s.b1=b1;
        s.b2=b2;
        s.ip=ip;
        s.inva=inva;
    else
        inva=s.inva;
        At=s.At;
        G=s.Gt;
        D=s.Dt;
        E=s.Et;
        b1=s.b1;
        b2=s.b2;
        ip=s.ip;

    end;

    b1=b1+va.*ip;

    ut=At\b1;

    % [A A G][]=[b1]
    % [A A G][]=[b1]
    % [D D E][]=[b2]

    %E=E-D*inva*G;
    b2=b2-D*ut;

    pt=E\b2;
    
    ua=ut-inva*G*pt;

    u=[ua;pt];
    us=u(1:nnodes,1);
    vs=u(1+nnodes:nnodes*2,1);
    ps=u(1+2*nnodes:2*nnodes+nvert,1);
    s.uant=u;
end;

s.us=us;
s.vs=vs;

%s.ps=ps;

% correcao na pressao
% SETUNCOUPLEDBC deve ser ajustado para b2 = 0
s.ps=s.ps+ps; % correcao na pressao

s.time=s.time+dt;

