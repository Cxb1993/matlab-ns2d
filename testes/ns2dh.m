path('/home/gustavo/sandbox/NS2d',path)
workDir='/gesar/LEN/gustavo/simulacao2dh/';
workDirVtk='/gesar/LEN/gustavo/simulacao2dh/vtk/';

Re=10000;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% utilizacao:                                               %             
% test: test (step,cavity,couette)                          %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

m1=Model2dh();
m1=test(m1,29,15,'step');

%modelOut(m1,workDirVtk,'mesh');
%show(m1);

s1=Simulator2dh(m1)
s1=init(s1);

cfl=1.00;
dt=cfl*sqrt((max(m1.Y)-min(m1.Y))*(max(m1.X)-min(m1.X))/s1.nvert)/max(s1.us);

s1=step(s1,dt,true,'uncoupled',Re);
saveDump(s1,workDir,'sim',1)

for i=1:60
    i
    s1=step(s1,dt,false,'uncoupled',Re);
    %saveSol(s1,workDir,'sim',i)
    show(s1)
    %vtkCompleteOut(s1,workDirVtk,'field',i)
    %vort(s1);
end;
