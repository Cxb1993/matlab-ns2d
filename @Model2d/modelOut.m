function status = modelOut(m,fileName)

%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%MODEL2D model class constructor.
%   m = Model2d() creates a mesho object

%Name: test
%Location: <path>/@Model2d
%Purpose: model method to run test problems

% modificado em 01/05/2007
% revisado   em 09/04/2007


idbcu = m.idbcu;
idbcv = m.idbcv;
idbcp = m.idbcp;
uc = m.uc;
vc = m.vc;
pc = m.pc;

icc=[idbcu idbcv idbcp];


icc=unique(icc);

nbc=size(icc,2);

fname = sprintf('%s.bc',fileName);
fid = fopen(fname, 'wt');

fprintf(fid, 'u v p\n');
fprintf(fid, 'BC_DATA 3 %d\n',nbc);

for k=1:nbc
    ut=2; vt=2;pt=2;
    if(length(find(idbcu==icc(k)))>0)
        ut=1;
    end;
    if(length(find(idbcv==icc(k)))>0)
        vt=1;
    end;
    if(length(find(idbcp==icc(k)))>0)
        pt=1;
    end;

    fprintf(fid, '%d %d %2.2f %d %2.2f %d %2.2f\n', icc(k)-1, ut, uc(icc(k)), vt, vc(icc(k)), pt, pc(icc(k)));
end;
fclose(fid);

status=1
