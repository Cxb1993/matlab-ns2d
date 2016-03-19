function n = vtkMeshOut(m,dir,name)
%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%MODEL2D model class constructor.
%   m = Model2d(m) creates a simulator object from the mesh object

%Name: vtkMeshOut
%Location: <path>/@Model2d
%Purpose: save mesh in VTK file format                         

% modificado em 20/05/2007
% revisado   em 20/05/2007

IEN = m.IEN;
X= m.X;
Y=m.Y;
Z=m.Z;
pc=m.pc;

nvert=size(pc,1);
nelem=size(IEN,1);
nnodes=nvert+nelem;

for i=1:size(IEN,1)
	for j=1:3
		IEN(i,j)=IEN(i,j)-1;
	end;
end;

fname = sprintf('%s%s.vtk',dir,name);
fid = fopen(fname, 'wt');

fprintf(fid, '# vtk DataFile Version 1.0\n');
fprintf(fid, 'Mesh Navier-Stokes 2d\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid, 'POINTS %d float\n',nvert);

for k=1:nvert
    fprintf(fid, '%2.5f %2.5f %2.2f\n', X(k), Y(k), 0);
end;

fprintf(fid, '\n',i);
fprintf(fid, 'CELLS %d %d\n',nelem,4*nelem);
for i=1:nelem
    fprintf(fid, '3 ');
    for j=1:3
        fprintf(fid, '%d ', IEN(i,j));
    end;
    fprintf(fid, '\n',i);
end;

fprintf(fid, '\n',i);
fprintf(fid, 'CELL_TYPES %d\n',nelem);
for i=1:nelem
    fprintf(fid, '5 ');
end;
fprintf(fid, '\n',i);

fclose(fid);

