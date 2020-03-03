function [normals, flag] = coneGen(normal, mu1, mu2, fdir1 = 0, nd = 4)
  
  if ~isvector(normal)
    flag = -1;
    printf("Error in input arguments");
  endif
  
  if norm(normal)~=1
    printf(["Warning! Input vector of the wanted surface are not normalized.\n",...
    "Check whether you put a right vector!\n"]);
  endif
  
  baseVec = [0; 0; 1];
  normal = normalizeVector(normal);
  normals = zeros(nd, 3);
  if ~isequal(normal, baseVec) 
    rot_angle = acos((baseVec')*normal);
    rot_axis = normalizeVector(cross(baseVec, normal));
    Transform = createRotation3dLineAngle(rot_axis, rot_angle);
    Transform = Transform(1:3, 1:3);
  else
    Transform = eye(3);
  endif
  
  if isempty(mu2)
    mu2 = mu1;
  endif
  
  Elipse = [0, 0, 1/mu1, 1/mu2, fdir1];
  [x, y] = ellipseAsPolygon(Elipse, nd);
  
  coneVec = normalizeVector([x, y, ones(nd+1,1)]); 
  
  for i=1:nd
    normals(i,:) = (Transform*cross(coneVec(i+1,:), coneVec(i,:))')';
    normals(i,:) = normalizeVector(normals(i,:));
  endfor
  
endfunction