function [IEN_mesh,IENall,pAll] = elementIEN(ele,ien_map,solution_ien,p_ien)
%get the IEN array and P list for a given element
IEN_mesh = ien_map{ele,2};
IENall = solution_ien(ele,:);
IENall = IENall(IENall~=0);
pAll = p_ien(ele,:);
end

