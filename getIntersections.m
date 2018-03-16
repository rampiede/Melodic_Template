function [ tInt ] = getIntersections( f1, f2 )
    [tInt,~,~,~] = intersections(f1(:,1),f1(:,2),f2(:,1),f2(:,2),1);
end

