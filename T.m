function T = T(theta)
    T = [cosd(theta)^2 sind(theta)^2 2*sind(theta)^2; sind(theta)^2 cosd(theta)^2 -2*sind(theta)*cosd(theta); -sind(theta)*cosd(theta) sind(theta)*cosd(theta) cosd(theta)^2-sind(theta)^2];
end