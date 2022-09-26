function Qij = Qij(Q11, Q12, Q22, Q66, theta)
    Qr11 = Q11*cosd(theta)^4 + Q22*sind(theta)^4 + 2*(Q12+2*Q66)*sind(theta)^2*cosd(theta)^2;
    Qr12 = (Q11+Q22-4*Q66)*sind(theta)^2*cosd(theta)^2 + Q12*(cosd(theta)^4+sind(theta)^4);
    Qr22 = Q11*sind(theta)^4+Q22*cosd(theta)^4 + 2*(Q12+2*Q66)*sind(theta)^2*cosd(theta)^2;
    Qr16 = (Q11-Q12-2*Q66)*cosd(theta)^3*sind(theta)-(Q22-Q12-2*Q66)*sind(theta)^3*cosd(theta);
    Qr26 = (Q11-Q12-2*Q66)*cosd(theta)*sind(theta)^3-(Q22-Q12-2*Q66)*cosd(theta)^3*sind(theta);
    Qr66 = (Q11+Q22-2*Q12-2*Q66)*sind(theta)^2*cosd(theta)^2 + Q66*(sind(theta)^4+cosd(theta)^4);
    Qij = [Qr11 Qr12 Qr16;Qr12 Qr22 Qr26;Qr16 Qr26 Qr66];
end
