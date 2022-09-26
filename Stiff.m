function Stiff = Stiff(Qij1, Qij2, Qij3, z)
    z0 = z(1);
    z1 = z(2);
    z2 = z(3);
    z3 = z(4);
    
    A = Qij1*(z1-z0) + Qij2*(z2-z1) + Qij3*(z3-z2);
    B = (1/2)*((Qij1*(z1^2-z0^2) + Qij2*(z2^2-z1^2) + Qij3*(z3^2-z2^2)));
    D = (1/3)*((Qij1*(z1^3-z0^3) + Qij2*(z2^3-z1^3) + Qij3*(z3^3-z2^3)));

    Stiff = [A B;B D];
end