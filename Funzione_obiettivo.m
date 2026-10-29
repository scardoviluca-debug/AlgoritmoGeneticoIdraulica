function costTot = Funzione_obiettivo(xD)

  global q H0 A10 A12 A21 L ks N n
  global costi Hmin Dprog
  global Hk1 Hk2 Qk1 Qk2
  D(:,1) = Dprog(xD)/1000;
  costo = sum(costi(xD).*L');

  %metodo del gradiente
  lambda = (2.*log10((3.71.*D)./ks)).^(-2);
  bazin = (lambda.*8)./(9.81*pi^2); %coeff. di Bazin
  alpha = bazin.*(abs(Qk1).^(n-1)).*L./(D.^5);
  A11 = diag(alpha);
  G = N*A11;
  G1 = inv(G);
  M = inv(A21*G1*A12);
  Hk2 = M*(A21*Qk1-q)-(M*A21*G1)*(A11*Qk1+A10*H0);
  Qk2 = Qk1-G1*(A11*Qk1+A10*H0)-G1*A12*Hk2;

  %Procedura iterativa per il calcolo delle portate e dei carichi

    while max(abs(Hk2-Hk1)) >= 10^(-3)
      Qk1 = Qk2;
      Hk1 = Hk2;
      alpha = bazin.*(abs(Qk1).^(n-1)).*L./(D.^5);
      A11 = diag(alpha);
      G = N*A11;
      G1 = inv(G);
      M = inv(A21*G1*A12);
      Hk2 = M*(A21*Qk1-q)-(M*A21*G1)*(A11*Qk1+A10*H0);
      Qk2 = Qk1-G1*(A11*Qk1+A10*H0)-G1*A12*Hk2;
  end
  violazioni = find(Hk2<Hmin);
  violTot = sum(Hmin-Hk2(violazioni));
  
  costTot = costo*(1+violTot);
end

      