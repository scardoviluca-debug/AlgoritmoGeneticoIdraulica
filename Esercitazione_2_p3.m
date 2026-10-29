clear
close all
clc

%% Variabili Globali
global NT NN q H0 A10 A12 A21 L ks N n
global costi Hmin Dprog
global Hk1 Hk2 Qk1 Qk2

%% Dati
Dprog = [75 100 150 200 250 300 350 400 500 600];
costi = [227 231 272 299 328 360 399 439 503 581];
nDprog = length(Dprog);
Hmin = 25;

%TRONCHI:
%colonna 1: numero del tronco
%colonna 2: nodo di uscita del tronco
%colonna 3: nodo di ingresso del tronco
%colonna 4: lunghezza in m
%colonna 5: diametro in mm
%colonna 6: coeff. di scabrezza del tronco in mm
%colonna 7: esponente della portata nella formula che esprime le perdite di carico
Tronchi = load('Tronchi.txt'); %carico la matrice dei tronchi

%NODI:
%colonna 1: numero del nodo
%colonna 2: 1 se nodo a carico imposto/noto; 0 altrimenti
%colonna 3: quota in m del nodo a carico imposto; 0 altrimenti
%colonna 4: domanda al nodo in l/s
%colonna 5: quota in m del nodo
Nodi = load('Nodi.txt'); %carico la matrice dei nodi

%Ricostruisco la matrice topologica della rete
NT = length(Tronchi); %numero dei tronchi
NN = length(Nodi); %numero totale dei nodi
AA = zeros(NT,NN); % matrice topologica di zeri

for i = 1:NT %ciclo che riempie la matrice topologica
     nodo_us = Tronchi(i,2);
     nodo_in = Tronchi(i,3);
     AA(i,nodo_us) = -1;
     AA(i,nodo_in) = 1;
 end

INS = find(Nodi(:,2)==1); %vettore degli indici dei nodi a carico imposto
IN = find(Nodi(:,2)==0); %vettore degli indici dei nodi a carico incognito

A10 = AA(:,INS); %matrice A10 dei nodi a carico noto
A12 = AA(:,IN); %matrice A12 dei nodi a carico incognito
A21 = A12'; %matrice A21

%Inizializzazione
Hk1 = ones(length(IN),1); %valori iniziali di prova del carico ai nodi
Qk1 = ones(NT,1); %valori iniziali di prova delle portate
Hk2 = ones(length(IN),1);
Qk2 = ones(NT,1);

q = Nodi(IN,4)*10^(-3); %vettore delle portate dei nodi a carico incognito [m^3/s]
H0 = Nodi(INS,3); %carico noto
L = Tronchi(:,4); %lunghezza delle condotte [m]
%D = Tronchi(:,5)*10^(-3); %diametri delle condotte [m]
ks = Tronchi(:,6)*10^(-3); % vettore delle scabrezze[m]
n = Tronchi(:,7); %vettore degli esponenti
N = diag(n); %matrice diagonale degli esponenti della portata

% Parametri algoritmo genetico
nind = 200;
ngen = 200;
npar = NT;
Pmut = 0.05; % probabilità del 5%

% Range parametri
par_min = ones(1,NT);
par_max = ones(1,NT)*nDprog;

% Popolazione iniziale
c = 0;
Chrom = zeros(nind,npar);
for i = 1:nind
     for j = 1:npar
         Chrom(i,j) = round(par_min(j)+rand*(par_max(j)-par_min(j)));        
     end
end

% Valuto ciascun individuo della popolazione iniziale
ValFO = zeros(1,nind);
for i = 1:nind
     ValFO(i) = Funzione_obiettivo(Chrom(i,:));
end
ChromFO = [Chrom,ValFO'];
bestFO = zeros(1,ngen);
for g = 1:ngen
    c = c+1;
    clc; disp(c)
    % associo a ciascun individuo fitness=1/FO
    FIT = 1./ChromFO(:,npar+1);
    ChromFO(:,npar+2) = FIT;

    % associo a ciascun individuo una probabilità P
    ChromFO = sortrows(ChromFO,npar+2); %riordina le righe in base alla colonna npar+2 (fitness) in senso crescente
    Prob = 2/nind*((1:nind)-1)/(nind-1); %fa si che l'area sottesa sia uguale a 1
    Probcum = cumsum(Prob);
    
    % genero la nuova popolazione
    ChromNEW = zeros(nind,npar);
    Figlio = zeros(2,npar);

    for i = 1:(nind/2)
        % seleziono una coppia di genitori
        x1 = rand;
        x2 = rand;
        posg1 = find(Probcum<x1, 1, 'last' );
        posg2 = find(Probcum<x2, 1, 'last' );
        while posg1 == posg2
            % i genitori devo essere distinti, altrimenti lo scambio di
            % geni non avviene
            x2 = rand;
            posg2 = find(Probcum<x2, 1, 'last' );
        end
        g1 = ChromFO(posg1,1:npar);
        g2 = ChromFO(posg2,1:npar);

        % Cross-over a pettine
        pettine = round(rand(1,npar));
        Figlio(1,:) = g1.*pettine+g2.*(1-pettine);
        Figlio(2,:) = g2.*pettine+g1.*(1-pettine);

        % mutazione
        for k = 1:2 % ciclo sui due figli
            mut = rand;
            if mut<=Pmut
                gene_mut = round(1+rand*(npar-1));
                Figlio(k,gene_mut) = round(par_min(gene_mut)+rand*(par_max(gene_mut)-par_min(gene_mut)));
            end
        end

        ChromNEW(i,:) = Figlio(1,:);
        ChromNEW(nind/2+i,:) = Figlio(2,:);

       
    end % chiusura ciclo generazione nuova popolazione


    % valuto ciascun individuo della nuova popolazione
    for i = 1:nind
        ValFO(i) = Funzione_obiettivo(ChromNEW(i,:));
    end
    ChromNEWFO = [ChromNEW,ValFO'];
    
    % elitismo
    if min(ChromNEWFO(:,npar+1))>min(ChromFO(:,npar+1))
        [val1,pos_pegg] = max(ChromNEWFO(:,npar+1));
        [val2,pos_migl] = min(ChromFO(:,npar+1));
 
        ChromNEWFO(pos_pegg,1:npar+1) = ChromFO(pos_migl,1:npar+1);
    end
    ChromFO = ChromNEWFO;
    bestFO(g) = min(ChromNEWFO(:,npar+1));
end % chiusura ciclo sulle generazioni


[val,pos_migl] = min(ChromFO(:,npar+1));
par_ottimale = ChromFO(pos_migl,1:npar);
 
% simulazione di verifica
CostoBest = Funzione_obiettivo(par_ottimale);

figure(1)
x = 1:1:ngen;
plot(x,bestFO)
grid on
xlabel('Generazioni')
ylabel('Costo [euro]')
    
    
    
    










