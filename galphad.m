(* ::Package:: *)

(* input parameters: m, dm, Rp, a *)
m = 10;
dmlist = Range[-5,0, 0.2];
a = 0.1;
denslist = Table[, {i,Length[dmlist]}];


(* general alpha-disk models *)
M = m msun ;
Kes = 0.4;
Ledd = 4 pi G M c/Kes;
dMedd = Ledd/(0.1 c^2);
(* dM = dm4 10^(-4) (msun/yr) *)

(* physical and math constants *)
pc = 3.08 10^(18);
msun = 2 10^(33);
mp = 1.67 10^(-24);
me = 9.11 10^(-28);
yr = 365.25 24 3600;
pi = 3.141592653589793;
k = 1.38 10^(-16);
s = 5.67 10^(-5);
A = 7.56 10^(-15);
G = 6.67 10^(-8);
c = 2.998 10^(10);


For [i=1, i<=Length[dmlist], i++, 

dm = 10^dmlist[[i]];
Print[dm];
dM = dm dMedd;

(* scalings *)
Rs = G M/c^2;
as = 1 ;
(* a = as ap *)		(* scaled alpha parameter *)
R = Rs Rp;

(* Planck mean, fit to HHe *)
KP = 10^(-4.05*Log[10,T] + Log[10,r] + 28.29);

(* Rosseland mean opacities from Bell & Lin *)

(* electron-scattering opacity *)
K0 = 0.348;
(* bb and ff opacity*)
K1 = 1.5 10^(20) r T^(-5/2);
(* H-scattering *)
K2 = 1.0 10^(-36) r^(1/3) T^(10);
(* Molecules *)
K3 = 10^(-8) r^(2/3) T^3;
(* evaporation of metal grains *)
K4 = 2 10^(81) r T^(-24);
(* metal grains *)
K5 = 0.1 T^(1/2);
(* Evaporation of ice grains *)
K6 = 2 10^(16) T^(-7);
(* ice grains *)
K7 = 2 10^(-4) T^2;

(* these relations are true in any opacity regime *)
W = (G M/R^3)^(1/2);
H = cs/W;
n = a cs H;  (* general version of nu*)
(* this commented version assumed gas pressure dominates *)
(* n = a k T/(mp W) *)
S = 2 H r;
t = S K/2;
F = (9/8) S n W^2;
Te = (F/s)^(1/4);

eq1 = T^4 / (3 t Te^4/8);
eq2 = n S / (dM/(3 pi));

(* first electron scattering, RP dom *)

cs := Sqrt[(4/9) A T^4/r];
mu := 0.6 mp;
K := K0;
	(* solve *)
eq1p = PowerExpand[eq1^(1/Exponent[eq1,T])];
T0 = PowerExpand[T /. Solve[eq1p == 1,T][[-1]]];
eq2p = PowerExpand[eq2 /. T -> T0];
eq2p = PowerExpand[eq2p^(1/Exponent[eq2p,r])];
r0 = PowerExpand[r /. Solve[eq2p == 1,r][[-1]]];
T0 = T0 /. r -> r0;
S0 = PowerExpand[S /. {r -> r0, T -> T0}];
cs0 = PowerExpand[cs /. {r -> r0, T -> T0}];
H0 = PowerExpand[H /. {r -> r0, T -> T0}];
F0 = PowerExpand[F /. {r -> r0, T -> T0}];
Q0 = cs0 W/(pi G S0);
betar0 = 3 r0 k T0/(mu A T0^4);
Print["zone 0 :"];Print[S0];Print[H0];

	(* outer boundary *)
r01a = PowerExpand[Rp /. Solve[betar0 == 1,Rp][[-1]]];
r01b = PowerExpand[Rp /. Solve[(K0 == K1) /. {r -> r0, T -> T0},Rp][[-1]]];
Print[r01a];
Print[r01b];

(* bf,ff , RP dom *)

cs := Sqrt[(4/9) A T^4/r];
mu := 0.6 mp;
K := K1;
	(* solve *)
eq1p = PowerExpand[eq1^(1/Exponent[eq1,T])];
T1 = PowerExpand[T /. Solve[eq1p == 1,T][[-1]]];
eq2p = PowerExpand[eq2 /. T -> T1];
eq2p = PowerExpand[eq2p^(1/Exponent[eq2p,r])];
r1 = PowerExpand[r /. Solve[eq2p == 1,r][[-1]]];
T1 = T1 /. r -> r1;
S1 = PowerExpand[S /. {r -> r1, T -> T1}];
cs1 = PowerExpand[cs /. {r -> r1, T -> T1}];
H1 = PowerExpand[H /. {r -> r1, T -> T1}];
F1 = PowerExpand[F /. {r -> r1, T -> T1}];
tau1 = PowerExpand[ S1 K1 /. {r -> r1, T -> T1}];
Q1 = cs1 W/(pi G S1);
betar1 = 3 r1 k T1/(mu A T1^4);
Print["zone 1 :"];Print[S1];Print[H1];

	(* outer boundary *)
r12a = PowerExpand[Rp /. Solve[betar1 == 1,Rp][[-1]]];
r12b = PowerExpand[Rp /. Solve[PowerExpand[K1/K2 /. 
	{r -> r1, T -> T1}] == 1,Rp][[-1]]]; 
Print[r12a]; 
Print[r12b]; 

(* bf,ff , GP dom *)

cs := Sqrt[k T/mu];
mu := 0.6 mp;
K := K1;
	(* solve *)
eq1p = PowerExpand[eq1^(1/Exponent[eq1,T])];
T2 = PowerExpand[T /. Solve[eq1p == 1,T][[-1]]];
eq2p = PowerExpand[eq2 /. T -> T2];
eq2p = PowerExpand[eq2p^(1/Exponent[eq2p,r])];
r2 = PowerExpand[r /. Solve[eq2p == 1,r][[-1]]];

T2 = T2 /. r -> r2;
S2 = PowerExpand[S /. {r -> r2, T -> T2}];
cs2 = PowerExpand[cs /. {r -> r2, T -> T2}];
H2 = PowerExpand[H /. {r -> r2, T -> T2}];
F2 = PowerExpand[F /. {r -> r2, T -> T2}];
tau2 = PowerExpand[ S2 K1 /. {r -> r2, T -> T2}];
Q2 = cs2 W/(pi G S2);
betar2 = 3 r2 k T2/(mu A T2^4);
Print["zone 2 :"];Print[r2];Print[T2];

	(* outer boundary *)

r23 = PowerExpand[Rp /. Solve[PowerExpand[K1/K2 /. 
	{r -> r2, T -> T2}] == 1,Rp][[-1]]];
r23a = PowerExpand[Rp /. Solve[betar2 == 1,Rp][[-1]]];	
Print[r23];
Print[r23a];

(* Hscat , GP dom *)

cs := Sqrt[k T/mu];
mu := mp;
K := K2;
	(* solve *)
eq1p = PowerExpand[eq1^(1/Exponent[eq1,T])];
T3 = PowerExpand[T /. Solve[eq1p == 1,T][[-1]]];
eq2p = PowerExpand[eq2 /. T -> T3];
eq2p = PowerExpand[eq2p^(1/Exponent[eq2p,r])];
r3 = PowerExpand[r /. Solve[eq2p == 1,r][[-1]]];

T3 = T3 /. r -> r3;
S3 = PowerExpand[S /. {r -> r3, T -> T3}];
cs3 = PowerExpand[cs /. {r -> r3, T -> T3}];
H3 = PowerExpand[H /. {r -> r3, T -> T3}];
F3 = PowerExpand[F /. {r -> r3, T -> T3}];
tau3 = PowerExpand[ S3 K2 /. {r -> r3, T -> T3}];
Q3 = cs3 W/(pi G S3);
betar3 = 3 r3 k T3/(mu A T3^4);
Print["zone 3:"];Print[r3];Print[T3];

	(* outer boundary *)
tmp = PowerExpand[K2/K3 /.  {r -> r3, T -> T3}];
r34 = PowerExpand[Rp /. Solve[
	PowerExpand[tmp^(1/Exponent[tmp,Rp])] == 1, Rp][[-1]]];
	
Print[r34];

(* molec , GP dom *)

cs := Sqrt[k T/mu];
mu := mp;
K := K3;
	(* solve *)
eq1p = PowerExpand[eq1^(1/Exponent[eq1,T])];
T4 = PowerExpand[T /. Solve[eq1p == 1,T][[-1]]];
eq2p = PowerExpand[eq2 /. T -> T4];
eq2p = PowerExpand[eq2p^(1/Exponent[eq2p,r])];
r4 = PowerExpand[r /. Solve[eq2p == 1,r][[-1]]];

T4 = T4 /. r -> r4;
S4 = PowerExpand[S /. {r -> r4, T -> T4}];
cs4 = PowerExpand[cs /. {r -> r4, T -> T4}];
H4 = PowerExpand[H /. {r -> r4, T -> T4}];
F4 = PowerExpand[F /. {r -> r4, T -> T4}];
tau4 = PowerExpand[ S4 K3 /. {r -> r4, T -> T4}];
Q4 = cs4 W/(pi G S4);
betar4 = 3 r4 k T4/(mu A T4^4);
Print["zone 4:"];Print[r4];Print[T4];

	(* outer boundary *)
tmp = PowerExpand[K3/K4 /.  {r -> r4, T -> T4}];
r45 = PowerExpand[Rp /. Solve[
	PowerExpand[tmp^(1/Exponent[tmp,Rp])] == 1, Rp][[-1]]];
Print[r45];

(* metal grain evap, GP dom *)

cs := Sqrt[k T/mu];
mu := mp;
K := K4;
	(* solve *)
eq1p = PowerExpand[eq1^(1/Exponent[eq1,T])];
T5 = PowerExpand[T /. Solve[eq1p == 1,T][[-1]]];
eq2p = PowerExpand[eq2 /. T -> T5];
eq2p = PowerExpand[eq2p^(1/Exponent[eq2p,r])];
r5 = PowerExpand[r /. Solve[eq2p == 1,r][[-1]]];

T5 = T5 /. r -> r5;
S5 = PowerExpand[S /. {r -> r5, T -> T5}];
cs5 = PowerExpand[cs /. {r -> r5, T -> T5}];
H5 = PowerExpand[H /. {r -> r5, T -> T5}];
F5 = PowerExpand[F /. {r -> r5, T -> T5}];
tau5 = PowerExpand[ S5 K4 /. {r -> r5, T -> T5}];
Q5 = cs5 W/(pi G S5);
betar5 = 3 r5 k T5/(mu A T5^4);
Print["zone 5:"];Print[r5];Print[T5];

	(* outer boundary *)
tmp = PowerExpand[K4/K5 /.  {r -> r5, T -> T5}];
r56 = PowerExpand[Rp /. Solve[
	PowerExpand[tmp^(1/Exponent[tmp,Rp])] == 1, Rp][[-1]]];
	
Print[r56];

(* metal grain, GP dom *)

cs := Sqrt[k T/mu];
mu := mp;
K := K5;
	(* solve *)
eq1p = PowerExpand[eq1^(1/Exponent[eq1,T])];
T6 = PowerExpand[T /. Solve[eq1p == 1,T][[-1]]];
eq2p = PowerExpand[eq2 /. T -> T6];
eq2p = PowerExpand[eq2p^(1/Exponent[eq2p,r])];
r6 = PowerExpand[r /. Solve[eq2p == 1,r][[-1]]];

T6 = T6 /. r -> r6;
S6 = PowerExpand[S /. {r -> r6, T -> T6}];
cs6 = PowerExpand[cs /. {r -> r6, T -> T6}];
H6 = PowerExpand[H /. {r -> r6, T -> T6}];
F6 = PowerExpand[F /. {r -> r6, T -> T6}];
tau6 = PowerExpand[ S6 K5 /. {r -> r6, T -> T6}];
Q6 = cs6 W/(pi G S6);
betar6 = 3 r6 k T6/(mu A T6^4);
Print["zone 6:"];Print[r6];Print[T6];

	(* outer boundary *)
tmp = PowerExpand[K5/K6 /.  {r -> r6, T -> T6}];
r67 = PowerExpand[Rp /. Solve[
	PowerExpand[tmp^(1/Exponent[tmp,Rp])] == 1, Rp][[-1]]];
	
Print[r67];

Rp = 10;

If [Rp < r01b, dens = a S0,
 If [r01b < Rp < r12a, dens= a S1,
   If [ r12a < Rp < r23, dens= a S2, 
     If [ r23 < Rp < r34, dens= a S3, 
       If [ r34 < Rp < r45, dens = a S4,
         If [r45 < Rp < r56, dens = a S5,
           If [r56 < Rp < r67, dens = a S6, 0]
 ]]]]]] ;
Print[dens];
(*AppendTo[denslist, dens]; *)
denslist[[i]] = dens;

(*ClearAll["Global`*"]; *)
ClearAll[Rp, dm, dM];
ClearAll[cs, mu, K, r, T]; 

] 


denslist
dmlist


(* ::Section:: *)
(*plots*)


(*ListLinePlot[Thread[{dmlist, denslist}], ScalingFunctions \[Rule] "Log"]*)
data = Transpose[{denslist, dmlist}];
ListLogLinearPlot[data,Joined -> True, PlotRange -> {{1, 10^4.5}, {-5, 0}}, FrameLabel -> {"\[Alpha] \[CapitalSigma]", "log" OverDot["m"]}, Frame -> True, ImageSize -> Large]


(* ::InheritFromParent:: *)
(*Graphics[{}, ContentSelectable -> True, ImageSize -> {480, 360}, PlotRange -> {{0, 480/360}, {0, 1}}]*)
