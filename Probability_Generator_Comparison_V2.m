(* ::Package:: *)

(*\[Theta]sample=0.4;*)
Ma\[Mu]eV=10;
Plasma=1;
(*\[Alpha]rad=  0.00000000001 ;*)
\[Eta]grav=0;
MagneticOn=1;


(*Witte!! Pulsar details and axion-photon coupling*)
gGEV=10^-14;
Psec= 2\[Pi];
BGauss=10^14;


c = 2.99792 10^10 CM/SEC;
hbar = 6.582119 10^-16 ev SEC;


sec=(SEC/.Solve[hbar==1,SEC][[1]])/.{ev-> 10^-9 GEV}//Simplify


cm=CM/.Solve[c==1,CM][[1]]/.{SEC->sec}


(*====Physical Constants====*)
Gauss=1.95 10^-20 GEV^2;
\[Alpha]EM=1/137;
eEM = Sqrt[4\[Pi] \[Alpha]EM];
Me =0.5 10^-3 GEV;

KM = 10^5 cm;
\[Mu]eV=10^-6 10^-9 GEV;


va=220 *(10^5 cm/sec)


G=132712000000.0  (KM / Msun) * (KM/sec)^2;
MNS= 1.0 Msun;(*NS mass in solar masses*)


(*===\[Equal]CONVERSION TO DIMENSIONLESS VARIABLES IN UNITS R=1 for  10 km NS =======*)
me=Me(10KM);
\[CapitalOmega]R=((2\[Pi])/(Psec sec))(10KM);
maR=(Ma\[Mu]eV \[Mu]eV)(10KM);
\[Beta]=(10KM)^2 (BGauss Gauss);(*magnetic flux through a surface size RNS^2*)
RS = \[Eta]grav (2 G MNS)/(10KM);
gGEVR = gGEV GEV^-1 (10KM)^-1;
na = cm^-3 (10 KM)^3; (*dimensionless number density*)
\[Rho]= 0.45 GEV/cm^3 (10 KM)^4;


eVperSec = (10^-9 GEV /sec)(10 KM)^2 (*in dimensionless units of R = 10 KM*);


rC = Simplify[(((BGauss Gauss) (2\[Pi])/(Psec sec) \[Alpha]EM 4\[Pi])/((Ma\[Mu]eV \[Mu]eV)^2 (0.5 10^-3 GEV)))^(1/3),GEV>0];(*rc stellar units - see eq. (61)of 1910.11907*);


(*\[CapitalDelta]b0=0.02 rC;*)(*in units of R - used in first serach*)


X={tt,rr,\[Theta]\[Theta],\[CurlyPhi]\[CurlyPhi]};(*x index up*)
K={kt,kr,k\[Theta],k\[CurlyPhi]};(*k index down*)

(*Golreich Julian Density*)

nGJ[t_,r_,\[Theta]_,\[CurlyPhi]_,\[Alpha]_,\[CapitalOmega]_]=(\[Beta] \[CapitalOmega] (Cos[\[Alpha]]+3 Cos[\[Alpha]] Cos[2 \[Theta]]+3 Cos[\[CurlyPhi] (*-t \[CapitalOmega]*)] Sin[\[Alpha]] Sin[2 \[Theta]]))/(2 eEM r^3);


\[Omega]p[t_,r_,\[Theta]_,\[CurlyPhi]_,\[Alpha]_,\[CapitalOmega]_]=(   ( 4 \[Pi] \[Alpha]EM nGJ[t,r,\[Theta],\[CurlyPhi],\[Alpha],\[CapitalOmega]]/me  )^2 )^(1/4);


(*k= va maR;*)


KHatXYZ =  {Sin[\[Theta]v]Cos[\[CurlyPhi]v], Sin[\[Theta]v]Sin[\[CurlyPhi]v], Cos[\[Theta]v]};


xhat={1,0,0};
yhat={0,1,0};
zhat={0,0,1};


rhat = Sin[\[Theta]]Cos[\[CurlyPhi]] xhat + Sin[\[Theta]]Sin[\[CurlyPhi]]yhat + Cos[\[Theta]]zhat;
\[Theta]hat = Cos[\[Theta]]Cos[\[CurlyPhi]] xhat + Cos[\[Theta]]Sin[\[CurlyPhi]]yhat - Sin[\[Theta]]zhat;
\[CurlyPhi]hat =-Sin[\[CurlyPhi]] xhat + Cos[\[CurlyPhi]]yhat;


kHatSpher = {rhat.KHatXYZ,\[Theta]hat.KHatXYZ,\[CurlyPhi]hat.KHatXYZ};


Brc = \[Beta]/r^3 (Cos[\[Alpha]] Cos[\[Theta]]+ Sin[\[Alpha]] Sin[\[Theta]] Cos[\[CurlyPhi]-t \[CapitalOmega]] );
B\[Theta]c = \[Beta]/(2 r^3) (Cos[\[Alpha]] Sin[\[Theta]]- Sin[\[Alpha]] Cos[\[Theta]] Cos[\[CurlyPhi]-t \[CapitalOmega]] );(*Note leroy has a typo - see appendix of this NB*)
B\[CurlyPhi]c = \[Beta]/(2 r^3)(Sin[\[Alpha]] Sin[\[CurlyPhi]-t \[CapitalOmega]]);


BNS= {Brc,B\[Theta]c,B\[CurlyPhi]c};


Cos\[Theta]B = (BNS.kHatSpher)/((BNS.BNS)^(1/2))//Simplify;


klo =\[Omega] Sqrt[\[Omega]^2 -WP^2]/Sqrt[\[Omega]^2 - WP^2 Cos\[Theta]B^2];


E\[Gamma] = Sqrt[1/2 (k^2+ \[Omega]P^2+Sqrt[\[Omega]P^4+k^4+2 k^2 \[Omega]P^2 (1 - 2 Cos\[Theta]B^2)])]/.{\[Omega]P-> \[Omega]p[t,r,\[Theta],\[CurlyPhi],\[Alpha],\[CapitalOmega]]};

GradE=(Grad[E\[Gamma] ,{r,\[Theta],\[CurlyPhi]},"Spherical"].Grad[E\[Gamma] ,{r,\[Theta],\[CurlyPhi]},"Spherical"])^(1/2)/.{k-> klo};

E\[Gamma]Iso = Sqrt[k^2+ \[Omega]P^2]/.{\[Omega]P-> \[Omega]p[t,r,\[Theta],\[CurlyPhi],\[Alpha],\[CapitalOmega]]};

GradEIso=(Grad[E\[Gamma]Iso ,{r,\[Theta],\[CurlyPhi]},"Spherical"].Grad[E\[Gamma]Iso ,{r,\[Theta],\[CurlyPhi]},"Spherical"])^(1/2)/.{k-> klo};


Cos\[Theta]vp=Cos[VectorAngle[Grad[E\[Gamma] ,{r,\[Theta],\[CurlyPhi]},"Spherical"],kHatSpher]]/.{k-> klo};
Cos\[Theta]vpIso=Cos[VectorAngle[Grad[E\[Gamma]Iso ,{r,\[Theta],\[CurlyPhi]},"Spherical"],kHatSpher]]/.{k-> Sqrt[\[Omega]^2 -WP^2]};


VpAniso=(k/\[Omega])/.{k-> klo};
VpIso=(k/\[Omega])/.{k-> Sqrt[\[Omega]^2 -WP^2]};


\[Theta]m=0.4;


rc[\[Alpha]_,\[Theta]_,\[CurlyPhi]_,\[CapitalOmega]_,t_]=((2 \[Pi])^(1/3) (\[Alpha]EM \[Beta] \[CapitalOmega]R(Abs[Cos[\[Alpha]]+3 Cos[\[Alpha]] Cos[2 \[Theta]]+3 Cos[\[CurlyPhi]-t \[CapitalOmega]] Sin[\[Alpha]] Sin[2 \[Theta]]]))^(1/3))/(eEM^(1/3) (maR)^(2/3) me^(1/3));


WP=\[Omega]p[t,r,\[Theta],\[CurlyPhi],\[Alpha],\[CapitalOmega]];


\[Omega]inf = maR (1 + va^2)^(1/2);


P[\[Alpha]_,\[Theta]_,\[CurlyPhi]_,\[CapitalOmega]_,t_,\[Theta]v_,\[CurlyPhi]v_,\[Omega]_] = (\[Pi]/2)  gGEVR^2 (BNS.BNS) \[Omega]^4 (1 - Cos\[Theta]B^2) / (Cos\[Theta]B^2 WP^2 (WP^2 - 2 \[Omega]^2) + \[Omega]^4) / Abs[VpAniso Cos\[Theta]vp GradE]/.{r-> rc[\[Alpha],\[Theta],\[CurlyPhi],\[CapitalOmega],t]};


Piso[\[Alpha]_,\[Theta]_,\[CurlyPhi]_,\[CapitalOmega]_,t_,\[Theta]v_,\[CurlyPhi]v_,\[Omega]_] = (\[Pi]/2)  gGEVR^2 (BNS.BNS)  ((1 - Cos\[Theta]B^2)/ Abs[VpIso Cos\[Theta]vpIso GradEIso])(\[Omega]/WP)/.{r-> rc[\[Alpha],\[Theta],\[CurlyPhi],\[CapitalOmega],t]};


PReg[\[Alpha]_,\[Theta]_,\[CurlyPhi]_,\[CapitalOmega]_,t_,\[Theta]v_,\[CurlyPhi]v_,\[Omega]_] = (\[Pi]/2)  gGEVR^2 (BNS.BNS) \[Omega]^4 (1 - Cos\[Theta]B^2) / (Cos\[Theta]B^2 WP^2 (WP^2 - 2 \[Omega]^2) + \[Omega]^4) / Abs[VpAniso GradE]/.{r-> rc[\[Alpha],\[Theta],\[CurlyPhi],\[CapitalOmega],t]};


PisoReg[\[Alpha]_,\[Theta]_,\[CurlyPhi]_,\[CapitalOmega]_,t_,\[Theta]v_,\[CurlyPhi]v_,\[Omega]_] = (\[Pi]/2)  gGEVR^2 (BNS.BNS)  (1 - Cos\[Theta]B^2)/ Abs[VpIso GradEIso](\[Omega]/WP)/.{r-> rc[\[Alpha],\[Theta],\[CurlyPhi],\[CapitalOmega],t]};


(*fiducial values*)
(*NS polar coordinates*)
\[Theta]0=0.5;
\[CurlyPhi]0=0.6 \[Pi];

(*fiducial axion polar velocity angle*)
\[Theta]v0=\[Pi]/3 ;
\[CurlyPhi]v0Reg=0.2;


rc\[Theta][\[Theta]_]=Interpolation[Table[{\[Theta],Abs[rc[\[Theta]m,\[Theta],\[CurlyPhi]0,\[CapitalOmega]R,0]]},{\[Theta],0,\[Pi]/2,\[Pi]/2/10000}],\[Theta]];


Show[PolarPlot[Boole[rc[\[Theta]m,\[Theta]-\[Pi]/2,\[CurlyPhi]0,\[CapitalOmega]R,0]>1]rc[\[Theta]m,\[Theta]-\[Pi]/2,\[CurlyPhi]0,\[CapitalOmega]R,0],{\[Theta],-\[Pi],\[Pi]}],Graphics[Circle[{0,0},1]]]


(*values of \[Theta] where the surface hits the star*)
\[Theta]val={};
NDSolveValue[{F'[\[Theta]]==rc\[Theta]'[\[Theta]],F[0]==rc\[Theta][0],WhenEvent[F[\[Theta]]==1, \[Theta]val=Append[\[Theta]val,\[Theta]]]},F[\[Theta]],{\[Theta],0,\[Pi]/2}];


(*Flux angle*)
CosAngle[\[Alpha]_,\[Theta]_,\[CurlyPhi]_,\[CapitalOmega]_,t_,\[Theta]v_,\[CurlyPhi]v_,\[Omega]_] = (Cos\[Theta]vp) /.{r-> rc[\[Alpha],\[Theta],\[CurlyPhi],\[CapitalOmega],t]};
cos\[Theta]vFunc[\[Theta]v_]=Interpolation[Table[{\[Theta]v,CosAngle[\[Theta]m,\[Theta]0,\[CurlyPhi]0,\[CapitalOmega]R,0,\[Theta]v,\[CurlyPhi]v0Reg,\[Omega]inf]},{\[Theta]v,0,2\[Pi],0.01}],\[Theta]v];


(*values of the poles*)
\[Theta]valPole={};
NDSolveValue[{F'[\[Theta]v]==cos\[Theta]vFunc'[\[Theta]v],F[0]==CosAngle[\[Theta]m,\[Theta]0,\[CurlyPhi]0,\[CapitalOmega]R,0,0,\[CurlyPhi]v0Reg,\[Omega]inf],WhenEvent[F[\[Theta]v]==0, \[Theta]valPole=Append[\[Theta]valPole,\[Theta]v]]},F[\[Theta]v],{\[Theta]v,0, 0.999 \[Pi]}];


Plot[CosAngle[\[Theta]m,\[Theta]0,\[CurlyPhi]0,\[CapitalOmega]R,0,\[Theta]v,\[CurlyPhi]v0Reg,\[Omega]inf],{\[Theta]v,0,2\[Pi]}]


(*P[\[Alpha]_,\[Theta]_,\[CurlyPhi]_,\[CapitalOmega]_,t_,\[Theta]v_,\[CurlyPhi]v_,\[Omega]_] *)
OneDvs3D=LogPlot[{P[\[Theta]m,\[Theta]0,\[CurlyPhi]0,\[CapitalOmega]R,0,\[Theta]v,\[CurlyPhi]v0Reg,\[Omega]inf],PReg[\[Theta]m,\[Theta]0,\[CurlyPhi]0,\[CapitalOmega]R,0,\[Theta]v,\[CurlyPhi]v0Reg,\[Omega]inf]}, {\[Theta]v,0,\[Pi]},
Frame-> False,
ImageSize->500,
Axes-> False,
PlotStyle->{Darker[Green],Darker[Blue]},FrameLabel->{"\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(v\)]\)","\!\(\*SubscriptBox[\(P\), \(a\[Gamma]\)]\)"},
PlotRange->{{0,\[Pi]},{10^-6,10^-2}},
PlotRangePadding-> None,
GridLines-> {\[Theta]valPole,None},
GridLinesStyle ->  Directive[Darker[Gray],Thickness[0.005], Dashing[0.01]],PlotRangePadding-> False]
Export["/home/mcdonaldj/Dropbox/Pete_Bjorn/Probabilites_Following_Millar_Comments/AnisotropicPlot.pdf",OneDvs3D];


Table[i/4 \[Pi],{i,0,4}]


Table[(1.0 i)/4 \[Pi],{i,0,4}]


(*P[\[Alpha]_,\[Theta]_,\[CurlyPhi]_,\[CapitalOmega]_,t_,\[Theta]v_,\[CurlyPhi]v_,\[Omega]_] *)
OneDvs3D=LogPlot[{Piso[\[Theta]m,\[Theta]0,\[CurlyPhi]0,\[CapitalOmega]R,0,\[Theta]v,\[CurlyPhi]v0Reg,\[Omega]inf],PisoReg[\[Theta]m,\[Theta]0,\[CurlyPhi]0,\[CapitalOmega]R,0,\[Theta]v,\[CurlyPhi]v0Reg,\[Omega]inf]}, {\[Theta]v,0,\[Pi]},
Frame-> False,
ImageSize->500,
Axes-> False,
PlotStyle->{Darker[Green],Darker[Blue]},FrameLabel->{"\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(v\)]\)","\!\(\*SubscriptBox[\(P\), \(a\[Gamma]\)]\)"},
PlotRange->{{0,\[Pi]},{10^-6,10^-2}},
PlotRangePadding-> None,
GridLines-> {\[Theta]valPole,None},
GridLinesStyle ->  Directive[Darker[Gray],Thickness[0.005], Dashing[0.01]],PlotRangePadding-> False]
Export["/home/mcdonaldj/Dropbox/Pete_Bjorn/Probabilites_Following_Millar_Comments/IsotropicPlot.pdf",OneDvs3D];


(*P[\[Alpha]_,\[Theta]_,\[CurlyPhi]_,\[CapitalOmega]_,t_,\[Theta]v_,\[CurlyPhi]v_,\[Omega]_] *)
OneDvs3D=LogPlot[{P[\[Theta]m,\[Theta]0,\[CurlyPhi]0,\[CapitalOmega]R,0,\[Theta]v0,\[CurlyPhi]v,\[Omega]inf],Piso[\[Theta]m,\[Theta]0,\[CurlyPhi]0,\[CapitalOmega]R,0,\[Theta]v0,\[CurlyPhi]v,\[Omega]inf]}, {\[CurlyPhi]v,0,2\[Pi]},
Frame-> False,
ImageSize->500,
Axes-> False,
PlotStyle->{Darker[Green],Darker[Blue]},FrameLabel->{"\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(v\)]\)","\!\(\*SubscriptBox[\(P\), \(a\[Gamma]\)]\)"},PlotRange->{{0,2\[Pi]},{10^-10,10^-4}},
PlotRangePadding-> None,
GridLines-> {\[Theta]valPole,None},
GridLinesStyle ->  Directive[Darker[Gray],Thickness[0.005], Dashing[0.01]],PlotRangePadding-> False]
Export["/home/mcdonaldj/Dropbox/Pete_Bjorn/Integrated_Probability/OneDvs3D.pdf",OneDvs3D];


(*P[\[Alpha]_,\[Theta]_,\[CurlyPhi]_,\[CapitalOmega]_,t_,\[Theta]v_,\[CurlyPhi]v_,\[Omega]_] *)
ProbRatio=LogPlot[P[\[Theta]m,\[Theta]0,\[CurlyPhi]0,\[CapitalOmega]R,0,\[Theta]v0,\[CurlyPhi]v,\[Omega]inf]/Piso[\[Theta]m,\[Theta]0,\[CurlyPhi]0,\[CapitalOmega]R,0,\[Theta]v0,\[CurlyPhi]v,\[Omega]inf], {\[CurlyPhi]v,0,2\[Pi]},Frame-> False,ImageSize->500,
(*PlotLegends->{"3D","1D"}*)PlotStyle->Black,FrameLabel->{"\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(v\)]\)","\!\(\*SubscriptBox[\(P\), \(a\[Gamma]\)]\)"},PlotRange->{{0,2\[Pi]},{0.1,10^2}},Axes-> False]
Export["/home/mcdonaldj/Dropbox/Pete_Bjorn/Integrated_Probability/ProbRat.pdf",ProbRatio];


CosPlot=Plot[{CosAngle[\[Theta]m,\[Theta]0,\[CurlyPhi]0,\[CapitalOmega]R,0,\[Theta]v0,\[CurlyPhi]vi,\[Omega]inf],0},{\[CurlyPhi]vi,0,2\[Pi]},PlotRange->{{0,2\[Pi]},{-1,1}},
Frame-> False,
Axes-> False,
PlotStyle->{Black,{Gray,Thin}},
PlotRangePadding->None,
GridLines-> {\[Theta]valPole,None},
GridLinesStyle ->  Directive[Darker[Gray],Thickness[0.005], Dashing[0.01]],PlotRangePadding-> False]
Export["/home/mcdonaldj/Dropbox/Pete_Bjorn/Integrated_Probability/CosPlot.pdf",CosPlot];


(*SamplePoints=Table[(Boole[rc[\[Alpha],\[Theta],\[CurlyPhi],\[CapitalOmega],t]>1] Abs[GradE]/maR)/.{r-> rc[\[Alpha],\[Theta],\[CurlyPhi],\[CapitalOmega],t]}/.{
\[Alpha]-> RandomReal[{0,\[Pi]/2}],
\[Theta]-> RandomReal[{0,\[Pi]}],
 \[CurlyPhi]-> RandomReal[{0,2\[Pi]}], 
\[Theta]v-> RandomReal[{0,\[Pi]}],
\[CurlyPhi]v-> RandomReal[{0,\[Pi]}],
\[CapitalOmega]-> \[CapitalOmega]R, 
\[Omega]\[Rule] \[Omega]inf,
t-> 0},{i,1,5 10^3}];*)


(*GradEMin=NMinimize[{Abs[GradE]/maR/.{r-> rc[\[Alpha],\[Theta],\[CurlyPhi],\[CapitalOmega],t]}/.{\[CapitalOmega]-> \[CapitalOmega]R,t-> 0},rc[\[Alpha],\[Theta],\[CurlyPhi],\[CapitalOmega]R,t]>= 1},{\[Alpha],\[Theta],\[CurlyPhi],\[Theta]v,\[CurlyPhi]v,t}]*)


(*ListLogPlot[{SamplePoints,Table[{i,GradEMin[[1]]},{i,1,Length[SamplePoints]}]},
Frame-> True,
FrameLabel->{Style["Rand[\[Alpha],\[Theta],\[CurlyPhi],\[Theta]v,\[CurlyPhi]v,t]",16],Style["R |\[Del]E|/\!\(\*SubscriptBox[\(m\), \(a\)]\)",16]},
PlotStyle->{Black,Red},
PlotLegends-> {"Random Pts","Numerical Minimum"}]*)


(*ContourPlot[(Abs[GradE]/maR)/.{r-> rc[\[Alpha],\[Theta],\[CurlyPhi],\[CapitalOmega],t]}/.{\[Alpha]-> \[Theta]m,\[Theta]-> 0.5, \[CurlyPhi]-> 0.6\[Pi],\[CapitalOmega]-> \[CapitalOmega]R,t-> 0},{\[Theta]v,0,\[Pi]},{\[CurlyPhi]v,0,2\[Pi]},PlotLegends->Automatic,Mesh-> None]*)


(*GraphicsRow[{
ContourPlot[Log[10,P[\[Theta]m,0.5,0.6 \[Pi],\[CapitalOmega]R,0,\[Theta]v,\[CurlyPhi]v]],{\[Theta]v,0,\[Pi]},{\[CurlyPhi]v,0,2\[Pi]},PlotLegends->Automatic], 
ContourPlot[Log[10,Piso[\[Theta]m,0.5,0.6 \[Pi],\[CapitalOmega]R,0,\[Theta]v,\[CurlyPhi]v]],{\[Theta]v,0,\[Pi]},{\[CurlyPhi]v,0,2\[Pi]},PlotLegends->Automatic]
}]*)


(*P[\[Alpha]_,\[Theta]_,\[CurlyPhi]_,\[CapitalOmega]_,t_,\[Theta]v_,\[CurlyPhi]v_]*)


(*this time we take values of \[CurlyPhi] = \[Pi] *)


rc\[Theta][\[Theta]_]=Interpolation[Table[{\[Theta],Abs[rc[\[Theta]m,\[Theta],\[Pi],\[CapitalOmega]R,0]]},{\[Theta],0,\[Pi],\[Pi]/2/10000}],\[Theta]];


(*values of the zeros*)
\[Theta]val={};
F\[Theta]=NDSolveValue[{F'[\[Theta]]==rc\[Theta]'[\[Theta]],F[0]==rc\[Theta][0],WhenEvent[F[\[Theta]]==1, \[Theta]val=Append[\[Theta]val,\[Theta]]]},F[\[Theta]],{\[Theta],0,\[Pi]}];


\[Theta]val


Plot[rc\[Theta][\[Theta]],{\[Theta],0,\[Pi]},GridLines-> {\[Theta]val,{1}},Frame-> True]


(*P[\[Theta]m,\[Theta]0,\[CurlyPhi]0,\[CapitalOmega]R,t,\[Theta]v0,\[CurlyPhi]v,\[Omega]inf]*)


PInt[\[Theta]_,\[CurlyPhi]_]:=NIntegrate[Sin[\[Theta]v] PReg[\[Theta]m,\[Theta],\[CurlyPhi],\[CapitalOmega]R,0,\[Theta]v,\[CurlyPhi]v,\[Omega]inf](*Boole[Abs[rc\[Theta][\[Theta]]]\[GreaterEqual] 1]Sqrt[ rc\[Theta][\[Theta]]^2 Sin[\[Theta]]^2(rc\[Theta][\[Theta]]^2+(rc\[Theta]'[\[Theta]])^2)] )*),{\[CurlyPhi]v,0,2\[Pi]},{\[Theta]v,0,\[Pi]}
,Method->"AdaptiveMonteCarlo"]


scale=10^6;
\[Theta]tab1=Table[{\[Theta],scale PInt[\[Theta],\[Pi]]},{\[Theta],0.01,\[Theta]val[[1]],0.0025}];
\[Theta]tab2=Table[{\[Theta],scale PInt[\[Theta],\[Pi]]},{\[Theta], \[Theta]val[[2]],\[Theta]val[[3]],0.0025}];
\[Theta]tab3=Table[{\[Theta],scale PInt[\[Theta],\[Pi]]},{\[Theta], \[Theta]val[[4]], \[Pi],0.0025}];


PIntIso[\[Theta]_,\[CurlyPhi]_]:=NIntegrate[Sin[\[Theta]v] PisoReg[\[Theta]m,\[Theta],\[CurlyPhi],\[CapitalOmega]R,0,\[Theta]v,\[CurlyPhi]v,\[Omega]inf](*Boole[Abs[rc\[Theta][\[Theta]]]\[GreaterEqual] 1]Sqrt[ rc\[Theta][\[Theta]]^2 Sin[\[Theta]]^2(rc\[Theta][\[Theta]]^2+(rc\[Theta]'[\[Theta]])^2)] )/.{BGauss\[Rule] 10^14}*),{\[CurlyPhi]v,0,2\[Pi]},{\[Theta]v,0,\[Pi]}
,Method->"AdaptiveMonteCarlo"]


\[Theta]tab1Iso=Table[{\[Theta],scale PIntIso[\[Theta],\[Pi]]},{\[Theta],0.01,\[Theta]val[[1]],0.0025}];
\[Theta]tab2Iso=Table[{\[Theta],scale PIntIso[\[Theta],\[Pi]]},{\[Theta], \[Theta]val[[2]],\[Theta]val[[3]],0.0025}];
\[Theta]tab3Iso=Table[{\[Theta],scale PIntIso[\[Theta],\[Pi]]},{\[Theta], \[Theta]val[[4]], \[Pi],0.0025}];


Pmax=7;
Pmin=0;


ShadeRegion=Show[
Plot[1Pmax,{\[Theta],\[Theta]val[[1]],\[Theta]val[[2]]},
PlotStyle->Pink,
PlotRange-> {{0,\[Pi]},{Pmin,Pmax}},
GridLines->{\[Theta]val,None},
Filling-> Axis,
Axes-> False,
Frame-> None,
GridLinesStyle ->  Directive[Darker[Red],Thickness[0.005], Dashing[0.01]],PlotRangePadding-> False
],
Plot[1Pmax,{\[Theta],\[Theta]val[[3]],\[Theta]val[[4]]},
PlotStyle->Pink,
PlotRange-> {{0,\[Pi]},{Pmin,Pmax}},
GridLines->{\[Theta]val,None},
Filling-> Axis,
Axes-> False,
GridLinesStyle ->  Directive[Darker[Red],Thickness[0.005], Dashing[0.01]],PlotRangePadding-> False]
];


IntProb=ListPlot[{\[Theta]tab1,\[Theta]tab2,\[Theta]tab3},
PlotRange-> {Automatic,{Pmin,Pmax}},Joined-> True,PlotStyle->Darker[Green],Frame-> False,Filling-> Axis,Axes-> None,PlotRangePadding->None];


IntProbIso=ListPlot[{\[Theta]tab1Iso,\[Theta]tab2Iso,\[Theta]tab3Iso},
PlotRange-> {Automatic,{Pmin,Pmax}},Joined-> True,PlotStyle->Darker[Blue],Frame-> False,Filling-> Axis,Axes-> None,PlotRangePadding->None];


PowerPlot=Show[ShadeRegion,IntProb,IntProbIso,ImageSize-> 500,GridLines->{\[Theta]val,None},PlotRangePadding-> None,AspectRatio-> 1/3,
PlotRange-> {Automatic,{0,7}},PlotRangePadding-> False]


Export["/home/mcdonaldj/Dropbox/Pete_Bjorn/Integrated_Probability/PowerPlot.pdf",PowerPlot]
