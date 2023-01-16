(* ::Package:: *)

BeginPackage["UHEigenMods`",{"CutOffs`"}];
EigenModsND;


Begin["`Private`"];
Needs["BasicFunc`"];


Int[\[Omega]_/;\[Omega]min<\[Omega]<=\[Omega]max]:=intmem[\[Omega]]=
NIntegrate[qplus[\[Omega],\[Xi]]-qminus[\[Omega],\[Xi]],{\[Xi],leftCutOff[\[Omega]],rightCutOff[\[Omega]]},Method->{"GlobalAdaptive",Method->"ClenshawCurtisRule"}];
Int[\[Omega]min]=NIntegrate[qplus[\[Omega]min,\[Xi]]-qminus[\[Omega]min,\[Xi]],{\[Xi],xmin,xOp},Method->{"GlobalAdaptive",Method->"ClenshawCurtisRule"}];
maxm=Floor[Int[\[Omega]min]/(2 \[Pi])-1/2];
{fmax,fmin}=\[Omega]2f/@{\[Omega]max,\[Omega]min};
FindModeND[m_]:=f/.FindRoot[{Int[f2\[Omega][f]]==2 \[Pi] m+\[Pi]},{f,fmin,0.9999999 fmax},Method->"Brent"];
EigenModsND[umstart_,umend_]:=With[{mstart=Min[Max[0,umstart],maxm],mend=Max[0,Min[umend,maxm]]},
Datafm=Table[{m, FindModeND[m]},{m,mstart,mend}]s
]

EigenModsND[]:=EigenModsND[0,maxm];


With[{coef=((Derivative[0,1,0][Duh][kmin,\[Omega]min,xmin])(-Derivative[2,0,0][Duh][kmin,\[Omega]min,xmin]Derivative[0,0,2][Duh][kmin,\[Omega]min,xmin]/4 )^(-1/2))^(1/2)},
\[Xi]0[\[Omega]_]:=Sqrt[(\[Omega]-\[Omega]min)]coef;]

Int[\[Omega]_/;\[Omega]min<=\[Omega]<=\[Omega]max,\[Gamma]_]:=int[\[Omega]]=NIntegrate[N[(qplus[\[Omega]+I \[Gamma],\[Xi]]-qminus[\[Omega]+I\[Gamma],\[Xi]])],{\[Xi],leftCutOff[\[Omega]],rightCutOff[\[Omega]]},Method->{"GlobalAdaptive",Method->"ClenshawCurtisRule"}];

rhsRe[\[Omega]_/;\[Omega]min<=\[Omega]<=\[Omega]max,m_]:=\[Pi] (2 m+1+\[Xi]0[\[Omega]]^2/(2 \[Pi]) Log[\[Xi]0[\[Omega]]^2/(2E)] -1/\[Pi] Arg[Gamma[1/2+I \[Xi]0[\[Omega]]^2/2]]);
zeroRe[\[Omega]_/;\[Omega]min<=\[Omega]<=\[Omega]max,\[Gamma]_,m_]:=Re[Int[\[Omega],\[Gamma]]]-rhsRe[\[Omega],m];






End[];
EndPackage[];
