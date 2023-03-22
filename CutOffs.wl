(* ::Package:: *)

BeginPackage["CutOffs`",{"BasicFunc`"}];

kmax;\[Omega]max;xmax;
kmin;\[Omega]min;xmin;
\[Omega]UH0;
(*\[CapitalDelta];*)
xOp;
leftCutOff::usage="leftCutOff[\[Omega]], return the left cut-off coordinate of trapped uh wave";
rightCutOff::usage="rightCutOff[\[Omega]], return the right cut-off coordinate of trapped uh wave";




Begin["`Private`"];
Needs["Profiles`"];
xmaxguess=8.5;
{kmax,\[Omega]max,xmax}={k,\[Omega],x}/.FindRoot[{Duh[k,\[Omega],x],Derivative[1,0,0][Duh][k,\[Omega],x],Derivative[0,0,1][Duh][k,\[Omega],x]},{{k,20},{\[Omega],f2\[Omega]@52.5 },{x,xmaxguess}}];
xminguess=2;
{kmin,\[Omega]min,xmin}={k,\[Omega],x}/.FindRoot[{Duh[k,\[Omega],x],Derivative[1,0,0][Duh][k,\[Omega],x],Derivative[0,0,1][Duh][k,\[Omega],x]},{{k,20},{\[Omega],f2\[Omega]@52.5 },{x,xminguess}}];

xOpguess=10;
\[CapitalDelta][\[Omega]_,x_]=\[CurlyEpsilon][\[Omega],x]^2-(((4 g[\[Omega],x]^2) \[Omega]^2) lt[\[Omega],x]^2)/c^2;
xOp=x/. FindRoot[\[CapitalDelta][\[Omega]min,x],{x,xOpguess}];

leftCutOff[\[Omega]_/;\[Omega]min<=\[Omega]<=\[Omega]max]:=$xl1[\[Omega]]=x/.FindRoot[\[CapitalDelta][\[Omega],x],{x,xmin,xmax},Method->"Brent"];
rightCutOff[\[Omega]_/;\[Omega]min<=\[Omega]<=\[Omega]max]:=$xr[\[Omega]]=x/.FindRoot[\[CapitalDelta][\[Omega],x],{x,xmax,xOp},Method->"Brent"];



(*SetAttributes[{leftCutOff,rightCutOff,rightRes,leftRes,\[Omega]UH0,kmax,\[Omega]max,xmax,kmin,\[Omega]min,xmin,xOp},{Protected}]*)


End[];
EndPackage[];
