(* ::Package:: *)

BeginPackage["IonDispCurves`",{"HotDisp`"}];
IonDispCurves::usage="IDC=IonDispCurves[\[Omega],xStart,xEnd] 
use it like this: IDC\[LeftDoubleBracket]branchIndex,rangeIndex\[RightDoubleBracket][#]&
Returns nested list of dispersion curves for different branches and ranges between ion cyclotron resonance (\[Omega]=\!\(\*SubscriptBox[\(n\[Omega]\), \(ci\)]\)). 
";
IntPiecewise;
\[Omega]ciLines::usage="\[Omega]ciLines[\[Omega]_,xL_,xR_]";
IonDispCurvesPlot::usage="IonDispCurvesPlot[IDC]";


Begin["Private`"];
Needs["BasicFunc`"]

IonDispCurves[\[Omega]_,xStart_,xEnd_]:=With[{nSteps=50},
 
 (*\:041e\:043f\:0440\:0435\:0434\:0435\:043b\:044f\:0435\:043c \:0434\:0438\:0430\:043f\:0430\:0437\:043e\:043d\:044b*)
 xRes=NSolve[{\[Omega]==n \[Omega]ci[x],xStart<x<xEnd,1<n<100,n\[Element]Integers},{x,n}];
 rngFlat=Join[{xStart},Replace[x,xRes],{xEnd}];
 rngL=Transpose[{rngFlat[[;;-2]],rngFlat[[2;;]]}];
 (*rngL - list of ranges *)
 
 (*\:0418\:0449\:0435\:043c \:0442\:043e\:0447\:043a\:0443 \:0442\:0440\:0430\:043d\:0441\:0444\:043e\:0440\:043c\:0430\:0446\:0438\:0438 \:0432 \:043a\:0430\:0436\:0434\:043e\:043c \:0434\:0438\:0430\:043f*)
 {xtrL,qtrL} = Transpose@Table[{x,q} /. FindRoot[{DIBq[q, \[Omega], x], \[CurlyEpsilon]Exact[q, \[Omega], x]}, {q, 60}, {x, Mean@rng}], {rng, rngL}];
 
 
 (*\:0418\:0449\:0435\:043c \:0434\:0438\:0441\:043f \:043a\:0440\:0438\:0432\:0443\:044e \:0432 \:043a\:0430\:0436\:0434\:043e\:043c \:0434\:0438\:0430\:043f*)

qData[x_,ql_,qr_]:=Block[{sol},sol=FindRoot[\[CurlyEpsilon]11Hot[q,\[Omega],x],{q,ql,qr},Method->"Brent"];If[Abs[\[CurlyEpsilon]11Hot[q,\[Omega],x]/.sol]<10^-3,{x,q/.sol},Nothing,Nothing]];


 qibL=qlhL=Table[#&,{ri,Length@rngL}];

Table[
rng=Table[1.001rngL[[ri,1]]+(Min[xtrL[[ri]],rngL[[ri,2]]]-1.001rngL[[ri,1]])Log[nSteps,i],{i,1,nSteps-1}];
qibL[[ri]]=Interpolation[Join[
Table[qData[x,qtrL[[ri]],1000],{x,rng}],
{If[xtrL[[ri]]<rngL[[ri,2]],{xtrL[[ri]],qtrL[[ri]]},Nothing]}]];
qlhL[[ri]]=Interpolation[Join[
Table[qData[x,1,qtrL[[ri]]],{x,rng}],
{If[xtrL[[ri]]<rngL[[ri,2]],{xtrL[[ri]],qtrL[[ri]]},Nothing]}]];,
{ri,Length@rngL}];


(*\:041d\:0430 \:0432\:044b\:0445\:043e\:0434 \:0438\:0434\:0451\:0442 \:0441\:043f\:0438\:0441\:043e\:043a \:0438\:0437 \:0441\:043f\:0438\:0441\:043a\:043e\:0432 InterpolatingFunction \:0434\:043b\:044f \:043a\:0430\:0436\:0434\:043e\:0439 \:043e\:0431\:043b\:0430\:0441\:0442\:0438*)
{qlhL,qibL}

 ]

IntPiecewise[IntList_]:=Piecewise[Table[{IntList[[ri]][x],IntList[[ri]][[1,1,1]]<x<IntList[[ri]][[1,1,2]]},{ri,Length[IntList]}],10^-8 I];

\[Omega]ciLines[\[Omega]_,xL_,xR_]:=Module[{xres=NSolve[{\[Omega]==n \[Omega]ci[x],xL<x<xR,1<n<100,n\[Element]Integers},{x,n}]},
Graphics[{Red,Dashed,InfiniteLine[{{x,0},{x,1}}]/.xres},Frame->True]
];

IonDispCurvesPlot[IonDispCurves_]:=Block[{qLH,qIB},
{qLH,qIB}=IntPiecewise/@IonDispCurves;
With[{xl=IonDispCurves[[1]][[1]][[1,1,1]],xr=IonDispCurves[[1]][[-1]][[1,1,2]]},
Plot[{qIB,qLH},{x,xl,xr},PlotRange->{Automatic,{Automatic,200}},PlotStyle->{Directive[RGBColor[0.368417, 0.506779, 0.709798]],Directive[RGBColor[0.560181, 0.691569, 0.194885]]}]
]
]


End[];
EndPackage[]
