pluck[{{f_, g_}, {umin_, umax_}}, (Nbordpts_)?
       (IntegerQ[#1] && Positive[#1] & ), {Ax_, Ay_}, Check_:True] /; 
     If[RegionMember[parametricregion[{{f, g}, {umin, umax}}], {Ax, Ay}], 
      True, Message[pluck::hitinregion, {Ax, Ay}]; Abort[]; ] := 
    Module[{hmax, \[Delta], \[Gamma], recurt, constructxy, oneline, datapts, 
      recurtpoints, mesh}, 
     If[Check, If[NMinValue[{Sqrt[(f[u] - Ax)^2 + (g[u] - Ay)^2], u > umin, 
           u <= umax}, u] <= 0.1, Message[pluck::nearedge, {Ax, Ay}], False], 
       False]; hmax = 1; \[Delta] = 0.05; \[Gamma][{x1_, y1_}, {x2_, y2_}][
        t_] := {(1 - t)*x1 + t*x2, (1 - t)*y1 + t*y2}; 
      recurt[{n_, t_, z_}] := {n + 1, 
        1 - ((1 - t)/(hmax + \[Delta]/(n + 1) - z))*((2*\[Delta])/(n + 1)), 
        hmax - \[Delta]/(n + 1)}; recurtpoints = 
       Join[NestList[recurt, {0, 0, 0}, 10], 
        Table[{0, (1 - (2*\[Delta])/(hmax + \[Delta]))*(i/4), 
          (hmax + \[Delta])*(1 - (2*\[Delta])/(hmax + \[Delta]))*(i/4)}, 
         {i, 1, 3}]]; constructxy[U_][{n_, t_, z_}] := 
       {\[Gamma][{f[U], g[U]}, {ax, ay}][t], -z}; oneline[Num_][U_] := 
       constructxy[U] /@ recurtpoints; datapts = 
       DeleteDuplicatesBy[Append[Flatten[oneline[10] /@ Range[umin, umax, 
             (umax - umin)/Nbordpts], 1], {{ax, ay}, -hmax}], First] /. 
        {ax -> Ax, ay -> Ay}; mesh = NDSolve`FEM`ToElementMesh[
        datapts[[All,1]] + RandomReal[{-10^(-9), 10^(-9)}, 
          {Length[datapts]}]]; NDSolve`FEM`ElementMeshInterpolation[{mesh}, 
       datapts[[All,2]], InterpolationOrder -> 1]]
 
pluck /: pluck::hitinregion = 
     "The point `1` is not inside the region defined."
 
pluck /: pluck::nearedge = "Input point `1` is very close to the border of \
the region. Expect interpolation to not be smooth."
 
pluck /: pluck::usage = "pluck[{{f,g},{umin,umax}},Nbordpts,{Ax,Ay},Check] \
produces an interpolation function that has a 'dimple' centered at the \
thump-point {Ax,Ay}\n-f,g are the *heads* of the border functions\n-Check is \
optional and decides whether to check if point is near border (default is \
True)"
 
f = 129.2652371397844
 
parametricregion[{{f_, g_}, {umin_, umax_}}] := Module[{pts, plot}, 
     plot = ParametricPlot[{f[u], g[u]}, {u, umin, umax}]; 
      pts = Flatten[Table[plot[[1,1,3,i,1]], 
         {i, 2, Length[Position[plot, Line]] + 1}], 1]; 
      pts = DeleteDuplicates[pts]; BoundaryMeshRegion[pts, 
       Line[Join[Range[Length[pts] - 1], {1}]]]]
 
parametricregion /: parametricregion::usage = "parametricregion[{{f,g},{umin,\
umax}}] produces a region from the input simple closed parametrized \
border\n-f,g are the *heads* of the border functions"
 
pts = {{0, 1}, {1, 0}, {-1, 0}, {0, -1}, {0.008000000000000007, 
     0.0050000000000001155}}
 
n = 5
 
peak[{{f_, g_}, {umin_, umax_}}, (\[Delta]_)?Positive, {x0_, y0_}, 
      check_:True] /; If[RegionMember[parametricregion[
        {{f, g}, {umin, umax}}], {x0, y0}], True, 
      Message[pluck::hitinregion, {x0, y0}]; Abort[]; ] := 
    Module[{A, \[Sigma]x, \[Sigma]y, Gaussian}, 
     If[check, If[NMinValue[{Sqrt[(f[u] - x0)^2 + (g[u] - y0)^2], u > umin, 
           u <= umax}, u] <= 2*\[Delta], Message[peak::nearedge, {x0, y0}], 
        False], False]; A = -1; \[Sigma]x = \[Delta]; \[Sigma]y = \[Delta]; 
      Gaussian[u_, v_] := A*Exp[-((u - x0)^2/(2*\[Sigma]x^2) + 
           (v - y0)^2/(2*\[Sigma]y^2))]; Gaussian]
 
peak /: peak::hitinregion = "The point `1` is not inside the region defined."
 
peak /: peak::nearedge = "Input point `1` is near the border of the region, \
so the produced function will not be equal to 0 on the border close to `1`."
 
peak /: peak::usage = "peak[{{f_,g_},{umin_,umax_}},\[Delta],{x0,y0},check] \
produces a spatial Gaussian with peak at {\!\(\*SubscriptBox[\(x\), \
\(0\)]\),\!\(\*SubscriptBox[\(y\), \(0\)]\)} and width \[Delta]\n-f,g are the \
*heads* of the border functions\n-check is an optional boolean used to see if \
the hitpoint {x0,y0} is near the border and the default is True"
 
egsystem[(Num_)?(IntegerQ[#1] && Positive[#1] & ), const_:1][
     (region_)?RegionQ] := Module[{f, freqs, vals, funcs}, 
     {vals, funcs} = NDEigensystem[{(-const)*Laplacian[f[x, y], {x, y}], 
         DirichletCondition[f[x, y] == 0, True]}, {f}, 
        Element[{x, y}, region], Num]; freqs = Sqrt[vals]; {freqs, funcs}]
 
egsystem /: egsystem::usage = "egsystem[Num,const][region] determines the \
eigensystem up to Num'th eigenvalue for the given region\n-const is the \
multiplied into Laplacian and is optional, default = 1"
 
QOtriangleApprox[efunclist_List] := Module[{emesh, indices}, 
     emesh = efunclist[[1]]["ElementMesh"]; 
      indices = emesh["MeshElements"] /. {NDSolve`FEM`TriangleElement[t_]} :> 
         Flatten[t[[All,4 ;; 6]]]; 
      Table[Total[Total[Partition[eigenfunction["ValuesOnGrid"][[indices]], 
           3], {2}]*(First[eigenfunction["ElementMesh"][
            "MeshElementMeasure"]]/3)], {eigenfunction, efunclist}]]
 
QOtriangleApprox[efunclist_List, hit_] := Module[{emesh, indices, hitvals}, 
     emesh = efunclist[[1]]["ElementMesh"]; 
      indices = emesh["MeshElements"] /. {NDSolve`FEM`TriangleElement[t_]} :> 
         Flatten[t[[All,4 ;; 6]]]; hitvals = (hit[#1[[1]], #1[[2]]] & ) /@ 
        efunclist[[1]]["Grid"][[indices]]; 
      Table[Total[Total[Partition[hitvals*eigenfunction["ValuesOnGrid"][[
             indices]], 3], {2}]*(First[eigenfunction["ElementMesh"][
            "MeshElementMeasure"]]/3)], {eigenfunction, efunclist}]]
 
QOtriangleApprox /: QOtriangleApprox::usage = "QOtriangleApprox[efunclist] is \
a VERY fast approximation that uses the Quadratic Order Triangle Rule \
(performs on each element on list \
input)\nQOtriangleApprox[efunclist,hitfunction] does the same fast \
approximation but has the integrand efunclist[[i]]*hit\n**Function list must \
be interpolation functions over the same mesh**\nWARNING: If interpolation \
functions' meshes are not made of triangles, this will give a false \
approximation."
 
freqsintens[{{f_, g_}, {umin_, umax_}}, (numberfreqs_)?
       (IntegerQ[#1] && Positive[#1] & ), method_:"Peak"][thumpoints_] := 
    Module[{region, funcs, freqs, thumps, coefs, intens}, 
     region = parametricregion[{{f, g}, {umin, umax}}]; 
      {freqs, funcs} = egsystem[numberfreqs][region]; 
      thumps = Switch[method, "Peak", Table[peak[{{f, g}, {umin, umax}}, 
          0.05, thumpoints[[i]]], {i, 1, Length[thumpoints]}], "Pluck", 
        Table[pluck[{{f, g}, {umin, umax}}, 100, thumpoints[[i]]], 
         {i, 1, Length[thumpoints]}], _, Message[freqsintens::fmethod, 
          OptionValue[method]]; Abort[]; ]; 
      coefs = Table[QOtriangleApprox[funcs, thumps[[i]]], 
        {i, Length[thumps]}]; intens = Table[(coefs[[i,All]]/coefs[[i,1]])^2, 
        {i, 1, Length[coefs]}]; Table[{freqs[[All]], intens[[i,All]]}, 
       {i, 1, Length[thumpoints]}]]
 
freqsintens /: freqsintens::fmethod = 
     "Method `1` is not a valid method. Valid methods are: {Peak,Pluck}"
 
freqsintens /: freqsintens::usage = "freqsintens[{{f,g},{umin,umax}},numberfr\
eqs,method][thumpoints] gives the first numberfreqs eigenfrequencies of the \
region {{f,g},{umin,umax}} with corresponding intensities when thumped at the \
points thumpoints\n-method is how the thumps are modeled and is optional, \
default = Peak\n-Uses parametricregion, egsystem, pluck or peak, and \
QOtriangleApprox2\n-Scales the intensities so that the fundamental has \
intensity 1\n-outputs a list in the \
form\n\t{{{freqlist,intensitylist1},{freqlists,intensitylist2},...},...}"
 
graphfi[{{f_, g_}, {umin_, umax_}}, (numberfreqs_)?
       (IntegerQ[#1] && Positive[#1] & ), method_:"Peak"][thumpoints_] := 
    Module[{list, fipairs}, list = freqsintens[{{f, g}, {umin, umax}}, 
         numberfreqs, method][thumpoints]; 
      fipairs = Table[{list[[i,1,x]], list[[i,2,x]]}, 
        {i, 1, Length[thumpoints]}, {x, 1, Length[list[[1,2]]]}]; 
      ListLogPlot[Table[Legended[fipairs[[i]], thumpoints[[i]]], 
        {i, 1, Length[thumpoints]}], Filling -> Axis, PlotRange -> All, 
       AxesLabel -> {"Frequency", "Intensity"}]]
 
graphfi /: graphfi::usage = "graphfi[{{f,g},{umin,umax}},numberfreqs,method][\
thumpts] outputs a frequency vs intensity graph for the given border and \
thump points\nUses a log scale on the y-axis, so it is easier to see all the \
frequencies\n**Uses freqsintens so will take just as long**"
 
keep[eigenfunctionlist_List, threshhold_:0] := Module[{ints, keeplist, i}, 
     ints = Abs[QOtriangleApprox[eigenfunctionlist]]; keeplist = {}; 
      For[i = 1, i <= Length[eigenfunctionlist], i++, 
       If[ints[[i]] >= threshhold, AppendTo[keeplist, True], 
        AppendTo[keeplist, False]]]; keeplist]
 
keep /: keep::usage = "keep[eigenfuncitonlist,threshhold] returns a list: \
True if Keep and False if Not Keep according to threshhold\n-threshhold is \
optional, default = 0\n-Requires eigenfunctions to be the head of function \
over region\n-Uses QOtriangleApprox"
 
loudfreqs[fun_List, threshhold_:0] := Module[{decision, keeps, i}, 
     decision = keep[fun, threshhold]; keeps = {}; 
      For[i = 1, i <= Length[fun], i++, If[decision[[i]], 
        AppendTo[keeps, i]; , False; ]]; keeps]
 
loudfreqs /: loudfreqs::usage = "loudfreqs[eigenfunctionlist, threshhold] \
returns list of indices of frequencies that do not 'cancel' themselves out \
according to threshhold\n-Uses keep"
 
graphweight[fun_List] := Module[{}, Clear[ints, lp, thp]; 
      ints = Abs[QOtriangleApprox[fun]]; 
      lp = ListPlot[ints, PlotRange -> All]; 
      thp[t_] := Plot[t, {x, 0, Length[ints] + 5}, Filling -> Axis, 
        FillingStyle -> Opacity[0.2], PlotStyle -> Thin]; 
      Manipulate[Show[lp, thp[threshhold]], {threshhold, 0, Max[ints], 
        (Max[ints] - Min[ints])/100}]]
 
graphweight /: graphweight::input = "[eigenfunctionlist]"
 
graphweight /: graphweight::usage = "graphweight[eigenfuncitonlist] produces \
a Manipulate[] command that can be used to find a good threshhold/cutoff for \
heard and un-heard frequencies\n-Uses QOtriangleApprox"
 
colorfreqs[fun_List, threshhold_:0] := Module[{decision, color}, 
     decision = keep[fun, threshhold]; Table[If[decision[[i]], Red, Blue], 
       {i, Length[fun]}]]
 
colorfreqs /: colorfreqs::usage = " colorfreqs[eigenfunctionlist,threshhold] \
assigns Red to indices that should be kept (given a threshhold) and Blue to \
indices that should not be kept\n-threshhold is optional, default = 0\n-Uses \
keep"
 
freqband[min_Numeric, max_Numeric, color_:RGBColor[1, 0, 0]] := 
    Module[{ave, thickness}, ave = Mean[{min, max}]; 
      thickness = (min - max)/2; {color, Thickness[thickness], Opacity[0.5], 
       InfiniteLine[{ave, 0}, {0, 1}]}]
 
freqband /: freqband::usage = "freqband[minimum, maximum, color] outputs a \
band for a single frequency in Graphics[]\n-color is optional, default = Red"
 
freqlines[minsystem_List, maxsystem_List, threshhold_:0] := 
    Module[{color, lines}, color = colorfreqs[minsystem[[2]], threshhold]; 
      lines[i_] := If[Length[maxsystem] == 0, {color[[i]], Thin, 
         Infinity[{minsystem[[1,i]]}, {0, 1}]}, freqband[minsystem[[1,i]], 
         maxsystem[[1,i]], color[[i]]]]; Graphics[Table[lines[i], 
        {i, Length[minsystem[[1]]]}], PlotRange -> 
        {{0, Max[minsystem[[1]]]*1.1}, {0, Max[minsystem[[1]]]*0.25}}, 
       Frame -> True]]
 
freqlines /: freqlines::usage = "freqlines[minsystem, maxsystem, threshhold] \
outputs graphic of frequency bands corresponding to minsystem and maxsystem \
with correct colors and thicknesses\n-threshhold is optional, default = \
0\n-Uses colorfreqs and freqband"
 
print3Dborder[{{x_, y_}, {umin_, umax_}}, (width_)?Positive, 
     (height_)?Positive] := ParametricPlot3D[{x[u], y[u], z}, 
     {u, umin, umax}, {z, 0, height}, PlotStyle -> Thickness[width], 
     PlotPoints -> 100, Boxed -> False, Mesh -> None, Axes -> False]
 
print3Dborder /: print3Dborder::usage = "print3Dborder[{{f,g},{umin,umax}},wi\
dth,height] produces a 3D border (rectangular cross-sections) that can be \
exported to .stl and 3D-printed\n-f,g need to be the heads of the border \
functions"
 
print3Dregion[{{f_, g_}, {umin_, umax_}}, (borderwidth_)?Positive, 
     (borderthickness_)?Positive, (membranethickness_)?Positive] := 
    Module[{region, membrane, border}, 
     region = parametricregion[{{f, g}, {umin, umax}}]; 
      membrane = ParametricPlot3D[{x, y, 0}, Element[{x, y}, region], 
        PlotStyle -> Thickness[membranethickness], Mesh -> None, 
        Boxed -> False, Axes -> False]; border = print3Dborder[
        {{f, g}, {umin, umax}}, borderwidth, borderthickness]; 
      Show[membrane, border]]
 
print3Dregion /: print3Dregion::usage = "print3Dregion[{{f,g},{umin,umax}},bo\
rderwidth,borderthickness,membranethickness] constructs the membrane and \
bottom border (possibly base) for a given parametrized border\n-f,g need to \
be the heads of the border functions\n-Uses parametricregion and \
print3Dborder"
 
p[V_List][t_] := Module[{\[ScriptL]}, 
     \[ScriptL][i_][t0_] := {(1 - t0)*V[[i,1]] + t0*V[[i + 1,1]], 
        (1 - t0)*V[[i,2]] + t0*V[[i + 1,2]]}; \[ScriptL][1][t] + 
       Sum[UnitStep[t - i + 1]*(\[ScriptL][i][t - i + 1] - 
          \[ScriptL][i - 1][t - i + 2]), {i, 2, Length[V] - 1}]]
 
p /: p::usage = "p[V][t] produces a parametrization of a polygonal line \
determined by points in list V in the order given"
 
vard[dmin_:0.1][\[CapitalOmega]_List] := Module[{\[CapitalOmega]star, i}, 
     \[CapitalOmega]star = {\[CapitalOmega][[2]] - \[CapitalOmega][[1]]}; 
      For[i = 2, i < Length[\[CapitalOmega]], i++, 
       AppendTo[\[CapitalOmega]star, \[CapitalOmega][[i + 1]] - 
         \[CapitalOmega][[i]]]]; For[i = 1, i <= Length[\[CapitalOmega]star], 
       i++, If[\[CapitalOmega]star[[i]] <= dmin, \[CapitalOmega]star = 
         Delete[\[CapitalOmega]star, i]]]; Variance[\[CapitalOmega]star]]
 
vard /: vard::usage = "vard[dmin][\[CapitalOmega]] outputs the variance in \
the differences of the input list of frequencies \[CapitalOmega]\n-dmin is \
optional, default = 0.1"
 
disvector[rmin_:1.01][\[CapitalOmega]_List] := 
    Module[{freqs, number, line, perf, displacement, i}, 
     freqs = {{1, \[CapitalOmega][[1]]}}; number = 1; 
      For[i = 2, i <= Length[\[CapitalOmega]], i++, 
       If[\[CapitalOmega][[i]]/\[CapitalOmega][[i - 1]] <= rmin, 
        AppendTo[freqs, {number, \[CapitalOmega][[i]]}], 
        AppendTo[freqs, {++number, \[CapitalOmega][[i]]}]]]; 
      line[t_] := Evaluate[Fit[freqs, {t}, t]]; 
      perf = line /@ freqs[[All,1]]; displacement = freqs[[All,2]] - perf; 
      displacement]
 
disvector /: disvector::usage = "disvector[rmin][\[CapitalOmega]] outputs a \
'displacement' vector where each element is the distance a frequency is from \
a corresponding perfect frequency\n-rmin is optional, default = 1.01"
 
soundanalysispg[filepath_, numhits_] := 
    Module[{bin, data, peak, peaks, peakfc, peaklocation, allpeaks}, 
     bin = Import[filepath, "Data"]; data = bin[[1,1 ;; -1 ;; 1]]; 
      peak = FindPeaks[data, 10, 0, 0.1]; peaks = Round[peak[[All,1]]]; 
      peakfc = FindClusters[peaks, numhits]; peaklocation = 
       Round[Mean /@ peakfc]; allpeaks = 
       Table[data[[peaklocation[[i]] - 300 ;; Min[{peaklocation[[i]] + 2000, 
            Length[data]}]]], {i, Length[peaklocation]}]; 
      Periodogram[allpeaks, SampleRate -> 11025, PlotStyle -> 
        {PlotTheme -> "Detailed"}, PlotRange -> {{0, 2000}, All}, 
       PlotLegends -> Table[IntegerString[i], {i, Length[allpeaks]}], 
       ImageSize -> Large]]
 
soundanalysispg /: soundanalysispg::usage = "soundanalysispg[filepath,numhits\
] produces the Periodogram of the sound file at filepath that has numhits \
hits in it"
 
comparetheovals[{theovals_, theocolor_}, pg_, amax_:1000] := 
    (ClearAll[vertlines]; vertlines[a_] := 
      Graphics[Table[{theocolor[[i]], Thin, InfiniteLine[{a*theovals[[i]], 
           0}, {0, 1}]}, {i, Length[theovals]}]]; 
     Manipulate[Show[pg, vertlines[a]], {a, 0, amax}])
 
comparetheovals /: comparetheovals::usage = "comparetheovals[{theovals,theoco\
lor},pg,amax] produces a Manipulate command that allows for lines at theovals \
(with each colored according to theocolor) to be scaled up to a factor of \
amax overlayed on pg."
 
vertlines[a$_] := Graphics[Table[{{RGBColor[1, 0, 0], RGBColor[0, 0, 1], 
         RGBColor[0, 0, 1], RGBColor[0, 0, 1], RGBColor[0, 0, 1], 
         RGBColor[1, 0, 0], RGBColor[0, 0, 1], RGBColor[0, 0, 1], 
         RGBColor[0, 0, 1], RGBColor[0, 0, 1], RGBColor[0, 0, 1], 
         RGBColor[0, 0, 1], RGBColor[0, 0, 1], RGBColor[0, 0, 1], 
         RGBColor[1, 0, 0], RGBColor[0, 0, 1], RGBColor[0, 0, 1], 
         RGBColor[0, 0, 1], RGBColor[0, 0, 1], RGBColor[0, 0, 1], 
         RGBColor[0, 0, 1], RGBColor[0, 0, 1], RGBColor[0, 0, 1], 
         RGBColor[0, 0, 1], RGBColor[0, 0, 1], RGBColor[0, 0, 1], 
         RGBColor[0, 0, 1], RGBColor[0, 0, 1], RGBColor[0, 0, 1], 
         RGBColor[1, 0, 0], RGBColor[0, 0, 1], RGBColor[0, 0, 1], 
         RGBColor[0, 0, 1], RGBColor[0, 0, 1], RGBColor[0, 0, 1], 
         RGBColor[0, 0, 1], RGBColor[0, 0, 1], RGBColor[0, 0, 1], 
         RGBColor[0, 0, 1], RGBColor[0, 0, 1], RGBColor[0, 0, 1], 
         RGBColor[0, 0, 1], RGBColor[0, 0, 1], RGBColor[0, 0, 1], 
         RGBColor[0, 0, 1], RGBColor[0, 0, 1], RGBColor[0, 0, 1], 
         RGBColor[0, 0, 1], RGBColor[0, 0, 1], RGBColor[0, 0, 1]}[[i]], Thin, 
       InfiniteLine[{a$*{2.4048339535892507, 3.831797113388671, 
            3.8317978287373955, 5.135990024261174, 5.136059096151318, 
            5.520754632240807, 6.3812950823735965, 6.381317681111972, 
            7.01775126007653, 7.0178224937838145, 7.590860601280539, 
            7.591072872794096, 8.421851138913722, 8.422710422186487, 
            8.660008044566908, 8.776751118278279, 8.776763030510674, 
            9.770470388828397, 9.771078698657918, 9.944801841701308, 
            9.946833640687165, 10.185705995432771, 10.185963561225158, 
            11.080011647362436, 11.084049565060361, 11.1030685286253, 
            11.103075194311309, 11.640891881099346, 11.646121888387668, 
            11.817807808363245, 12.251309264077156, 12.25325905283967, 
            12.366718068351073, 12.366959116079368, 13.053236153386042, 
            13.055560045790816, 13.367143379262203, 13.368658190586606, 
            13.396623457287058, 13.396635517477424, 13.631379740710331, 
            13.633949350140442, 14.427320132391074, 14.437548955448134, 
            14.533727179772093, 14.543976243054741, 14.863466256293904, 
            14.864300393737688, 14.885525269152437, 14.88564071984852}[[i]], 
         0}, {0, 1}]}, {i, Length[{2.4048339535892507, 3.831797113388671, 
         3.8317978287373955, 5.135990024261174, 5.136059096151318, 
         5.520754632240807, 6.3812950823735965, 6.381317681111972, 
         7.01775126007653, 7.0178224937838145, 7.590860601280539, 
         7.591072872794096, 8.421851138913722, 8.422710422186487, 
         8.660008044566908, 8.776751118278279, 8.776763030510674, 
         9.770470388828397, 9.771078698657918, 9.944801841701308, 
         9.946833640687165, 10.185705995432771, 10.185963561225158, 
         11.080011647362436, 11.084049565060361, 11.1030685286253, 
         11.103075194311309, 11.640891881099346, 11.646121888387668, 
         11.817807808363245, 12.251309264077156, 12.25325905283967, 
         12.366718068351073, 12.366959116079368, 13.053236153386042, 
         13.055560045790816, 13.367143379262203, 13.368658190586606, 
         13.396623457287058, 13.396635517477424, 13.631379740710331, 
         13.633949350140442, 14.427320132391074, 14.437548955448134, 
         14.533727179772093, 14.543976243054741, 14.863466256293904, 
         14.864300393737688, 14.885525269152437, 14.88564071984852}]}]]
 
Attributes[a$] = {Temporary}
 
necessaryfunctionslist = {pluck, peak, parametricregion, egsystem, 
     QOtriangleApprox, freqsintens, graphfi, keep, loudfreqs, graphweight, 
     colorfreqs, freqband, freqlines, print3Dborder, print3Dregion, p, vard, 
     disvector, soundanalysispg, comparetheovals}
 
necessaryfunctionslist /: necessaryfunctionslist::usage = "{pluck, peak, \
parametricregion, egsystem, QOtriangleApprox, freqsintens, graphfi, keep, \
loudfreqs, graphweight, colorfreqs, freqband, freqlines, print3Dborder, \
print3Dregion, p, vard, disvector,soundanalysispg,comparetheovals}"
