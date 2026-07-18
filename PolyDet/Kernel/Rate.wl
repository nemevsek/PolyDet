(* ::Package:: *)

(* PolyDet` -- unified multi-sector decay rate.  Loaded by Kernel/PolyDet.wl LAST
   (after Scalar/Multiscalar/Fermion/Gauge), since it assembles all of them.

   PolyDetDecayRate[Vpp, bf, extra] returns ln(Gamma/Volume): the driving (Higgs)
   scalar (from Vpp + bf, translational zero modes removed) plus a list of extra
   field sectors, combined with the correct one-loop signs:
     ln det_tot = Sigma^Higgs + Sum_scalars Sigma + Sum_gauge Sigma_gauge
                  - Sum_fermions Sigma^psi     (the (-1)^F),
   then  ln(Gamma/V) = (D/2) ln(S0/2pi) - (1/2) ln det_tot - S0 + zero-mode measure.
   All sectors share one MS-bar scale mu and the Simpson-averaged shell masses. *)

(* generic per-segment Simpson-averaged shell masses for ANY mass function fn(phi),
   structurally identical to PolyDetInputs (plateau = fn at the end nodes,
   wall = per-segment Simpson average).  Returns {shellMasses, radii}. *)
sectorMasses[fn_, bf_] := Module[{nodes, radii, av},
  nodes = Flatten[bf["Path"]];  radii = bf["Radii"];
  av = Table[(fn[nodes[[i]]] + 4 fn[(nodes[[i]] + nodes[[i + 1]])/2] + fn[nodes[[i + 1]]])/6,
     {i, Length[nodes] - 1}];
  If[Abs[N@radii[[1]]] < 10^-10,
    {Join[av, {fn[nodes[[-1]]]}], Rest[radii]},
    {Join[{fn[nodes[[1]]]}, av, {fn[nodes[[-1]]]}], radii}]];

(* shell FIELD values phi_s + radii (the gauge gradient jump is dmA = g Delta phi_s) *)
shellFields[bf_] := Module[{nodes, radii, seg},
  nodes = Flatten[bf["Path"]];  radii = bf["Radii"];  seg = (Most[nodes] + Rest[nodes])/2;
  If[Abs[N@radii[[1]]] < 10^-10,
    {Join[seg, {nodes[[-1]]}], Rest[radii]},
    {Join[{nodes[[1]]}, seg, {nodes[[-1]]}], radii}]];

(* one extra sector's signed contribution to ln det_tot.  tgt = the "Target"
   accuracy knob, threaded to the scalar/fermion sectors (gauge uses "Nsum"). *)
sectorContribution[spec_, bf_, phiS_, rE_, dim_, mu_, ig_, tgt_] := Module[
  {type = spec["Type"], mult = Lookup[spec, "Mult", 1], m2, r, mm, g, vaa, mA2, ma2, dmA},
  Switch[type,
    "Scalar",  {m2, r} = sectorMasses[spec["Mass2"], bf];
               mult PolyDetRen[m2, r, "Dimension" -> dim, "Scale" -> mu,
                  "ZeroRemoval" -> False, "Target" -> tgt],
    "Fermion", {mm, r} = sectorMasses[spec["Mass"], bf];          (* signed real mass M(phi) *)
               -mult PolyDetRenPsi[mm, r, "Scale" -> mu, "Target" -> tgt,
                  "Dimension" -> dim],                            (* the (-1)^F fermion sign *)
    "Gauge",   g = spec["Coupling"];  vaa = Lookup[spec, "Vaa", (0 &)];
               mA2 = g^2 phiS^2;  ma2 = g^2 phiS^2 + (vaa /@ phiS);  dmA = g Differences[phiS];
               mult PolyDetSigmaGauge[mA2, ma2, dmA, rE, "Ig" -> ig,
                  "Dimension" -> dim]["Sigma"],
    _, Message[PolyDetDecayRate::sector, type]; 0]];

PolyDetDecayRate::sector = "Unknown sector type `1`; use \"Scalar\", \"Fermion\", or \"Gauge\".";
Options[PolyDetDecayRate] = {"Scale" -> 1, "ZeroRemoval" -> True, "Goldstones" -> 0,
   "ChargeNorm" -> 1, "OrbitVolume" -> 2 Pi, "Ig" -> Automatic, "Target" -> Automatic,
   "Breakdown" -> False};
PolyDetDecayRate[Vpp_, bf_, extra_List : {}, OptionsPattern[]] := Module[
  {dim, S0, mu, tgt, m2H, rH, lndetH, phiS, rE, contribs, lndetTot, rate},
  dim = bf["Dimension"];  S0 = bf["Action"];  mu = OptionValue["Scale"];
  tgt = OptionValue["Target"];
  {m2H, rH} = PolyDetInputs[Vpp, bf];                                    (* driving Higgs *)
  lndetH = PolyDetRen[m2H, rH, "Dimension" -> dim, "Scale" -> mu,
     "ZeroRemoval" -> OptionValue["ZeroRemoval"], "Target" -> tgt];
  {phiS, rE} = shellFields[bf];
  contribs = sectorContribution[#, bf, phiS, rE, dim, mu, OptionValue["Ig"], tgt] & /@ extra;
  lndetTot = lndetH + Total[contribs];
  rate = PolyDetRate[lndetTot, S0, "Dimension" -> dim, "Goldstones" -> OptionValue["Goldstones"],
     "ChargeNorm" -> OptionValue["ChargeNorm"], "OrbitVolume" -> OptionValue["OrbitVolume"]];
  If[TrueQ[OptionValue["Breakdown"]],
     <|"Rate" -> rate, "lndetTotal" -> lndetTot, "Higgs" -> lndetH, "Sectors" -> contribs,
       "S0" -> S0, "Dimension" -> dim|>, rate]];
