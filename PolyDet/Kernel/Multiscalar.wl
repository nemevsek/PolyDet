(* ::Package:: *)

(* PolyDet` -- coupled multi-scalar + rate section.  Loaded by Kernel/PolyDet.wl from inside
   PolyDet`Private` (Get'd by the loader); NOT a standalone package -- it has no
   BeginPackage/Begin of its own and relies on the loader's context.
   Matrix Bessel blocks (Cayley-Hamilton), PolyDet*Mat, matrix zero removal
   PolyDetR1PrimeMat, and the nucleation-rate assembly PolyDetRate.
   Depends on the scalar section (Ip/GP, nuOf/dDeg/ZetaD, numinAuto, PolyDetLnR). *)

(* matrix Bessel building blocks via Cayley-Hamilton (MatrixFunction;
   eigenvalues only, no eigenvectors): I_nu(M R), K_nu(M R), M I_{nu-1}(M R),
   M K_{nu-1}(M R) with M = Sqrt[m2]. *)
(* MatrixFunction::valtlrg is a conservative "value >> argument" warning that fires
   on the growing-mode (large-nu / large mR) Bessels; the per-mode WP scaling in
   PolyDetSigmaMat already compensates, so it is handled noise -- Quiet it
   (otherwise it spams the user on every matrix evaluation). *)
matBlk[nu_, R_, m2_] := Quiet[<|
   "I"   -> MatrixFunction[BesselI[nu, Sqrt[#] R] &, m2],
   "K"   -> MatrixFunction[BesselK[nu, Sqrt[#] R] &, m2],
   "MI1" -> MatrixFunction[Sqrt[#] BesselI[nu - 1, Sqrt[#] R] &, m2],
   "MK1" -> MatrixFunction[Sqrt[#] BesselK[nu - 1, Sqrt[#] R] &, m2]|>, MatrixFunction::valtlrg];

(* matrix transfer T^s = G_{s+1}^{-1} G_s.  Each block term is kept as
   (fn of m2_{s+1}) . (fn of m2_s): M' and M do NOT commute, so the scalar
   factoring (ms1 I0 Kf) would reorder them and is WRONG for matrices.
   det_{2N_f}(T^s) = 1 preserved (per-segment Wronskian is scalar). *)
TsMat[nu_, R_, a2_, b2_] := Module[{A = matBlk[nu, R, a2], B = matBlk[nu, R, b2]},
  R ArrayFlatten[{
     {B["MK1"] . A["I"] + B["K"] . A["MI1"],   B["MK1"] . A["K"] - B["K"] . A["MK1"]},
     {B["MI1"] . A["I"] - B["I"] . A["MI1"],   B["MI1"] . A["K"] + B["I"] . A["MK1"]}}]];

(* ---- matrix translational zero-mode removal (R_1') -------------------------
   Matrix promotion of the scalar dTs source: d/d(mu^2) of the matrix Bessel blocks
   under the ISOTROPIC offset m^2 -> m^2 + mu^2*1.  Differentiate each scalar Bessel
   per eigenvalue (MatrixFunction does the symbolic d/dx; eigenvalues only). *)
matBlkD[nu_, R_, m2_] := Module[{x},
  Quiet[<|"I"   -> MatrixFunction[Function[xx, Evaluate[D[BesselI[nu, Sqrt[x] R], x] /. x -> xx]], m2],
    "K"   -> MatrixFunction[Function[xx, Evaluate[D[BesselK[nu, Sqrt[x] R], x] /. x -> xx]], m2],
    "MI1" -> MatrixFunction[Function[xx, Evaluate[D[Sqrt[x] BesselI[nu - 1, Sqrt[x] R], x] /. x -> xx]], m2],
    "MK1" -> MatrixFunction[Function[xx, Evaluate[D[Sqrt[x] BesselK[nu - 1, Sqrt[x] R], x] /. x -> xx]], m2]|>,
    MatrixFunction::valtlrg]];

(* d T^s / d(mu^2): slot 1 differentiates the A (inner m^2_s) blocks of TsMat, slot 2
   the B (outer m^2_{s+1}) blocks -- same block structure/ordering as TsMat.
   det T^s = 1 is preserved (d/d(mu^2) det T^s = 0). *)
dTsMat[nu_, R_, a2_, b2_, slot_] := Module[{A = matBlk[nu, R, a2], B = matBlk[nu, R, b2], dA, dB},
  If[slot === 1, dA = matBlkD[nu, R, a2];
   R ArrayFlatten[{{B["MK1"] . dA["I"] + B["K"] . dA["MI1"], B["MK1"] . dA["K"] - B["K"] . dA["MK1"]},
                   {B["MI1"] . dA["I"] - B["I"] . dA["MI1"], B["MI1"] . dA["K"] + B["I"] . dA["MK1"]}}],
   dB = matBlkD[nu, R, b2];
   R ArrayFlatten[{{dB["MK1"] . A["I"] + dB["K"] . A["MI1"], dB["MK1"] . A["K"] - dB["K"] . A["MK1"]},
                   {dB["MI1"] . A["I"] - dB["I"] . A["MI1"], dB["MI1"] . A["K"] + dB["I"] . A["MK1"]}}]]];

(* general adjugate via Faddeev-LeVerrier: division-free (only by the integer k),
   stable for the rank-(N_f-1) near-singular T11; reduces to Tr[A]1 - A for N_f=2.
   Works for any N_f. *)
adjMat[A_] := Module[{n = Length[A], id = IdentityMatrix[Length[A]], M},
  M = id;  Do[M = A . M - (Tr[A . M]/k) id, {k, 1, n - 1}];  (-1)^(n - 1) M];
(* right null vector of the near-singular T (the zero-mode direction): the columns of
   adj(T) span ker T when rank = N_f-1, so take the largest-norm column.  Exact, stable,
   any N_f; only the direction is used (ref is a Rayleigh quotient, sign-free). *)
nullVec[adjT_] := Module[{cols = Transpose[adjT]}, cols[[First@Ordering[Norm /@ cols, -1]]]];

(* matrix R_1' : reduced l=1 (nu=D/2) determinant.  Offset interior shells
   m^2_s -> m^2_s + mu^2*1 (FV fixed); det R_{l=1}(0)=0 (ONE zero eigenvalue along the
   bounce tangent, for any N_f), so R_1' is the linear coefficient -- d/d(mu^2) det R_{l=1}
   via Jacobi/adjugate, scaled by ref = FV mass projected on ker(T11) (option a). *)
r1pMat[m2_, radii_, dim_, prec_, l_:1] := Module[
  {nu = l + dim/2 - 1, nf = Length[m2[[1]]], mm, rr, M, Tlist, Tt, dTt, T11, dT11, m1, mN, adjT, detderiv, vk, ref},
  mm = SetPrecision[m2, prec];  rr = SetPrecision[radii, prec];  M = Length[radii];
  m1 = mm[[1]];  mN = mm[[-1]];
  Tlist = Table[TsMat[nu, rr[[s]], mm[[s]], mm[[s + 1]]], {s, M}];  Tt = Dot @@ Reverse[Tlist];
  dTt = Sum[Dot @@ Reverse[ReplacePart[Tlist, s ->
        dTsMat[nu, rr[[s]], mm[[s]], mm[[s + 1]], 1] +
        If[s < M, dTsMat[nu, rr[[s]], mm[[s]], mm[[s + 1]], 2], 0]]], {s, M}];
  T11 = Tt[[1 ;; nf, 1 ;; nf]];  dT11 = dTt[[1 ;; nf, 1 ;; nf]];
  adjT = adjMat[T11];                                   (* general adjugate, reused for ker *)
  detderiv = (Det[mN]/Det[m1])^(nu/2) (Tr[adjT . dT11] - (nu/2) Tr[Inverse[m1]] Det[T11]);
  vk = nullVec[adjT];  ref = Re[(vk . mN . vk)/(vk . vk)];
  <|"R1p" -> Re[ref detderiv], "ref" -> ref|>];

Options[PolyDetR1PrimeMat] = {"Dimension" -> 4, "Multipole" -> 1, WorkingPrecision -> 80};
PolyDetR1PrimeMat[m2_, radii_, OptionsPattern[]] :=
  r1pMat[m2, radii, OptionValue["Dimension"], OptionValue[WorkingPrecision],
    OptionValue["Multipole"]]["R1p"];

PolyDetLnRMat[nu_, m2_, radii_, prec_: 60] := Module[{m2p, rr, nf, Tlist, Tt, T11, m1, mN},
  m2p = SetPrecision[m2, prec];  rr = SetPrecision[radii, prec];  nf = Length[m2p[[1]]];
  Tlist = Table[TsMat[nu, rr[[s]], m2p[[s]], m2p[[s + 1]]], {s, Length[rr]}];
  Tt = Dot @@ Reverse[Tlist];  T11 = Tt[[1 ;; nf, 1 ;; nf]];
  m1 = m2p[[1]];  mN = m2p[[-1]];
  Re[(nu/2) (Log[Det[mN]] - Log[Det[m1]]) + Log[Det[T11]]]];

PolyDetXMaxMat[m2_, radii_] := Max[Table[
   radii[[s]] Sqrt[Max[Abs[Eigenvalues[N[m2[[s]]]]], Abs[Eigenvalues[N[m2[[s + 1]]]]]]],
   {s, Length[radii]}]];

(* matrix I_p: the scalar table GP under a flavour trace with GROUPED matrix
   powers, tr(UX^i . UY^j) (cyclic-safe; UX = m_s^2 R^2/4, UY = m_{s+1}^2 R^2/4).
   Divergent orders p=1,3 are single-matrix powers (tr dm^2, tr dm^4 --
   basis-independent); higher orders feed the optimal-truncation accel + zeta
   add-back. *)
IpMat[p_, UX_, UY_] := Total[
   (#[[2]] Tr[MatrixPower[UX, #[[1, 1]]] . MatrixPower[UY, #[[1, 2]]]]) & /@ GP[p]];
PolyDetIpMat[m2_, radii_, pmax_] := Table[
   Sum[IpMat[p, m2[[s]] radii[[s]]^2/4, m2[[s + 1]] radii[[s]]^2/4], {s, Length[radii]}],
   {p, 1, pmax}];

(* ---- FrozenWall matrix band ----
   Case-b thin wall, N_f coupled fields.  Key collapse: the frozen box operator
   K^2 = m^2_s + (nu^2 - 1/4) c^2 * Id SHARES the shell mass matrix's eigenvectors, so
   each shell propagates DIAGONALLY in its own (rho-independent, nu-independent,
   precomputed) eigenbasis -- N_f scalar frozen boxes per shell -- and the rotating
   eigenbasis enters ONLY as the frame rotation V_{s+1}.V_s^T at each matching (psi and
   psi' are continuous; no jumps).  Exact matrix content: N_f plateau Bessel
   log-derivatives (eigenbasis of m^2_1) and N_f FV Wronskian rows (eigenbasis of
   m^2_N), 4 N_f scalar Bessels per mode, N-independent, machine band.  Requires
   positive-definite plateau and FV mass matrices (case-b guard); tachyonic wall
   eigenvalues just make kappa imaginary (exact per box). *)
lnRMatFrozen[a___] := Quiet[lnRMatFrozenC[a], {General::munfl}];
lnRMatFrozenC[nu_, m2_, radii_, nsub_, precB_: 40] := Module[
  {nf, rr, nsh, lamV, lam1, V1, lamN, VN, R0, RN, x0, xN, sig0, lnI0, sigK, lnKN,
   U, Up, Q, kap, ch, sh, a, b, h, c, mag, lnscale = 0.},
  nf = Length[m2[[1]]];  rr = N[radii];  nsh = Length[radii];
  lamV = Table[Eigensystem[N[m2[[s]]]], {s, Length[m2]}];   (* {vals, vecs(rows)} per shell *)
  {lam1, V1} = lamV[[1]];  {lamN, VN} = lamV[[-1]];
  R0 = SetPrecision[radii[[1]], precB];  RN = SetPrecision[radii[[-1]], precB];
  x0 = Sqrt[SetPrecision[lam1, precB]] R0;  xN = Sqrt[SetPrecision[lamN, precB]] RN;
  sig0 = N[1/(2 R0) + Sqrt[SetPrecision[lam1, precB]] BesselI[nu - 1, x0]/BesselI[nu, x0] - nu/R0];
  lnI0 = N[Log[Sqrt[R0]] + Log[BesselI[nu, x0]]];
  sigK = N[1/(2 RN) - Sqrt[SetPrecision[lamN, precB]] BesselK[nu - 1, xN]/BesselK[nu, xN] - nu/RN];
  lnKN = N[Log[Sqrt[RN]] + Log[BesselK[nu, xN]]];
  U = N[IdentityMatrix[nf]];  Up = DiagonalMatrix[sig0];    (* frame: eigenbasis of shell 1 *)
  Do[
    Q = lamV[[s + 1, 2]] . Transpose[lamV[[s, 2]]];          (* rotate frame s -> s+1 at rr[[s]] *)
    U = Q . U;  Up = Q . Up;
    Do[
      a = rr[[s]] + (rr[[s + 1]] - rr[[s]]) (q - 1)/nsub;
      b = rr[[s]] + (rr[[s + 1]] - rr[[s]]) q/nsub;  h = b - a;  c = Log[b/a]/h;
      kap = Sqrt[lamV[[s + 1, 1]] + (nu^2 - 1/4) c^2 + 0. I];
      ch = Cosh[kap h];  sh = Sinh[kap h]/kap;
      {U, Up} = {ch U + sh Up, (kap^2 sh) U + ch Up},        (* diagonal: per-row scaling *)
      {q, nsub}];
    mag = Max[Abs[U], Abs[Up]];
    If[mag > 10.^15, U /= mag; Up /= mag; lnscale += Log[mag]],
    {s, nsh - 1}];
  Q = VN . Transpose[lamV[[-2, 2]]];  U = Q . U;  Up = Q . Up;   (* FV frame at rr[[-1]] *)
  Re[(nu/2) (Total[Log[SetPrecision[lamN, precB]]] - Total[Log[SetPrecision[lam1, precB]]]) +
     Total[lnI0] + Total[lnKN] + nf lnscale + Log[Det[Up - sigK U]]]];

Options[PolyDetSigmaMat] = Options[PolyDetSigma];
PolyDetSigmaMat::frozen = "\"Method\" -> \"FrozenWall\" needs a case-b-like chain (positive-definite plateau and FV mass matrices, ordered radii); falling back to \"Exact\".";
PolyDetSigmaMat[m2_, radii_, OptionsPattern[]] := Module[
  {dim, precOpt, nseg, precOf, pmax, tgt, xs, Ns, meth, lstar, caseb, Ip, sub, lnRv, direct, addback, nu},
  dim = OptionValue["Dimension"];  pmax = OptionValue["Pmax"];  tgt = OptionValue["Target"];
  xs = PolyDetXMaxMat[m2, radii];  nseg = Length[radii];
  Ns = Replace[OptionValue["Nsum"], Automatic :>
       If[NumericQ[tgt], Max[Ceiling[xs] + 2, Ceiling[xs tgt^(-1/(pmax + 2 - dim))]],
          Max[100, Ceiling[5 xs]]]];
  (* Matrices need MUCH more precision than the scalar chain.  At MachinePrecision
     high-nu Bessels underflow to an EXACT zero matrix (singular chain -> Det=0 ->
     Log[0]); at fixed arbitrary precision the matrix product loses digits set by the
     accumulated dynamic range of the chained MatrixFunctions (spanning e^{+-nu eta}).
     That range grows with the segment count AND with how large the Bessel arguments
     get (the turning point x_max): a stiff, tachyonic, large-x_max chain (e.g. a
     coupled Higgs-singlet bounce) loses ~2-3x more than nu*nseg/5.  So scale WP with
     nu, the segment count, AND x_max, PER MODE: prec(nu) ~ nu*(nseg + 3 x_max)/5 + 30
     (>=40).  An explicit WorkingPrecision overrides for all modes.
     ("Method" -> "FrozenWall" sidesteps ALL of this for the l >= 2 band: the frozen
     evaluator is machine-precision by construction -- the worst matrix-WP stress case
     removed; l = 0, 1 stay on the exact chain.) *)
  precOpt = OptionValue[WorkingPrecision];
  precOf[v_] := Replace[precOpt, Automatic :> Ceiling[Max[40, v (nseg + 3 xs)/5 + 30]]];
  caseb := Min[Eigenvalues[N[m2[[1]]]]] > 0 && Min[Eigenvalues[N[m2[[-1]]]]] > 0 &&
           N[radii[[1]]/radii[[-1]]] >= 1/2;
  meth = Replace[OptionValue["Method"], Automatic :>
     If[NumericQ[tgt] && tgt >= 1/1000 && xs >= 30 && Length[radii] >= 3 && caseb,
        "FrozenWall", "Exact"]];
  If[meth === "FrozenWall" && ! (Min[Eigenvalues[N[m2[[1]]]]] > 0 && Min[Eigenvalues[N[m2[[-1]]]]] > 0),
     Message[PolyDetSigmaMat::frozen]; meth = "Exact"];
  Ip = PolyDetIpMat[m2, radii, pmax];
  sub[v_] := Sum[Ip[[p]]/v^p, {p, 1, pmax}];
  (* offset (Hurwitz) assembly, as the scalar "SubtractFrom" (the matrix I_p are traces,
     so the split is identical; kills the x_max^(pmax-dim+1) cancellation) *)
  lstar = Replace[OptionValue["SubtractFrom"], {
     Automatic :> If[xs >= 30, lminOf[Ceiling[xs], dim], 0],
     n_?NumericQ /; n > 0 :> lminOf[n, dim], _ :> 0}];
  lstar = Min[lstar, Ns];
  lnRv = If[meth === "FrozenWall",
     SetPrecision[
       Join[Table[PolyDetLnRMat[nuOf[l, dim], m2, radii, precOf[nuOf[l, dim]]], {l, 0, Min[1, Ns - 1]}],
            Table[lnRMatFrozen[nuOf[l, dim], m2, radii, OptionValue["NSub"]], {l, 2, Ns - 1}]],
       If[lstar > 0, 40, Max[40, (pmax - dim + 1) Log10[Max[xs, 2]] + 25]]],
     Table[PolyDetLnRMat[nuOf[l, dim], m2, radii, precOf[nuOf[l, dim]]], {l, 0, Ns - 1}]];
  direct = Sum[nu = nuOf[l, dim];
     dDeg[l, dim] (lnRv[[l + 1]] - If[l >= lstar, sub[nu], 0]), {l, 0, Ns - 1}];
  addback = If[lstar == 0,
     Sum[Ip[[p]] ZetaD[p, dim], {p, dim, pmax}],
     Sum[Ip[[p]] ZetaDOff[p, dim, lstar], {p, dim, pmax}] -
       Sum[Ip[[p]] Sum[dDeg[l, dim]/nuOf[l, dim]^p, {l, 0, lstar - 1}], {p, 1, dim - 1}]];
  direct + addback];

i3LogMat[m2_, radii_] := -(1/32) Sum[
   radii[[s]]^4 Tr[MatrixPower[m2[[s]], 2] - MatrixPower[m2[[s + 1]], 2]] Log[radii[[s]]],
   {s, Length[radii]}];
PolyDetRenMat::gold = "\"Goldstones\"->`1`: only one (corank-1) Goldstone is supported \
(the rank-1 adjugate); for several equivalent Goldstones remove each as a scalar l=0 \
mode on the Goldstone mass chain.";
Options[PolyDetRenMat] = Join[{"Scale" -> 1, "ZeroRemoval" -> True, "Goldstones" -> 0}, Options[PolyDetSigma]];
PolyDetRenMat[m2_, radii_, opts : OptionsPattern[]] := Module[
  {dim, mu, prec, sig, i3, muT, ct, zr, zrG, zrPrec, info, infoG, nG},
  dim = OptionValue["Dimension"];  mu = OptionValue["Scale"];  nG = OptionValue["Goldstones"];
  prec = Replace[OptionValue[WorkingPrecision], Automatic -> MachinePrecision];
  zrPrec = Max[60, Ceiling[2 PolyDetXMaxMat[m2, radii]]];
  sig = PolyDetSigmaMat[m2, radii, FilterRules[{opts}, Options[PolyDetSigma]]];
  (* l=1 translational zero modes: PolyDetSigmaMat includes the (artifact) term
     d_1 ln|det R_{l=1}| (-> -inf as the polygon refines).  Replace it with the
     removed-zero d_1 ln|R_1'/ref|, R_1' the matrix reduced l=1 determinant
     (PolyDetR1PrimeMat), ref = FV mass on the zero-mode direction. *)
  zr = If[TrueQ[OptionValue["ZeroRemoval"]],
     info = r1pMat[m2, radii, dim, zrPrec];
     dDeg[1, dim] (Log[Abs[info["R1p"]/info["ref"]]] - PolyDetLnRMat[dim/2, m2, radii, zrPrec]),
     0];
  (* l=0 GOLDSTONE zero modes (broken global symmetry): same removal at nu=D/2-1,
     d_0 ln|R_0'/ref|, R_0' the l=0 reduced determinant (kernel = transverse Goldstone
     direction; the radial negative mode stays inside R_0').  nG = #broken generators. *)
  zrG = If[nG > 0,
     If[nG > 1, Message[PolyDetRenMat::gold, nG]];
     infoG = r1pMat[m2, radii, dim, zrPrec, 0];
     nG dDeg[0, dim] (Log[Abs[infoG["R1p"]/infoG["ref"]]] - PolyDetLnRMat[dim/2 - 1, m2, radii, zrPrec]),
     0];
  ct = Which[
    dim == 4, i3 = PolyDetIpMat[m2, radii, 3][[3]];
      muT = Sqrt[4 Pi] Exp[-EulerGamma/2] SetPrecision[mu, prec];
      i3 (3/4 + EulerGamma + Log[muT/2]) + i3LogMat[m2, radii],
    dim == 3, 0,
    True, Message[PolyDetRen::dim, dim]; Return[$Failed]];
  sig + ct + zr + zrG];

(* ---- collective-coordinate assembly: the bounce-nucleation rate prefactor -------
   Callan-Coleman:  Gamma/V = (S0/2pi)^{D/2} |det'/det_FV|^{-1/2} e^{-S0},  
   det' the reduced determinant with the D translational zero modes removed.  
   
   By the McKane-Tarlie construction, PolyDet's lndetRen
   (= ln|det'/det_FV|, ZeroRemoval->True) IS that det', so
       ln(Gamma/V) = (D/2) ln(S0/2pi) - (1/2) lndetRen - S0 .
       
   Each broken-symmetry Goldstone (removed in lndetRen, "Goldstones"->nG) adds the
   orbit measure  (Q/2pi)^{1/2} Vol(G/H)  per mode, where Q MUST be the FULL
   D-dimensional zero-mode norm
       Q = int d^Dx (T^a phibar)^2 = Omega_{D-1} * int rho^{D-1} (T^a phibar)^2 drho ,
   Omega_{D-1} = 2 pi^{D/2}/Gamma(D/2)  (2 pi^2 in D=4, 4 pi in D=3), with the
   generator normalisation matching the supplied "OrbitVolume" (e.g. T^a = sigma^a/2
   with Vol(S^3) = 2 pi^2).  R_0' inside lndetRen is a pure eigenvalue-product ratio
   and carries NO norm, so the measure norm is NOT a matter of chain convention --
   an earlier version of this comment claimed Q "cancels against R_0'" and could be
   passed in any consistent normalisation; that was WRONG (adjudicated 2026-07-08,
   python/verify_goldstone_measure.py: R_0' == reduced eigenvalue product to 0.2%).
   This is the same full-volume convention as the translational S0 =
   (Omega_{D-1}/D) int rho^{D-1} phidot^2 (Derrick), and as AFS 1707.08124 (J_G).
   
   The single negative mode (l=0 radial) is in |det'|; Coleman's negative-mode
   contour-1/2 and the 2 of Gamma=2 Im CANCEL exactly, so this is the COMPLETE rate
   (no extra factor) for the standard one-negative-mode bounce. *)
   
Options[PolyDetRate] = {
  "Dimension" -> 4, 
  "Goldstones" -> 0, 
  "ChargeNorm" -> 1, 
  "OrbitVolume" -> 2 Pi};

PolyDetRate[lndetRen_, S0_, OptionsPattern[]] := Module[{dim, nG, Q, vol},
  dim = OptionValue["Dimension"];  
  nG = OptionValue["Goldstones"];
  Q = OptionValue["ChargeNorm"];  
  vol = OptionValue["OrbitVolume"];
  
  (dim/2) Log[S0/(2 Pi)] - lndetRen/2 - S0 + If[nG > 0, nG (1/2 Log[Q/(2 Pi)]) + Log[vol], 0]];

