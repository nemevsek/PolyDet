#!/usr/bin/env wolframscript
(* Build.wls -- build the documentation and pack the PolyDet paclet into a
   distributable .paclet archive (FindBounce-style GitHub release).
   Run:  wolframscript -file mathematica/Build.wls
   Output: mathematica/PolyDet/build/PolyDet-<version>.paclet  (with built docs).
   Source doc notebooks are authored by notebooks/build_polydet_docs.wls. *)
Needs["PacletTools`"];
dir = DirectoryName[$InputFileName];
src = FileNameJoin[{dir, "PolyDet"}];
out = FileNameJoin[{src, "build"}];
If[! DirectoryQ[out], CreateDirectory[out]];

(* PacletBuild builds the "Documentation" extension (-> indexed F1/search pages) AND
   assembles the archive with the BUILT docs.  It returns a Success[...] object and
   emits a few benign internal messages, so DO NOT wrap it in Check (that would force
   $Failed); just Quiet them and read PacletArchive off the result. *)
res = Quiet@PacletBuild[src, out];
archive = Which[
   MatchQ[res, _Success], res["PacletArchive"],
   AssociationQ[res] && KeyExistsQ[res, "PacletArchive"], res["PacletArchive"],
   StringQ[res] && FileExistsQ[res], res,
   True, $Failed];
If[! (StringQ[archive] && FileExistsQ[archive]),
   Print["[FAIL] PacletBuild returned ", Short[res]]; Exit[1]];
Print["[OK] built ", archive];

(* sanity: install, load the kernel, confirm the docs extension shipped *)
pac = PacletInstall[archive, ForceVersionInstall -> True];
Needs["PolyDet`"];
nbs = FileNames["*.nb", FileNameJoin[{pac["Location"], "Documentation"}], Infinity];
Print["installed ", pac["Name"], " ", pac["Version"],
   " ; PolyDetLnR ok: ", NumericQ[PolyDetLnR[3, {1.0, 2.0, 1.5}, {2.0, 3.0}]],
   " ; documentation notebooks: ", Length@nbs];
PacletUninstall[pac];   (* don't leave the dev build installed system-wide *)
Exit[0];
