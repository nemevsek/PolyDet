#!/usr/bin/env wolframscript
(* Run the PolyDet` test suite and report PASS/FAIL.
   Usage: wolframscript -file mathematica/PolyDet/Tests/run.wls
   (exit code is the MKL-shutdown artifact 3 on this machine -- rely on the printed line). *)
rep = TestReport[FileNameJoin[{DirectoryName[$InputFileName], "PolyDet.wlt"}]];
ok = rep["AllTestsSucceeded"];
np = rep["TestsSucceededCount"]; nf = rep["TestsFailedCount"];
Print[If[ok, "[ALL PASS] " <> ToString[np] <> " PolyDet tests",
   "[FAIL] " <> ToString[nf] <> "/" <> ToString[np + nf] <> " -- " <>
   ToString[Keys[rep["TestsFailed"]]]]];
Exit[If[ok, 0, 1]];
