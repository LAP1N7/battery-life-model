# D2 summary

rule: top of charge, thickness median at rest time 400-580 s; reference = first check-up (BOL).

- rest-rule sensitivity: 199 check-ups have a ~1 h rest; |X(400-580 s) - X(end of rest)|: median 0.10 um, 95th pct 0.60, max 1.20 um (irreversible changes are 20-500 um, so the rest rule is not a limiting factor).
- short rests (< 580 s): 0 check-ups (they are excluded from X_top and show NaN).
- artifact cells flagged automatically: 03 (step_event), 09 (baseline_drop), 13 (sensor_nan_or_frozen), 16 (sensor_nan_or_frozen), 21 (baseline_drop,step_event). Previously known: 03, 09, 13, 21.
  - step event: cell 03 at check-up 8 (684 Ah)
  - step event: cell 21 at check-up 9 (1746 Ah)
- vs team out/expansion.csv (x_irr = bottom-state thickness, first 20 charge samples, no rest rule): matched 223 check-ups after BOL, Spearman(top-state dX_irr, x_irr) 0.917; difference median -28.5 um, 95th pct |diff| 98.0 um. My X_bottom (first 20 samples of the C/20 charge) differs from team x_dis by more than 1 um in 40 check-ups (cells 03, 06, 09, 12, 15, 18, 21; median offset 99 um): the team's selection 100<I<400 mA picks up the CV-taper rows of a preceding 1C charge at the top of charge, so for these check-ups x_dis is a TOP-state thickness while for the others it is a bottom-state thickness (mixed states in x_irr).

- last VALID check-up dX_irr_um by temperature group (median [min, max]):
  - C: 192 [118, 533] um (n=7)
  - H: 163 [34, 348] um (n=7)
  - R: 184 [131, 284] um (n=7)