Pseudobulk differential expression, one folder per contrast.

B4_vs_Interrupt/      paired within participant, participant as a blocking
                      factor. This is the contrast the design supports.
B4_vs_Control/        unpaired.
Interrupt_vs_Control/ unpaired.

The control arm is a single donor. BC_control_A and BC_control_B are two GEM
wells from one blood draw, so they are technical replicates. Both control
contrasts are 3 donors against 1, with dispersion estimated as if there were
2 control replicates; their p-values are anticonservative.

Donor, GEM well and dissociation stress are confounded with group. The
stress and immediate-early programme is about five times higher in the
control libraries than in the ATI libraries.

Columns added for this reason:
  tech_floor_logFC    absolute logFC between BC_control_A and BC_control_B
                      for the same gene and cell type. Same blood, so this is
                      the change the protocol produces with nothing differing.
  exceeds_tech_floor  abs(logFC) greater than tech_floor_logFC.
  stress_programme    gene is in the stress and immediate-early set.
  n_control_donors    1 for the control contrasts, recorded per row.
