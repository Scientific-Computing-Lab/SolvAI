# Cover letter

Dear Editors,

We submit “Distilling solvation physics into AI for simulation-free hydration
free energies” for consideration as an Article in *Nature*.

Hydration free energy is a basic molecular quantity, but high accuracy can
require explicit solvent and path-integral sampling for every new solute. We
show a different route: physical calculations are used once to define reusable
solvent-response supervision, which is distilled into a system that accepts
only molecular structure at inference. On the same 85-solute reference set used
to establish the ARROW/PIMD8 result, SolvAI reaches a strict five-fold OOF MAE
of 0.197 kcal/mol and a five-partition mean of 0.204 ± 0.005 kcal/mol, while
running no MD, PIMD or probe calculation at deployment.

The conceptual advance extends beyond this particular thermodynamic endpoint.
It suggests that expensive simulation can train transferable physical response
coordinates and thereby be amortized across future molecules. The complete
model, held-out predictions, leakage audit, provenance record, figure source
and manuscript build are supplied in a reproducible public repository.

The manuscript does not claim generic FreeSolv state of the art or robust
sub-0.20 accuracy. It makes a direct, leakage-controlled comparison with the
ARROW/PIMD8 reference chemistry and reports the harder family and scaffold
tests.

Sincerely,

Gal Oren and Michael Levitt
