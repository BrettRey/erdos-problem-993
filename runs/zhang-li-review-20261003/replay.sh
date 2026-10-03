#!/bin/bash
# Local replay of github.com/selfreferencing/erdos993-lean at v1.0-claim (865e814).
set -u
cd "$(dirname "$0")/lean"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
log "HEAD $(git rev-parse HEAD)"; log "toolchain $(cat lean-toolchain)"
log "cache get"; lake exe cache get 2>&1 | tail -3
log "build default Erdos993Lean"; lake build Erdos993Lean 2>&1 | tail -5
for m in C02_20 C21_30 C31_39 C40 C41 C42 C43 C44 C45 C46 C47 C48 C49 C50 C51 C52 C53 C54 C55 C56 C57 C58 C59 C60; do
  log "kernel module $m"; lake build Erdos993Lean.ZhangKernel.Checks.$m 2>&1 | tail -2
done
log "build Erdos993LeanTheorem"; lake build Erdos993LeanTheorem 2>&1 | tee ../theorem_build.log | tail -8
log "sorry warnings in theorem build: $(grep -c "declaration uses 'sorry'" ../theorem_build.log)"
cat > ../AxiomCheck.lean <<'LEAN'
import Erdos993Lean.Analytic.Erdos993
import Erdos993Lean.ZhangKernel.Main
#print axioms Erdos993Lean.Analytic.erdos993
#print axioms Erdos993Lean.Zhang.forest_unimodal_of_card_le_sixty_kernel
#print axioms Erdos993Lean.Zhang.certificatesSound_kernel
LEAN
log "axioms"; lake env lean ../AxiomCheck.lean 2>&1
log "done"
