# Pair-AB Push Fast2 Scan (k9/g6/s30)

## Source artifact
- Downloaded JSON: `~/Downloads/AB_attachment_search_result.json`
- Imported as: `results/fast2_k9_g6_s30.json`

## Reported command (from JSON)
- `python3 /mnt/data/pair_AB_attachment_search_push_fast2.py --gadget-max-n 9 --max-gadgets 6 --max-added-size 30 --m-target 5 --lambda-target 7/9 --out /mnt/data/fast2_k9_g6_s30.json`

## Outcome
- `witness_found=false`
- shared key set: `[(deltaN=1, m=5, lambda=7/9)]`
- frontier top miss:
  - `rhoA=22901648/33160527`
  - `rhoB=8295005488/16144619103`
  - `|rhoA-rhoB|=10519070037889952/59484897518860809`

## Bounded accounting (from certificate)
- gadget universe size: `124`
- patterns enumerated: `3759303`
- `good_patterns_A_count=1`, `good_patterns_B_count=1`
- bottleneck: enumeration

## Integrity check
- Recomputed `|rhoA-rhoB|` from reported fractions and confirmed exact match to JSON `rho_gap`.

## Missing companion script
- The generator script path in JSON is `/mnt/data/pair_AB_attachment_search_push_fast2.py`.
- That `.py` file is not present in `~/Downloads` yet, so only the output artifact was recoverable this round.
