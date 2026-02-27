# Pair A/B Multiset-Subset Scan (reported by C, 2026-02-28)

## Reported script/command
- Script: `/mnt/data/pair_AB_attachment_search_multiset_subset.py`
- Command:
  - `python3 /mnt/data/pair_AB_attachment_search_multiset_subset.py --gadget-max-n 9 --max-gadgets 5 --max-added-size 30 --m-min 4 --subset-size 75 --out /mnt/data/pair_AB_attachment_search_multiset_sub75_k9_g5_s30.json`

## Reported artifact
- `/mnt/data/pair_AB_attachment_search_multiset_sub75_k9_g5_s30.json`

## Reported outcome
- `witness_found=false`
- Shared-key frontier contains only:
  - `(ΔN=0, m=5, λ=7/9)`
- Best/only reported rho gap at that shared key:
  - `rhoA=22901648/33160527`
  - `rhoB=8295005488/16144619103`
  - `|rhoA-rhoB|=10519070037889952/59484897518860809`

## Sync status
- This artifact is not yet present in repo `results/`.
- Local corroboration currently available in:
  - `results/frontier_adjacent_overlap_attachment_search_k8_g4_s24_local.json`
  - `notes/frontier_adjacent_overlap_attachment_search_2026-02-28.md`
