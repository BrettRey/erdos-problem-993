# Pair-AB Attachment Search (subset 75, k9/g5/s30)

## Script and command
- Script: `scripts/pair_ab_attachment_search_multiset_subset.py`
- Command:
  - `python3 scripts/pair_ab_attachment_search_multiset_subset.py --gadget-max-n 9 --max-gadgets 5 --max-added-size 30 --m-min 4 --subset-size 75 --seed 0 --out results/pair_AB_attachment_search_multiset_sub75_k9_g5_s30.json`

## Artifact
- `results/pair_AB_attachment_search_multiset_sub75_k9_g5_s30.json`

## Totals
- `checked_total=2087324`
- `checked_A=1043662`
- `checked_B=1043662`
- `kept_A=1043662`
- `kept_B=1043662`
- `shared_keys=1`
- `witness_found=false`

## Frontier (top 1)
- `m=5`, `lambda=7/9`
- `NA=13`, `NB=14`, `deltaN=1`
- `rhoA=22901648/33160527`
- `rhoB=8295005488/16144619103`
- `|rhoA-rhoB|=10519070037889952/59484897518860809`

## Reporter
- `python3 scripts/report_pair_ab_frontier.py --in results/pair_AB_attachment_search_multiset_sub75_k9_g5_s30.json --top 5`
