import json
with open('results/analysis_n26.json') as f:
    data = json.load(f)
for failure in data['lc_failures']:
    print('LC failure poly:', failure['poly'][:15])
    print('  n:', failure['n'], 'lc_ratio:', failure['lc_ratio'])
