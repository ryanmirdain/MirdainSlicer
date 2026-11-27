"""Leaf-symbol histogram per thread for a time window of a perf run.

    perf-any script -i perf.data -G -F comm,tid,time,ip,sym --no-inline | python3 perf_window.py [START_S] [END_S]

Times are seconds from the first sample. Prints the busiest threads in the window and the top leaf
symbols of the busiest one (main-thread phases whose call stacks perf cannot unwind, e.g. 3MF/XML).
"""
import collections
import re
import sys

start, end = (float(sys.argv[1]) if len(sys.argv) > 1 else 0.), (float(sys.argv[2]) if len(sys.argv) > 2 else 1e9)
first, threads, symbols = None, collections.Counter(), collections.defaultdict(collections.Counter)
for line in sys.stdin:
    m = re.match(r'\s*(\S+)\s+(\d+)\s+([\d.]+):\s+(.*)', line)
    if not m:
        continue
    t = float(m[3])
    first = t if first is None else first
    if not start <= t - first <= end:
        continue
    key = (m[1], int(m[2]))
    threads[key] += 1
    symbols[key][m[4].strip().split(' ', 1)[-1][:120]] += 1
for key, n in threads.most_common(5):
    print(key, n)
main = threads.most_common(1)[0][0]
print('-- top symbols of', main)
for name, n in symbols[main].most_common(30):
    print(f'{n / threads[main]:6.1%} {name}')
