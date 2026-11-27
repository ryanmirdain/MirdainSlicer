"""Break down perf samples under one function.

    perf-any script -i perf.data --no-inline -F comm,tid,ip,sym | python3 perf_callers.py [MATCH]

MATCH (default 'unwarp_gcode') is a substring of the function of interest. Prints, for samples whose
stack contains it: its direct callees, its 3-deep caller chains, and the self symbols.
"""

import collections
import sys

match = sys.argv[1] if len(sys.argv) > 1 else 'unwarp_gcode'
self_c, child_c, caller_c, total = collections.Counter(), collections.Counter(), collections.Counter(), 0
stack = []


def flush():
    global total
    idx = next((i for i, n in enumerate(stack) if match in n), None)
    if idx is None:
        return
    total += 1
    self_c[stack[0][:110]] += 1
    child_c[(stack[idx - 1] if idx > 0 else 'SELF')[:110]] += 1
    caller_c[' <- '.join(n[:60] for n in stack[idx + 1:idx + 4])] += 1


for line in sys.stdin:
    line = line.rstrip('\n')
    if not line.strip():
        flush()
        stack.clear()
        continue
    if line.startswith(('\t', ' ')):
        parts = line.strip().split(' ', 1)
        if len(parts) == 2:
            stack.append(parts[1].rsplit(' (', 1)[0])
flush()
print('samples under', match, total)
for title, counter, count in (('direct callees', child_c, 15), ('callers (3 deep)', caller_c, 10),
                              ('self symbols', self_c, 25)):
    print('--', title)
    for name, c in counter.most_common(count):
        print(f'{c / max(total, 1):6.1%} {name}')
