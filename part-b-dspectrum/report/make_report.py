"""Builds report.md for Part B. All result tables are read from ../run_output.txt (no manual copying)."""
import re, os

HERE = os.path.dirname(os.path.abspath(__file__))
txt = open(os.path.join(HERE, '..', 'results', 'run_output.txt'), encoding='utf-8').read()


def block(name):
    return txt.split(name)[1].split('\n\n')[0].splitlines()


# Table 1
t1 = []
for line in block('Table 1'):
    m = re.match(r'\s*(\d+) \|\s*([\d.]+) \|\s*([\d.]+) \|\s*([\d.]+)', line)
    if m:
        t1.append(m.groups())
# Tables 2, 3
def ptable(name):
    rows = []
    for line in block(name):
        m = re.match(r'\s*([\d.]+) \|\s*([\d.]+) \|\s*([\d.]+) \|\s*([\d.]+)', line)
        if m:
            rows.append(m.groups())
    return rows
t2, t3 = ptable('Table 2'), ptable('Table 3')
# Table 4
t4, stats = [], {}
for line in block('Table 4'):
    m = re.match(r'\s*(\d+) \|\s*([\d.]+) \|\s*([\d.]+)', line)
    if m:
        t4.append(m.groups())
    m = re.match(r'\s*(mean|s|r\.e\.)\|\s*([\d.e+-]+) \|\s*([\d.e+-]+)', line)
    if m:
        stats[m.group(1)] = (float(m.group(2)), float(m.group(3)))
assert len(t1) == 30 and len(t2) == 11 and len(t3) == 11 and len(t4) == 10 and len(stats) == 3


def md_table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + ':---:|' * len(header)]
    out += ['| ' + ' | '.join(r) + ' |' for r in rows]
    return '\n'.join(out)


T1 = md_table(['i', 'M=1000 (f_i)', 'M=10000 (f_i)', 'M=50000 (f_i)'], t1)
T2 = md_table(['p', 'M=1000 (R)', 'M=10000 (R)', 'M=50000 (R)'], t2)
T3 = md_table(['p', 'M=1000 (R)', 'M=10000 (R)', 'M=50000 (R)'], t3)
T4 = md_table(['i', 'Destruction Spectrum', 'CMC'],
              [(i, f'{float(a):.6f}', f'{float(b):.4f}') for i, a, b in t4]
              + [('ממוצע', f"{stats['mean'][0]:.6f}", f"{stats['mean'][1]:.4f}"),
                 ('s', f"{stats['s'][0]:.2e}", f"{stats['s'][1]:.2e}"),
                 ('**r.e.**', f"**{stats['r.e.'][0]:.2e}**", f"**{stats['r.e.'][1]:.2e}**")])

ds_mean, cmc_mean = stats['mean']
ds_s, cmc_s = stats['s']
ds_re, cmc_re = stats['r.e.']
ratio = cmc_re / ds_re
n_ones = sum(1 for _, _, b in t4 if float(b) == 1.0)

tpl = open(os.path.join(HERE, 'report_template.md'), encoding='utf-8').read()
out = (tpl.replace('{{T1}}', T1).replace('{{T2}}', T2).replace('{{T3}}', T3).replace('{{T4}}', T4)
          .replace('{{DS_MEAN}}', f'{ds_mean:.6f}').replace('{{CMC_MEAN}}', f'{cmc_mean:.4f}')
          .replace('{{DS_S}}', f'{ds_s:.2e}').replace('{{CMC_S}}', f'{cmc_s:.2e}')
          .replace('{{DS_RE}}', f'{ds_re:.2e}').replace('{{CMC_RE}}', f'{cmc_re:.2e}')
          .replace('{{DS_REP}}', f'{100 * ds_re:.3f}%').replace('{{CMC_REP}}', f'{100 * cmc_re:.3f}%')
          .replace('{{RATIO}}', f'{ratio:.1f}').replace('{{N_ONES}}', str(n_ones)))
assert '{{' not in out
open(os.path.join(HERE, 'report.md'), 'w', encoding='utf-8').write(out)
print('report.md written;', 'ratio', round(ratio, 2), 'ones', n_ones)
