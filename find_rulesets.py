"""Find 4-state ring automata that synchronize for a specified set of
lengths."""

from itertools import product
from os import fsync
from time import clock
from ruleset_utils import *


def defsim_ring(outfile, init_ruleset, lengths, nonF_states,
                sym=False, opt=False):
  f = open(outfile, 'w')
  f.close()

  cnt = 1
  def save_ruleset(ruleset):
    f = open(outfile, 'a')
    f.write(str(ruleset) + ',\n')
    """Merely closing the file does not guarantee writing the
    data to disk immediately."""
    f.flush()
    fsync(f.fileno())
    f.close()
    nonlocal cnt
    print('Ruleset ' + str(cnt))
    cnt += 1
  
  def defsim_cfg(cfg, cfgs, ruleset, lidx):
    """Define and simulate <seg>."""
    seg = cfg + cfg[:2] # Unroll ring array.
    udef_segs, def_states = get_udef_segs_def_states(seg,
                                                     ruleset)
    # length-n time-optimal rings must fire at time-step n
    if opt and len(cfgs) == len(cfg):
      if def_states - {'F'}:
        return None
      ruleset.update(dict.fromkeys(udef_segs, 'F'))
      return sim_seg(seg, ruleset)
    if 'F' in def_states:
      if (def_states - {'F'}) or (len(cfgs) < len(cfg)):
        return None # Misfired.
      """All defined states are 'F'; set any remaining
      undefined states to 'F' to avoid misfire."""
      ruleset.update(dict.fromkeys(udef_segs, 'F'))
      return sim_seg(seg, ruleset)
    if not udef_segs: return sim_seg(seg, ruleset)

    def update_ruleset(new_rules):
      new_ruleset = dict(ruleset)
      if sym:
        new_rules.update({r[::-1]:new_rules[r]
                          for r in new_rules})
      new_ruleset.update(new_rules)
      defsim(cfg, set(cfgs), new_ruleset, lidx)

    if sym:
      udef_segs = list(set(min(s, s[::-1]) for s in udef_segs))
    if (not def_states) and (len(cfgs) >= len(cfg)):
      update_ruleset(dict.fromkeys(udef_segs, 'F'))
    for s_tup in product(nonF_states, repeat=len(udef_segs)):
      update_ruleset(dict(zip(udef_segs, s_tup)))
    return None

  def defsim(cfg, cfgs, ruleset, lidx):
    n = len(cfg)
    allF = n * 'F'
    while cfg != allF:
      cfg = defsim_cfg(cfg, cfgs, ruleset, lidx)
      if (not cfg) or (cfg in cfgs): return
      cfgs.add(cfg)
    lidx += 1
    if lidx < len(lengths):
      n = lengths[lidx]
      cfg = 'G' + (n-1)*'Q'
      cfgs = set([cfg])
      defsim(cfg, set([cfg]), ruleset, lidx)
      return
    save_ruleset(ruleset)

  n = lengths[0]
  cfg = 'G'+(n-1)*'Q'
  defsim(cfg, set([cfg]), init_ruleset, 0)


def main():
  bgn_tm = clock()
  defsim_ring('power_of_2_sols.txt', {'QQQ':'Q'}, [2,4,8,16,32],
              ['Q','G','A'])
  print('elapsed_time = ' + str(clock() - bgn_tm) + ' secs')


if __name__ == '__main__': main()
