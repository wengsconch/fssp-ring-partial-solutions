"""Find exponential rulesets."""

from itertools import (product, permutations, combinations)
from math import (ceil, log)
from os import fsync
from time import clock
from ruleset_utils import *


def find_exp_rulesets_for_seg(init_seg, init_ruleset, base,
                              nonQF_states):

  """Find base-<base> exponential rulesets for <init_seg> using
  states from <nonQF_states>, Q and F. The resulting rulesets will
  be consistent with <init_ruleset> and enable simulating
  <init_seg> to one time step.

  Notation and terminology:
      state -- A string comprising one character.
      Q -- Quiescent state.
      F -- Firing state.
      all_states -- <nonQF_states> together with Q and F.
      segment -- A string of states taken from all_states
      ruleset -- A dict mapping segments to states, where the
                 segments are of the same length as those in
                 <init_ruleset>.

  Arguments:
      init_seg -- An initial segment for which to find exponential
                  rulesets.
      init_ruleset -- A dict denoting a non-empty initial ruleset.
                      This ruleset need not be exponential.
      base -- An integer greater than 1.
      nonQF_states -- A list of states that are neither Q nor F.

  Returns:
      A list, possibly empty, of exponential rulesets.

  """
  mseg_len = get_mseg_len(init_ruleset)
  null = mseg_len - 1
  nonF_states = nonQF_states + ['Q']
  found_rulesets = []

  def recurse(ruleset):

    def defsim_seg(seg, nofire=False):
      """Define and simulate <seg>.

      <nofire> -- if True, <seg> should not fire.
                  if False (default), <seg> may fire or not.
      """
      udef_segs, def_states = get_udef_segs_def_states(
          seg, ruleset, null)
      if 'F' in def_states:
        # Misfired.
        if nofire or (def_states - {'F'}): return False
        """All defined states are 'F'; set any remaining
        undefined states to 'F' to avoid misfire."""
        ruleset.update(dict.fromkeys(udef_segs, 'F'))
        return True
      if not udef_segs: return True # All msegs are defined.

      def add_new_rules(new_rules):
        new_ruleset = dict(ruleset)
        new_ruleset.update(new_rules)
        recurse(new_ruleset)

      if not (nofire or def_states):
        # Define <seg> to fire in next time step.
        add_new_rules(dict.fromkeys(udef_segs, 'F'))
      for s_tup in product(nonF_states, repeat=len(udef_segs)):
        # Define <seg> to not fire in the next time step.
        add_new_rules(dict(zip(udef_segs, s_tup)))
      return False

    def defsim_end_seg(seg, s):
      """Define and simulate the penultimate image <seg> of the
      exponentiated beginning segment <bgn_seg> given that
      <ruleset>[<bgn_seg>] = <s>."""
      e = seg[:mseg_len]
      # Conflicting definitions.
      if e in ruleset and ruleset[e] != s: return False
      ruleset[e] = s
      x = 'F' if s == 'F' else 'Q'
      udef_segs, def_states = get_udef_segs_def_states(
          seg[1:], ruleset, null)
      # Non-exponential form.
      if (def_states - {x}): return False
      ruleset.update(dict.fromkeys(udef_segs, x))
      return True

    """Start off the search procedure by defining and simulating
    the initial segment (which is not necessarily minimal)."""
    if not defsim_seg(init_seg): return
    mseg_queue = MsegQueue(null)
    mseg_queue.enqueue(init_seg)

    """While there are queued msegs, dequeue an mseg <seg>,
    exponentiate it to base <base>, and simulate it <base> time
    steps. <seg> must be defined under <ruleset>, say with rule
    <ruleset>[<seg>] = <s>. Backtrack by returning to the last
    call of recurse if a contradiction arises at any point of the
    simulation. If any undefined msegs are encountered, define
    them and deepen the search by calling recurse(). With each
    successfully simulated segment, add its as yet unencountered
    msegs to the queue. If the simulation makes it through to
    the end of the loop body, we end up with
    <ruleset>[<seg>] = <s> being exponential i.e. <ruleset> maps
    exponentiated <seg> to exponentiated <s>, provided that any
    rules used in the above-mentioned simulation are ultimately
    shown to be exponential as well."""
    while not mseg_queue.empty():
      seg = mseg_queue.dequeue()
      assert seg in ruleset
      s = ruleset[seg]
      seg = exp_seg(seg, base)
      for i in range(base - 1):
        if not defsim_seg(seg, nofire=True): return
        mseg_queue.enqueue(seg)
        seg = sim_seg(seg, ruleset)
      if not defsim_end_seg(seg, s): return
      mseg_queue.enqueue(seg)
    found_rulesets.append(dict(ruleset))

    recurse(init_ruleset)
    return found_rulesets


def find_firing_exp_rulesets_for_rings(base, null, init_len,
                                       nonQF_states,
                                       init_ruleset={}):
  """Find firing exponential rulesets for ring arrays.

  Notation and terminology:
      state -- A string comprising one character.
      Q -- Quiescent state.
      F -- Firing state.
      all_states -- <nonQF_states> together with Q and F.
      segment -- A String of states taken from all_states
      ruleset -- A dict mapping segments to states, where the
                 segments are of length <null>+1.

  Arguments:
      base -- An integer greater than 1.
      null -- A positive integer denoting the nullity.
      init_len -- A non-negative integer denoting the initial ring
                  from which to start searching.
      nonQF_states -- A list of states that are neither Q nor F.
                      The first state in the list is the initiator
                      state.
      init_ruleset -- A dict denoting the initial ruleset. The
                      initial ruleset need not be exponential. It
                      is empty by default.

  Returns:
      A generator yielding resulting rulesets.

  """
  init_ruleset[(null + 1) * 'Q'] = 'Q'
  m = int(ceil(null / init_len))
  allF = init_len * 'F'
  initiator = nonQF_states[0]

  cfgs = [initiator + (init_len - 1) * 'Q']
  cfgs_set = {cfgs[-1]}
  def recurse(ruleset):
    unrolled_cfg = cfgs[-1] + (m * cfgs[-1])[:null]
    new_rulesets = find_exp_rulesets_for_seg(
        unrolled_cfg, ruleset, base, nonQF_states)
    for new_ruleset in new_rulesets:
      new_ruleset.update(ruleset)
      cfg = sim_seg(unrolled_cfg, new_ruleset)
      if cfg == allF: # Fired.
        yield new_ruleset
        continue
      if cfg in cfgs_set: continue # Non-terminal.
      cfgs.append(cfg)
      cfgs_set.add(cfg)
      for r in recurse(new_ruleset):
        yield r
      cfgs_set.remove(cfg)
      cfgs.pop()

  return (r for r in recurse(init_ruleset))


def record_firing_exp_rulesets_for_rings(base, null, init_len,
                                         nonQF_states,
                                         init_ruleset={}, id=''):
  """Record firing exponential rulesets for ring arrays.

  Arguments:
      See find_firing_exp_rulesets_for_rings().
      id -- identifier to distinguish runs based on init_ruleset

  Side effects:
      Each resulting ruleset is immediately written to a file. The
      file name is determined by <base>, <null>, <exp> and
      <nonQF_states>.

  Returns:
      None

  """
  init_ruleset[(null + 1) * 'Q'] = 'Q'
  outstr = ('base%d_null%d_init_len%d_%s' %
            (base, null, init_len, ''.join(nonQF_states)))
  outfile = outstr + '_' + str(id) + '.txt'
  f = open(outfile, 'w')
  f.write(outstr + ' = dict(\n')
  f.write('base=%d, null=%d, init_len=%d, nonQF_states=%s,\n' %
          (base, null, init_len, str(nonQF_states)))
  f.write('init_ruleset=%s,\n' % str(init_ruleset))
  f.write('rulesets=[\n')
  f.close()
  print('find_firing_exp_rulesets_for_rings:')
  print('base -- %d' % base)
  print('nullity -- %d' % null)
  print('init_len -- %d' % init_len)
  print('nonQF_states -- ' + str(nonQF_states))
  print('init_ruleset -- ' + str(init_ruleset) + '\n')
  bgn_tm = clock()
  rulesets = find_firing_exp_rulesets_for_rings(
    base, null, init_len, nonQF_states,
    init_ruleset=init_ruleset)
  cnt = 1
  for ruleset in rulesets:
    if (get_ruleset_states(ruleset) - {'Q', 'F'}
        < set(nonQF_states)):
      print('Ruleset used less states than provided.')
      continue
    f = open(outfile, 'a')
    f.write(str(ruleset) + ',\n')
    """Merely closing the file does not guarantee writing the
    data to disk immediately."""
    f.flush()
    fsync(f.fileno())
    f.close()
    print('!!!Found ruleset %s' % (str(cnt)))
    cnt += 1
    f = open(outfile, 'a')
    f.write(']\n)\n')
    f.close()
    print('elapsed_time = ' + str(clock() - bgn_tm) + ' secs')


def main():
  #base = 2; null = 3; init_len = 2; nonQF_states = ['R']
  #base = 2; null = 4; init_len = 2; nonQF_states = ['R']
  #base = 2; null = 5; init_len = 2; nonQF_states = ['R']
  #base = 2; null = 6; init_len = 2; nonQF_states = ['R']
  base = 2; null = 2; init_len = 4; nonQF_states = ['R', 'S']
  #base = 2; null = 3; init_len = 2; nonQF_states = ['R', 'S']
  #base = 2; null = 2; init_len = 2; nonQF_states = ['R', 'S', 'T']

  record_firing_exp_rulesets_for_rings(base, null, init_len,
                                       nonQF_states)


if __name__ == '__main__': main()
