from ruleset_utils import *

def find_exp_from_base_ruleset(tsegs, base_ruleset, stride, gran,
                               rename_states, base):
  """Find base-<base> exponential ruleset embedded in
  <base_ruleset> originating from initial segments <tsegs>.
  Aggregate states are of granularity <gran> in terms of the base
  states of <base_ruleset> and are renamed according to states
  supplied in <rename_states>. Exponential ruleset is <stride>-
  power of <base_ruleset>.

  Notation and terminology:
      state -- A string comprising one character.
      Q -- Quiescent state.
      F -- Firing state.
      segment -- A string of states.
      ruleset -- A dict mapping segments to states, where the
                 segments are of the same length as those in
                 <base_ruleset>.

  Arguments:
      base_ruleset -- A dict denoting a non-empty ruleset.
      stride -- A positive integer denoting the power to which
                <base_ruleset> is raised.
      gran -- A positive integer denoting the granularity of
              aggregate states.
      tsegs -- A list of initial segments from which to start
               finding an exponential ruleset.
      rename_states -- A list of states used to rename aggregate
                       states encountered.
      base -- An integer greater than 1.

  Returns:
      An exponential ruleset if one exists else None.

  """
  null = (get_mseg_len(base_ruleset) - 1) * stride
  assert not (null % gran)
  rename_states_gen = (s for s in rename_states)
  rename_map = {gran * 'Q': 'Q', gran * 'F': 'F'}
  exp_ruleset = {}

  def rename_seg(seg):
    renamed_seg = []
    for s in split_seg(seg, gran):
      if s not in rename_map:
        rename_map[s] = next(rename_states_gen)
      renamed_seg.append(rename_map[s])
    return ''.join(renamed_seg)

  mseg_queue = MsegQueue(null, gran)
  [mseg_queue.enqueue(s) for s in tsegs]
  while not mseg_queue.empty():
    bgn_seg = mseg_queue.dequeue()
    seg = exp_seg(bgn_seg, base, gran)
    mseg_queue.enqueue(seg)
    for i in range(base):
      seg = sim_seg(seg, base_ruleset, stride)
      if not seg:
        return None
      mseg_queue.enqueue(seg)
    end_seg = sim_seg(bgn_seg, base_ruleset, stride)
    if not (end_seg and seg == exp_seg(end_seg, base, gran)):
      return None
    exp_ruleset[rename_seg(bgn_seg)] = rename_seg(end_seg)
    return exp_ruleset, rename_map
