"""Common utilities for rulesets."""
 
from sys import version_info
 
 
def enforce_python_version(version):
    if version_info.major < version:
        raise BaseException('Must use Python version ' +
                            str(version) + ' or greater.')
enforce_python_version(3)
 
 
def get_mseg_len(ruleset):
    """Get the minimal segment length from non-empty <ruleset>."""
    return len(list(ruleset.keys())[0])
 
 
def get_msegs(seg, null, gran=1):
    """Get minimal segments (length <null> + 1) of <seg> with respect
    to granularity <gran> (default 1)."""
    assert not (len(seg)%gran or null%gran)
    mseg_len = null + gran
    return set([seg[i:i + mseg_len]
                for i in range(0, len(seg) - null, gran)])
 
 
def get_udef_segs_def_states(seg, ruleset, null=2):
    """Get undefined minimal segments and defined states for <seg>."""
    msegs = get_msegs(seg, null)
    udef_segs = msegs - ruleset.keys()
    def_states = {ruleset[t] for t in (msegs - udef_segs)}
    return list(udef_segs), def_states
 
 
class MsegQueue:
    """Queue for scheduling minimal segments to be exponentiated."""
    def __init__(self, null, gran=1):
        self.null = null
        self.gran = gran
        self.queued_msegs = []
        self.seen_msegs = set()
 
    def empty(self):
        return not self.queued_msegs
 
    def enqueue(self, seg):
        """Break <seg> up into minimal segments and enqueue previously
        unencountered ones."""
        new_msegs = (get_msegs(seg, self.null, self.gran) -
                     self.seen_msegs)
        self.queued_msegs.extend(list(new_msegs))
        self.seen_msegs.update(new_msegs)
 
    def dequeue(self):
        return self.queued_msegs.pop(0)
 
    def seen_msegs_iter(self):
        return self.seen_msegs.__iter__()
 
 
def sim_seg(seg, ruleset, tskip=1):
    """Simulate <seg> for <tskip> (default 1) time-steps using
    <ruleset>."""
    mseg_len = get_mseg_len(ruleset)
    try:
        for t in range(tskip):
            seg = ''.join(ruleset[seg[i:i + mseg_len]]
                          for i in range(len(seg) - mseg_len + 1))
    except KeyError:  # Unmapped minimal preimage segments
        print(len(seg))
        seg = None
    return seg
 
def split_seg(seg, gran=1):
    """Return generator expression for <gran>-sized (default 1) chunks
    of <seg>."""
    return (seg[i:i + gran] for i in range(0, len(seg), gran))
 
 
def exp_seg(seg, base, gran=1):
    """Exponentiate <seg> to base <base> with respect to granularity
    <gran> (default 1).
 
    If non-empty <seg> is firing then it exponentiates with compound
    state <gran>*F, else with <gran>*Q.
    Example 1 -- FFF base-2 exponentiates to FFFFFF for <gran>=1.
    Example 2 -- RST base-3 exponentiates to RQQSQQTQQ for <gran>=1.
    Example 3 -- FFFF base-2 exponentiates to FFFFFFFF for <gran>=2.
    Example 4 -- RSTU base-2 exponentiates to RSQQTUQQ for <gran>=2.
 
    """
    gF, gQ = [gran * s for s in ('F', 'Q')]
    bgX = (base - 1) * (gF if seg[:gran] == gF else gQ)
    return ''.join((s + bgX) for s in split_seg(seg, gran))
 
 
def inv_exp_seg(seg, base, gran=1):
    return ''.join(seg[i:i+gran] for i in range(0, len(seg), gran*base))
 
 
def repr_ruleset(ruleset):
    """Return unambiguous string representation of <ruleset>."""
    return '|'.join((t + ':' + ruleset[t])
                    for t in sorted(ruleset.keys()))
 
 
def print_ruleset(ruleset):
    dist = [[] for i in range(get_mseg_len(ruleset) + 1)]
    [dist[mseg.count('Q')].append(mseg) for mseg in ruleset]
    for i,segs in enumerate(dist):
        print('[' + str(i) + '] ', end='')
        print(','.join(s + ':' + ruleset[s] for s in sorted(segs)))
 
 
def get_ruleset_states(ruleset):
    """Get states used in <ruleset>."""
    return {s for k,v in ruleset.items() for s in (k + v)}
