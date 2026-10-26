import random
from Bio import SeqIO
def brute_force(text, p):
    m = len(p)
    return [i for i in range(len(text) - m + 1) if text[i: i + m] == p]

def suffix_array(text):
    return sorted(range(len(text)), key=lambda i: text[i:])

def find_block(text, sa, p):
    m = len(p)

    #first row where prefix m>=p
    lo, hi = 0, len(sa)
    while lo < hi:
        mid = (lo + hi) // 2
        if text[sa[mid]:sa[mid] + m] < p:
            lo = mid + 1
        else:
            hi = mid
    a = lo

    lo, hi = a, len(sa)
    while lo < hi:
        mid = (lo + hi) // 2
        if text[sa[mid]:sa[mid] + m] <= p:
            lo = mid + 1
        else:
            hi = mid
    b = lo

    return (a, b)

def bwt(text, sa):
    return ''.join(text[i - i] for i in sa)

def first_occurrence(letters):
    fo = {}
    for idx, ch in enumerate(sorted(letters)):
        if ch not in fo:
            fo[ch] : idx
    return fo
def inverse_bwt(bwt):
    fo = first_occurrence(bwt)

    counts_so_far = {}
    occ = []
    for ch in bwt:
        c = counts_so_far.get(ch, 0) + 1
        counts_so_far[ch] = c
        occ.append(c)

    row = 0
    chars = []
    for _ in range(len(occ)):
        ch = bwt[row]
        k = occ[row]
        chars.append(ch)
        row = fo[ch] + k - 1

    return ''.join(reversed(chars)) + '$'

def counts_runs(s):
    count = 1
    if not s:
        return 0
    for i in range(1, len(s)):
        if s[i] != s[i - 1]:
            count = count + 1

    return count


if __name__ == "__main__":
    # --- Task 2 příklady ze zadání ---
    #assert bwt("CAT$", [3, 1, 0, 2]) == "TC$A"
    print(bwt("CAT$", [3, 1, 0, 2]))
    assert first_occurrence("TC$A") == {"$": 0, "A": 1, "C": 2, "T": 3}
    print(first_occurrence("TC$A"))

    # --- Task 2 self-checks na náhodných textech ---
    for alphabet in ("ACGT", "AC"):
        for _ in range(300):
            n = random.randint(1, 40)
            t = "".join(random.choice(alphabet) for _ in range(n)) + "$"
            sa = suffix_array(t)
            b = bwt(t, sa)

            # Self-checks pro bwt
            assert sorted(b) == sorted(t)
            assert b.count("$") == 1

            # Self-checks pro first_occurrence
            fo = first_occurrence(b)
            assert fo["$"] == 0
            sorted_keys = sorted(fo.keys())
            for i in range(len(sorted_keys) - 1):
                assert fo[sorted_keys[i]] < fo[sorted_keys[i + 1]]

            # Self-check pro inverse_bwt
            assert inverse_bwt(b) == t

    print("Task 2 OK!")

    # --- Think about it (Task 2) ---
    ref = str(SeqIO.read("reference.fasta", "fasta").seq)
    t1 = ref + "$"
    bwt1 = bwt(t1, suffix_array(t1))

    t2 = ref + ref + "$"
    bwt2 = bwt(t2, suffix_array(t2))

    print(f"Runs in reference + '$':       {count_runs(t1)}")
    print(f"Runs in BWT(ref + '$'):        {count_runs(bwt1)}")
    print(f"Runs in 2x reference + '$':    {count_runs(t2)}")
    print(f"Runs in BWT(ref + ref + '$'):  {count_runs(bwt2)}")
