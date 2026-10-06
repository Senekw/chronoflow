"""Build a local Tectonic bundle directory (texbundle/files) by byte-range download from the default bundle, for
environments where tectonic's own network access fails.  ~1.8 GB (macros, fonts, metrics, hyphenation); ~2 min.
Usage: python fetch_texbundle.py ; then  tectonic -X compile --bundle ../texbundle/files chronoflow.tex"""
import requests, gzip, os, re, time
sess = requests.Session()
def resolve(name):
    r = sess.head("https://relay.fullyjustified.net/" + name, allow_redirects=False, timeout=30)
    return r.headers.get("Location", "https://relay.fullyjustified.net/" + name)
idx_url = resolve("default_bundle_v33.tar.index.gz"); tar_url = resolve("default_bundle_v33.tar")
idx = gzip.decompress(sess.get(idx_url, timeout=60).content).decode()
entries = {}
for ln in idx.strip().split("\n"):
    name, off, ln_ = ln.rsplit(" ", 2); entries[name] = (int(off), int(ln_))
os.makedirs("texbundle/files", exist_ok=True)
open("texbundle/files/SHA256SUM", "w").write("a" * 64)      # tectonic only checks the digest's format for a directory bundle
ext = {"sty","cls","def","cfg","clo","fd","ltx","tex","ini","dat","lua","enc","map","tcx","pool","txt","bst","tec","pat","hyp","lic","chr","ec","dic","tfm","vf","ofm","ovf","ovp","pfb"}
want = [n for n in entries if n != "SHA256SUM" and ((n.rsplit(".", 1)[-1].lower() in ext) or ("." not in n) or ("hyph" in n.lower())
        or (n.lower().endswith((".otf", ".ttf")) and re.match(r"^(lm|latinmodern|texgyre|xits|libertinus|dejavu)", n, re.I)))]
segs = sorted((entries[n][0], entries[n][0] + entries[n][1]) for n in want); ranges = []
for a, b in segs:
    if ranges and a - ranges[-1][1] <= 1e6: ranges[-1][1] = max(ranges[-1][1], b)
    else: ranges.append([a, b])
t0 = time.time()
for i, (a, b) in enumerate(ranges):
    for attempt in range(4):
        try:
            r = sess.get(tar_url, headers={"Range": f"bytes={a}-{b-1}"}, timeout=600); assert r.status_code == 206; buf = r.content; break
        except Exception:
            if attempt == 3: raise
            time.sleep(3)
    for n in want:
        o, l = entries[n]
        if o >= a and o + l <= b:
            open(f"texbundle/files/{n}", "wb").write(buf[o-a:o-a+l])
print(len(os.listdir("texbundle/files")), "files", f"{time.time()-t0:.0f}s")
