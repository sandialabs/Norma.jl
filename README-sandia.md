# Running Norma.jl on the Sandia Network (macOS & Linux)

Unofficial setup notes for getting Norma.jl running on a Sandia-managed machine.

Norma's standard install (see `README.md`) is just:

```bash
git clone https://github.com/sandialabs/Norma.jl.git
cd Norma.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. src/Norma.jl examples/<case>/input.yaml
```

That works **everywhere except** the two points where Sandia's TLS inspection
(ZScaler on laptops / the SRN proxy) breaks tools that ship their **own** CA
store. macOS- and Linux-native tools trust Sandia's interception via the system
trust store, but **Julia and conda do not** — they carry their own CA bundles
and will fail with `certificate verify failed` until pointed at a bundle that
includes the **Sandia Root CA**.

So the only Sandia-specific work is **Section A** below (set up CA trust once).
Sections B–D are the normal Julia workflow.

---

## A. One-time: make Julia and conda trust Sandia's CA

### A.1 — Get a CA bundle that includes the Sandia root

The goal is a single PEM file (we'll use `~/ssl/my-ca-bundle.pem`) containing the
public root CAs **plus** Sandia's CA.

**macOS** — export it straight from the system keychain (Sandia's root is
installed there by MDM):

```bash
mkdir -p ~/ssl
{ security find-certificate -a -p /System/Library/Keychains/SystemRootCertificates.keychain
  security find-certificate -a -p /Library/Keychains/System.keychain
} > ~/ssl/my-ca-bundle.pem
```

**Linux** — a managed Sandia box already has the Sandia CA in the system trust
store, so copy that bundle:

```bash
mkdir -p ~/ssl
# Debian / Ubuntu:
cp /etc/ssl/certs/ca-certificates.crt ~/ssl/my-ca-bundle.pem
# RHEL / Fedora / CentOS:
# cp /etc/pki/tls/certs/ca-bundle.crt ~/ssl/my-ca-bundle.pem
```

If the Sandia root is **not** already in your Linux trust store, get the Sandia
root CA `.pem` from IT (or copy `~/ssl/my-ca-bundle.pem` from a working Mac) and
install it system-wide:

```bash
# Debian / Ubuntu:
sudo cp sandia-root.crt /usr/local/share/ca-certificates/ && sudo update-ca-certificates
# RHEL / Fedora:
sudo cp sandia-root.pem /etc/pki/ca-trust/source/anchors/ && sudo update-ca-trust
```

Verify the bundle works:

```bash
curl --cacert ~/ssl/my-ca-bundle.pem -sSI https://pkg.julialang.org | head -1   # expect HTTP 200/301
```

### A.2 — Point Julia (and other self-bundled tools) at it

Add to your shell startup (`~/.bashrc`, `~/.zshrc`, …):

```bash
export SSL_CERT_FILE="$HOME/ssl/my-ca-bundle.pem"      # Julia Pkg, curl, OpenSSL
export REQUESTS_CA_BUNDLE="$HOME/ssl/my-ca-bundle.pem" # Python / pip
export NODE_EXTRA_CA_CERTS="$HOME/ssl/my-ca-bundle.pem" # Node / npm (optional)
```

Open a new shell (or `source` the file) so these take effect.

### A.3 — Point conda at it (unblocks the PyCall build)

Norma depends on **PyCall**, whose build step runs `conda install numpy` using a
conda it installs into `~/.julia/conda`. That conda's Python ignores
`SSL_CERT_FILE`, so configure it explicitly in `~/.condarc`:

```yaml
ssl_verify: /home/<you>/ssl/my-ca-bundle.pem   # macOS: /Users/<you>/ssl/my-ca-bundle.pem
```

(Use your real absolute path — `~` is not expanded in `.condarc`.)

### A.4 — If you reach the network through `proxy.sandia.gov` (often on VPN)

When you are **not** on the SRN-wired / ZScaler path, you may need the proxy:

```bash
export http_proxy=http://proxy.sandia.gov:80
export https_proxy=http://proxy.sandia.gov:80
export no_proxy=localhost,127.0.0.1,.sandia.gov
```

---

## B. Install Julia (juliaup)

Same on macOS and Linux:

```bash
curl -fsSL https://install.julialang.org | sh
```

This installs the `juliaup` version manager and a `julia` launcher under
`~/.juliaup/bin` (make sure it's on your `PATH`). Manage versions with
`juliaup add <channel>`, `juliaup default <channel>`, `juliaup status`.

---

## C. Get Norma.jl

Use **HTTPS** — it's the least-friction option on the corporate network. (SSH
`git@github.com:…` needs a GitHub key *and* outbound port 22, which ZScaler may
block; if you prefer SSH, configure it over port 443.)

```bash
git clone https://github.com/sandialabs/Norma.jl.git ~/Repos/Norma.jl
cd ~/Repos/Norma.jl
```

---

## D. Instantiate and run

Use **`instantiate`**, not `update` — it reproduces the exact dependency
versions pinned in the committed `Manifest.toml` (faster and avoids surprise
version bumps):

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

This downloads the dependencies and **builds PyCall** (the `conda install numpy`
step that needs Section A.3). Then run a model:

```bash
julia --project=. src/Norma.jl examples/<case>/input.yaml
```

Validate the install with the test suite:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

---

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `Pkg.update()`/first op hangs forever | General registry bootstrap stalled on a transient network hiccup | Re-run; or `julia -e 'using Pkg; Pkg.Registry.add("General")'` then retry. Not a config problem. |
| `Pkg` fails fetching registry/packages with a TLS/cert error | Julia's bundled CA doesn't trust Sandia's intercepted cert | Section A.1–A.2 (`SSL_CERT_FILE` → CA bundle) |
| Build error in **PyCall** → `CondaSSLError: CERTIFICATE_VERIFY_FAILED ... unable to get local issuer certificate` | conda's Python doesn't trust Sandia's cert | Section A.3 (`~/.condarc` `ssl_verify`), then `julia --project=. -e 'using Pkg; Pkg.build("PyCall")'` |
| Connections time out off-site | Need the Sandia proxy | Section A.4 |
| `git clone` over SSH hangs | Outbound port 22 blocked | Use the HTTPS URL (Section C) |

### Quick sanity checks

```bash
echo "$SSL_CERT_FILE"                                   # should be set
curl --cacert "$SSL_CERT_FILE" -sSI https://pkg.julialang.org | head -1
julia -e 'using Downloads; Downloads.download("https://pkg.julialang.org/registries", tempname()); println("pkg server OK")'
```

---

*These are community setup notes, not official Sandia IT guidance. The only
Sandia-specific requirement is CA trust (Section A); everything else is the
standard juliaup → clone → instantiate → run workflow.*
