# Security Policy

## Supported Versions

Proteus is pre-1.0. Security fixes are applied to the `main` branch and the
latest tagged release.

## Reporting A Vulnerability

Please report suspected vulnerabilities privately by opening a GitHub security
advisory for this repository. If GitHub advisories are unavailable, contact the
maintainer through the email listed in the GitHub profile for `jvogan`.

Include:

- affected script or workflow
- command and input that triggered the issue
- expected impact
- whether a public API, local file, or external binary was involved

Do not open a public issue for secrets exposure, command injection, unsafe file
handling, or other exploitable behavior until a fix is available.

## Security Model

Proteus helper scripts are local command-line tools. They do not run a server,
collect credentials, or execute remote code by design. Some workflows call
public biology APIs or local visualization binaries such as PyMOL and ChimeraX.
Treat untrusted structure files as untrusted input, keep local tools updated,
and run workflows in a controlled workspace when analyzing files from unknown
sources.
