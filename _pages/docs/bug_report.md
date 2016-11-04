---
layout: single
title : "Documentation"
author_profile: false
excerpt: "Reporting a bug"
header:
  overlay_color: "#C70039"
permalink: /docs/bug_report/
sidebar:
  nav : docs
---


{% octicon alert height:32 class:"right left" aria-label:hi %} 

__To help us tracking and solving more efficiently issues with Grid, please report problems using the [issue system of GitHub](https://github.com/paboyle/Grid/issues) rather than sending emails to Grid developers.__

We also suggest to have a brief look at the [closed issues pages](https://github.com/paboyle/Grid/issues?q=is%3Aissue+is%3Aclosed) and check whether the problem has been addressed already.

{% capture notice-text %}
1. Check that the code is pointing to the `HEAD` of `develop` or any commit in `master` which is tagged with a version number. 
2. Give a description of the target platform (CPU, network, compiler).
3. Give the exact `configure` command used.
4. Attach `config.log`.
5. Attach `config.summary`.
6. Attach the output of `make V=1`.
7. Describe the issue and any previous attempt to solve it. If relevant, show how to reproduce the issue using a minimal working example.
{% endcapture %}

<div class="notice--warning">
  <h4>When you file an issue, please go though the following checklist:</h4>
  {{ notice-text | markdownify }}
</div>

