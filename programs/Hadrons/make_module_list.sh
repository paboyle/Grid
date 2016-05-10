#!/usr/bin/env bash

echo 'modules =\' > modules.inc
find Modules -name '*.cc' -type f -print | sed 's/^/  /;$q;s/$/ \\/' >> modules.inc