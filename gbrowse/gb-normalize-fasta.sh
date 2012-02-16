#! /bin/bash

perl -wlpaF'\t' -E 's/^>(.*)$/>\L$1\E/' $1
