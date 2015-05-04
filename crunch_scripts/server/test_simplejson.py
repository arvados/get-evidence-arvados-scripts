#!/usr/bin/python

import os
import sys
import simplejson as json
import re
import gzip

print sys.argv[1]

getev_flatfile = sys.argv[1]

f_in = getev_flatfile
if isinstance(getev_flatfile, str):
  if re.search(r'\.gz$', getev_flatfile):
    f_in = gzip.open(getev_flatfile)
  else:
    f_in = open(getev_flatfile)

count = 0
for line in f_in:
  data = json.loads(line)

  print count
  count += 1

