#! /usr/bin/env python

import re
import sys
from functools import reduce

re_expr = re.compile('.*acc-ratios.*?([0-9,]*) b-read.*?([0-9,]*) b-written.*')

def get_stats(str):
    match = re_expr.match(str)
    if match:
        return [int(number.replace(',', '')) for number in match.group(1, 2)]

def parse(path):
    with open(path) as file:
        return reduce(
            lambda first, second: (first[0] + second[0], first[1] + second[1]),
            filter(lambda obj: obj, map(get_stats, file.readlines()))
        )

if __name__ == '__main__':
    if len(sys.argv) > 1:
        print('bytes read: {stats[0]}\nbites written: {stats[1]}'.format(stats=parse(sys.argv[1])))
    else:
        print("Pass the path to the file with results of DHAT as an argument")
